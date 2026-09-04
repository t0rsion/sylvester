//! Python bindings for the `sylvester` crate: reduced Gröbner bases under
//! grevlex, over a prime field or over the rational numbers.
//!
//! The Rust library carries the coefficient domain as a type parameter.
//! Python learns the domain at run time, so each pair of Rust
//! instantiations collapses to one Python class holding an enum over the
//! two, and the ring constructor picks the arm. `docs/rational-design.md` section 9 fixes the
//! surface, the error mapping, and the rule around the GIL.
//!
//! Every computation releases the GIL. The rule is one order: convert and
//! clone every input into an owned Rust value, release the GIL, run, take
//! the GIL back, and only then build Python objects. No borrowed Python
//! value crosses into the released region.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::num::NonZeroUsize;
use std::time::Duration;

use num_bigint::BigInt;
use num_rational::BigRational;
use pyo3::exceptions::{PyException, PyIndexError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::sync::GILOnceCell;
use pyo3::types::{PyBytes, PyDict, PyIterator, PyList, PySlice, PyTuple, PyType};

use sylvester::verify::{self, VerifyError};
use sylvester::{
    Backend, BasisError, Budget, CertifiedGroebnerBasis, CertifyError, Coefficient, ComputeError,
    ComputeOptions, ComputeReport, Established, F4Counters, GroebnerBasis, HilbertError,
    HilbertSeries, Ideal, ModularLift, NormalFormError, ParseError, Polynomial, PolynomialRing,
    PrimeField, RationalOptions, RationalStop, Rationals, RingError,
};

/// The largest exponent one variable holds, which is what the engines pack.
const MAX_EXPONENT: i64 = u16::MAX as i64;

// Rust 1.88 marks this const-evaluated helper dead even though the block
// below calls it.
#[allow(dead_code)]
const fn assert_send_sync<T: Send + Sync>() {}

// Releasing the GIL around a computation is sound only while every value
// the closure holds crosses threads. A type that stops being `Send + Sync`
// fails the build here rather than at the call site.
const _: () = {
    assert_send_sync::<PolynomialRing<PrimeField>>();
    assert_send_sync::<PolynomialRing<Rationals>>();
    assert_send_sync::<Polynomial<PrimeField>>();
    assert_send_sync::<Polynomial<Rationals>>();
    assert_send_sync::<Ideal<PrimeField>>();
    assert_send_sync::<Ideal<Rationals>>();
    assert_send_sync::<GroebnerBasis<PrimeField>>();
    assert_send_sync::<GroebnerBasis<Rationals>>();
    assert_send_sync::<CertifiedGroebnerBasis>();
};

/// The exception classes of the package, built once at module import.
///
/// Each class carries two bases: [`SylvesterError`], which catches
/// everything this package raises, and the built-in class the design maps
/// the error to. Bad input is a `ValueError` and an exhausted budget or a
/// defect is a `RuntimeError`, so an ordinary `except ValueError` still
/// works.
struct Exceptions {
    base: Py<PyType>,
    ring: Py<PyType>,
    parse: Py<PyType>,
    basis: Py<PyType>,
    certificate_invalid: Py<PyType>,
    budget_exhausted: Py<PyType>,
    timeout: Py<PyType>,
    memory_limit: Py<PyType>,
    limit_exceeded: Py<PyType>,
    internal_defect: Py<PyType>,
}

static EXCEPTIONS: GILOnceCell<Exceptions> = GILOnceCell::new();

/// Build one exception class.
///
/// Python's own `type` builds it, because a class with two bases has no
/// pyo3 macro form.
fn new_exception<'py>(
    py: Python<'py>,
    name: &str,
    bases: &[&Bound<'py, PyType>],
    doc: &str,
) -> PyResult<Py<PyType>> {
    let namespace = PyDict::new(py);
    namespace.set_item("__doc__", doc)?;
    namespace.set_item("__module__", "sylvester")?;
    let bases = PyTuple::new(py, bases)?;
    let class = py.get_type::<PyType>().call1((name, bases, namespace))?;
    Ok(class.downcast_into::<PyType>()?.unbind())
}

struct LeafExceptions {
    ring: Py<PyType>,
    parse: Py<PyType>,
    basis: Py<PyType>,
    certificate_invalid: Py<PyType>,
    timeout: Py<PyType>,
    memory_limit: Py<PyType>,
    limit_exceeded: Py<PyType>,
    internal_defect: Py<PyType>,
}

fn build_exceptions(py: Python<'_>) -> PyResult<Exceptions> {
    let exception = py.get_type::<PyException>();
    let runtime_error = py.get_type::<PyRuntimeError>();
    let base = new_exception(
        py,
        "SylvesterError",
        &[&exception],
        "The base of every exception this package raises. Its subclasses \
         also inherit ValueError for bad input and RuntimeError for an \
         exhausted budget, a structural limit, or a defect.",
    )?;
    let base_type = base.bind(py).clone();
    let budget_exhausted = new_exception(
        py,
        "BudgetExhausted",
        &[&base_type, &runtime_error],
        "A computation stopped on its own budget. Nothing about the ideal \
         follows. Timeout and MemoryLimitExceeded are the two subclasses.",
    )?;
    let budget_type = budget_exhausted.bind(py).clone();
    let leaves = build_leaf_exceptions(py, &base_type, &budget_type)?;
    Ok(Exceptions {
        ring: leaves.ring,
        parse: leaves.parse,
        basis: leaves.basis,
        certificate_invalid: leaves.certificate_invalid,
        timeout: leaves.timeout,
        memory_limit: leaves.memory_limit,
        limit_exceeded: leaves.limit_exceeded,
        internal_defect: leaves.internal_defect,
        base,
        budget_exhausted,
    })
}

fn build_leaf_exceptions(
    py: Python<'_>,
    base: &Bound<'_, PyType>,
    budget: &Bound<'_, PyType>,
) -> PyResult<LeafExceptions> {
    let value_error = py.get_type::<PyValueError>();
    let runtime_error = py.get_type::<PyRuntimeError>();
    Ok(LeafExceptions {
        ring: new_exception(
            py,
            "RingError",
            &[base, &value_error],
            "A ring construction failed. The value does not meet the ring contract.",
        )?,
        parse: new_exception(
            py,
            "ParseError",
            &[base, &value_error],
            "The text is not a polynomial of the ring. The message gives the byte position.",
        )?,
        basis: new_exception(
            py,
            "BasisError",
            &[base, &value_error],
            "A supplied list is not a reduced Gröbner basis.",
        )?,
        certificate_invalid: new_exception(
            py,
            "CertificateInvalid",
            &[base, &value_error],
            "The verifier rejected the certificate bytes.",
        )?,
        timeout: new_exception(py, "Timeout", &[budget], "The deadline passed.")?,
        memory_limit: new_exception(
            py,
            "MemoryLimitExceeded",
            &[budget],
            "The live data passed the memory limit.",
        )?,
        limit_exceeded: new_exception(
            py,
            "LimitExceeded",
            &[base, &runtime_error],
            "The computation reached a structural limit.",
        )?,
        internal_defect: new_exception(
            py,
            "InternalDefect",
            &[base, &runtime_error],
            "Sylvester contradicted an internal invariant.",
        )?,
    })
}

fn exceptions(py: Python<'_>) -> &Exceptions {
    EXCEPTIONS
        .get(py)
        .expect("the module import built the exception classes")
}

/// The error as an instance of `class`, or the failure raised while
/// building it.
fn raise(py: Python<'_>, class: &Py<PyType>, message: String) -> PyErr {
    match class.bind(py).call1((message,)) {
        Ok(value) => PyErr::from_value(value),
        Err(failure) => failure,
    }
}

/// The error with its payload attached, so the variant survives the
/// boundary as attributes and not as a message alone.
fn raise_with(
    py: Python<'_>,
    class: &Py<PyType>,
    message: String,
    attributes: impl FnOnce(&Bound<'_, PyAny>) -> PyResult<()>,
) -> PyErr {
    match class.bind(py).call1((message,)) {
        Ok(value) => match attributes(&value) {
            Ok(()) => PyErr::from_value(value),
            Err(failure) => failure,
        },
        Err(failure) => failure,
    }
}

fn of_ring(py: Python<'_>, error: RingError) -> PyErr {
    match error {
        RingError::Parse(inner) => of_parse(py, inner),
        _ => raise(py, &exceptions(py).ring, error.to_string()),
    }
}

fn of_parse(py: Python<'_>, error: ParseError) -> PyErr {
    raise(py, &exceptions(py).parse, error.to_string())
}

fn of_compute(py: Python<'_>, error: ComputeError) -> PyErr {
    let classes = exceptions(py);
    let message = error.to_string();
    match error {
        ComputeError::Timeout => raise(py, &classes.timeout, message),
        ComputeError::MemoryLimitExceeded => raise(py, &classes.memory_limit, message),
        ComputeError::DegreeLimit { limit } | ComputeError::ExponentLimit { limit } => {
            raise_with(py, &classes.limit_exceeded, message, |value| {
                value.setattr("limit", limit)
            })
        }
        ComputeError::TableFull | ComputeError::PrimesExhausted => {
            raise_with(py, &classes.limit_exceeded, message, |value| {
                value.setattr("limit", py.None())
            })
        }
    }
}

fn of_normal_form(py: Python<'_>, error: NormalFormError) -> PyErr {
    let classes = exceptions(py);
    let message = error.to_string();
    match error {
        NormalFormError::RingMismatch => raise(py, &classes.ring, message),
        NormalFormError::ExponentLimit { limit } => {
            raise_with(py, &classes.limit_exceeded, message, |value| {
                value.setattr("limit", limit)
            })
        }
        NormalFormError::Timeout => raise(py, &classes.timeout, message),
        NormalFormError::MemoryLimitExceeded => raise(py, &classes.memory_limit, message),
    }
}

/// Raise `BasisError` with all three attributes set.
///
/// `index` names the one element that failed a shape check, and `pair`
/// the two whose S-polynomial did not reduce to zero. The attributes the
/// failure does not name are `None`, so reading one never raises
/// `AttributeError`.
fn raise_basis(
    py: Python<'_>,
    message: String,
    index: Option<usize>,
    pair: Option<(usize, usize)>,
) -> PyErr {
    raise_with(py, &exceptions(py).basis, message, |value| {
        value.setattr("index", index)?;
        value.setattr("left", pair.map(|(left, _)| left))?;
        value.setattr("right", pair.map(|(_, right)| right))
    })
}

fn of_basis(py: Python<'_>, error: BasisError) -> PyErr {
    let classes = exceptions(py);
    let message = error.to_string();
    match error {
        BasisError::RingMismatch => raise_basis(py, message, None, None),
        BasisError::ZeroPolynomial { index }
        | BasisError::NotMonic { index }
        | BasisError::NotSorted { index }
        | BasisError::NotInterreduced { index } => raise_basis(py, message, Some(index), None),
        BasisError::NotGroebner { left, right } => {
            raise_basis(py, message, None, Some((left, right)))
        }
        BasisError::ExponentLimit { limit } => {
            raise_with(py, &classes.limit_exceeded, message, |value| {
                value.setattr("limit", limit)
            })
        }
        BasisError::Timeout => raise(py, &classes.timeout, message),
        BasisError::MemoryLimitExceeded => raise(py, &classes.memory_limit, message),
    }
}

fn of_hilbert(py: Python<'_>, error: HilbertError) -> PyErr {
    let classes = exceptions(py);
    let message = error.to_string();
    match error {
        HilbertError::Timeout => raise(py, &classes.timeout, message),
        HilbertError::MemoryLimitExceeded => raise(py, &classes.memory_limit, message),
    }
}

fn of_certify(py: Python<'_>, error: CertifyError) -> PyErr {
    let classes = exceptions(py);
    let message = error.to_string();
    match error {
        CertifyError::Engine(inner) | CertifyError::WriterExhausted(inner) => of_compute(py, inner),
        CertifyError::VerifierExhausted(inner) => of_exhausted_verifier(py, inner),
        CertifyError::CapExceeded { .. } => {
            raise_with(py, &classes.limit_exceeded, message, |value| {
                value.setattr("limit", py.None())
            })
        }
        CertifyError::Emitter(_) | CertifyError::InputMismatch | CertifyError::Rejected(_) => {
            raise(py, &classes.internal_defect, message)
        }
    }
}

/// A verifier that ran out of its own budget, which says nothing about the
/// certificate.
fn of_exhausted_verifier(py: Python<'_>, error: VerifyError) -> PyErr {
    let classes = exceptions(py);
    let message = error.to_string();
    match error {
        VerifyError::DeadlineExceeded => raise(py, &classes.timeout, message),
        VerifyError::CapExceeded { limit, .. } => {
            raise_with(py, &classes.limit_exceeded, message, |value| {
                value.setattr("limit", limit)
            })
        }
        _ => raise(py, &classes.internal_defect, message),
    }
}

/// A verdict on untrusted bytes. A rejection is bad input; an exhausted cap
/// or deadline is not.
fn of_verify(py: Python<'_>, error: VerifyError) -> PyErr {
    if error.is_exhaustion() {
        return of_exhausted_verifier(py, error);
    }
    raise(py, &exceptions(py).certificate_invalid, error.to_string())
}

/// A polynomial of another ring, reported where the domains differ and the
/// library call cannot be made at all.
fn foreign_polynomial(py: Python<'_>, index: usize) -> PyErr {
    raise(
        py,
        &exceptions(py).ring,
        format!("polynomial {index} belongs to another ring"),
    )
}

/// The one polynomial an operation takes belongs to another ring.
fn foreign_argument(py: Python<'_>) -> PyErr {
    raise(
        py,
        &exceptions(py).ring,
        "the polynomial belongs to another ring".to_string(),
    )
}

/// The seconds as a duration.
///
/// A negative, an infinite, a NaN, and a value past what a duration holds
/// are rejected here rather than panicking inside the conversion.
fn duration_of(name: &str, seconds: f64) -> PyResult<Duration> {
    Duration::try_from_secs_f64(seconds).map_err(|_| {
        PyValueError::new_err(format!(
            "{name} is a count of seconds that is finite, not negative, and fits a duration"
        ))
    })
}

fn budget_of(timeout: Option<f64>, memory_limit: Option<usize>) -> PyResult<Budget> {
    let mut budget = Budget::new();
    if let Some(seconds) = timeout {
        budget = budget.timeout(duration_of("timeout", seconds)?);
    }
    if let Some(bytes) = memory_limit {
        budget = budget.memory_limit(bytes);
    }
    Ok(budget)
}

fn backend_of(name: &str) -> PyResult<Backend> {
    match name {
        "f4" => Ok(Backend::F4),
        "classic" => Ok(Backend::Classic),
        _ => Err(PyValueError::new_err(
            "backend is \"f4\" or \"classic\"".to_string(),
        )),
    }
}

fn stop_of(stop: Option<&str>, extra_primes: Option<usize>) -> PyResult<RationalStop> {
    let extra = match extra_primes {
        Some(count) => NonZeroUsize::new(count).ok_or_else(|| {
            PyValueError::new_err("extra_primes is at least 1: zero primes confirm nothing")
        })?,
        None => match RationalStop::default() {
            RationalStop::Unchanged { extra } | RationalStop::ContainsInput { extra } => extra,
        },
    };
    match stop {
        Some("unchanged") => Ok(RationalStop::Unchanged { extra }),
        None | Some("contains_input") => Ok(RationalStop::ContainsInput { extra }),
        Some(_) => Err(PyValueError::new_err(
            "stop is \"unchanged\" or \"contains_input\"".to_string(),
        )),
    }
}

/// The compute options of a prime-field call.
///
/// `stop` and `extra_primes` describe the multimodular engine, which runs
/// over the rational numbers alone. Passing either here is a `ValueError`
/// and not a silent no-op.
fn compute_options(
    backend: Option<&str>,
    timeout: Option<f64>,
    memory_limit: Option<usize>,
    threads: Option<usize>,
    stop: Option<&str>,
    extra_primes: Option<usize>,
) -> PyResult<ComputeOptions> {
    for (name, given) in [
        ("stop", stop.is_some()),
        ("extra_primes", extra_primes.is_some()),
    ] {
        if given {
            return Err(PyValueError::new_err(format!(
                "{name} describes the rational engine, and the ring is a prime field"
            )));
        }
    }
    let mut options = ComputeOptions::new().budget(budget_of(timeout, memory_limit)?);
    if let Some(name) = backend {
        options = options.backend(backend_of(name)?);
    }
    if let Some(count) = threads {
        options = options.threads(count);
    }
    Ok(options)
}

fn rational_options(
    backend: Option<&str>,
    timeout: Option<f64>,
    memory_limit: Option<usize>,
    threads: Option<usize>,
    stop: Option<&str>,
    extra_primes: Option<usize>,
) -> PyResult<RationalOptions> {
    let mut options = ComputeOptions::new().budget(budget_of(timeout, memory_limit)?);
    if let Some(name) = backend {
        options = options.backend(backend_of(name)?);
    }
    if let Some(count) = threads {
        options = options.threads(count);
    }
    Ok(RationalOptions::new()
        .compute(options)
        .stop(stop_of(stop, extra_primes)?))
}

/// A coefficient a caller passes to `PolynomialRing.polynomial`.
fn coefficient_of(value: &Bound<'_, PyAny>) -> PyResult<Coefficient> {
    if let Ok((numerator, denominator)) = value.extract::<(BigInt, BigInt)>() {
        return Ok(Coefficient::Fraction {
            numerator,
            denominator,
        });
    }
    if let Ok(integer) = value.extract::<BigInt>() {
        return Ok(Coefficient::Integer(integer));
    }
    if let Ok(rational) = value.extract::<BigRational>() {
        return Ok(Coefficient::Rational(rational));
    }
    Err(PyValueError::new_err(
        "a coefficient is an int, a fractions.Fraction, or a (numerator, denominator) pair"
            .to_string(),
    ))
}

fn exponents_of(exponents: &[i64]) -> PyResult<Vec<u16>> {
    exponents
        .iter()
        .map(|&exponent| {
            if !(0..=MAX_EXPONENT).contains(&exponent) {
                return Err(PyValueError::new_err(format!(
                    "an exponent is between 0 and {MAX_EXPONENT}"
                )));
            }
            Ok(exponent as u16)
        })
        .collect()
}

/// A ring over one of the two domains.
#[derive(Clone, Debug, PartialEq, Eq)]
enum AnyRing {
    Prime(PolynomialRing<PrimeField>),
    Rational(PolynomialRing<Rationals>),
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum AnyPolynomial {
    Prime(Polynomial<PrimeField>),
    Rational(Polynomial<Rationals>),
}

#[derive(Clone, Debug)]
enum AnyIdeal {
    Prime(Ideal<PrimeField>),
    Rational(Ideal<Rationals>),
}

#[derive(Clone, Debug)]
enum AnyBasis {
    Prime(GroebnerBasis<PrimeField>),
    Rational(GroebnerBasis<Rationals>),
}

impl AnyBasis {
    fn len(&self) -> usize {
        match self {
            AnyBasis::Prime(basis) => basis.len(),
            AnyBasis::Rational(basis) => basis.len(),
        }
    }

    fn get(&self, index: usize) -> AnyPolynomial {
        match self {
            AnyBasis::Prime(basis) => AnyPolynomial::Prime(basis[index].clone()),
            AnyBasis::Rational(basis) => AnyPolynomial::Rational(basis[index].clone()),
        }
    }
}

fn hash_of(polynomial: &AnyPolynomial) -> isize {
    let mut hasher = DefaultHasher::new();
    match polynomial {
        AnyPolynomial::Prime(f) => {
            f.ring().hash(&mut hasher);
            for (coeff, exps) in f.terms() {
                coeff.hash(&mut hasher);
                exps.hash(&mut hasher);
            }
        }
        AnyPolynomial::Rational(f) => {
            f.ring().hash(&mut hasher);
            for (coeff, exps) in f.terms() {
                coeff.hash(&mut hasher);
                exps.hash(&mut hasher);
            }
        }
    }
    hasher.finish() as isize
}

/// A polynomial ring under the grevlex order, over `F_p` or over `Q`.
///
/// Build one with `PolynomialRing.prime_field(p, variables)` or
/// `PolynomialRing.rationals(variables)`. The ring fixes the domain, the
/// variable names, and the variable order, and every polynomial, ideal, and
/// basis comes from it. Two rings are equal, and hash equally, when they
/// hold the same domain and the same variables in the same order.
#[pyclass(frozen, name = "PolynomialRing", module = "sylvester")]
struct PyRing {
    inner: AnyRing,
}

#[pymethods]
impl PyRing {
    /// The ring `F_p[variables]`.
    ///
    /// `p` is a prime of at most 2^31 - 1, and the ring takes at most 256
    /// variables. Each name starts with an ASCII letter or an underscore
    /// and holds only ASCII letters, digits, and underscores. Names must be
    /// distinct, and their order is the variable order. Raises RingError
    /// otherwise.
    #[staticmethod]
    #[pyo3(text_signature = "(p, variables)")]
    fn prime_field(py: Python<'_>, p: u64, variables: Vec<String>) -> PyResult<Self> {
        let ring = PolynomialRing::prime_field(p, variables).map_err(|e| of_ring(py, e))?;
        Ok(PyRing {
            inner: AnyRing::Prime(ring),
        })
    }

    /// The ring `Q[variables]`.
    ///
    /// The name rules and the variable cap are the ones `prime_field`
    /// applies. A rational ring has no certified path:
    /// `Ideal.groebner_basis_certified` raises ValueError on one.
    #[staticmethod]
    #[pyo3(text_signature = "(variables)")]
    fn rationals(py: Python<'_>, variables: Vec<String>) -> PyResult<Self> {
        let ring = PolynomialRing::rationals(variables).map_err(|e| of_ring(py, e))?;
        Ok(PyRing {
            inner: AnyRing::Rational(ring),
        })
    }

    /// The prime the coefficients live in, or None over the rational
    /// numbers.
    #[getter]
    fn modulus(&self) -> Option<u64> {
        match &self.inner {
            AnyRing::Prime(ring) => Some(ring.modulus()),
            AnyRing::Rational(_) => None,
        }
    }

    /// The variable names, largest variable first.
    #[getter]
    fn variables(&self) -> Vec<String> {
        match &self.inner {
            AnyRing::Prime(ring) => ring.variables().to_vec(),
            AnyRing::Rational(ring) => ring.variables().to_vec(),
        }
    }

    /// Read a polynomial from text.
    ///
    /// The syntax is the one `str(polynomial)` writes: terms separated by
    /// `+` and `-`, each a coefficient and a product of powers, as in
    /// `x^2*y - 3*z + 1`. Over `Q` a coefficient may name a fraction
    /// (`1/2*x`). Raises ParseError on text the ring cannot read, and
    /// RingError on a coefficient the domain has no value for.
    #[pyo3(text_signature = "($self, text)")]
    fn parse(&self, py: Python<'_>, text: &str) -> PyResult<PyPolynomial> {
        let inner = match &self.inner {
            AnyRing::Prime(ring) => {
                AnyPolynomial::Prime(ring.parse_polynomial(text).map_err(|e| of_ring(py, e))?)
            }
            AnyRing::Rational(ring) => {
                AnyPolynomial::Rational(ring.parse_polynomial(text).map_err(|e| of_ring(py, e))?)
            }
        };
        Ok(PyPolynomial { inner })
    }

    /// Build a polynomial from `(coefficient, exponents)` pairs.
    ///
    /// A coefficient is an int, a `fractions.Fraction`, or a
    /// `(numerator, denominator)` pair. Exponents hold one entry per
    /// variable, each between 0 and 65535. Repeated monomials add up and a
    /// term that reduces to zero drops out.
    #[pyo3(text_signature = "($self, terms)")]
    fn polynomial(
        &self,
        py: Python<'_>,
        terms: Vec<(Py<PyAny>, Vec<i64>)>,
    ) -> PyResult<PyPolynomial> {
        let mut read: Vec<(Coefficient, Vec<u16>)> = Vec::with_capacity(terms.len());
        for (coefficient, exponents) in &terms {
            read.push((
                coefficient_of(coefficient.bind(py))?,
                exponents_of(exponents)?,
            ));
        }
        let inner = match &self.inner {
            AnyRing::Prime(ring) => {
                AnyPolynomial::Prime(ring.polynomial(read).map_err(|e| of_ring(py, e))?)
            }
            AnyRing::Rational(ring) => {
                AnyPolynomial::Rational(ring.polynomial(read).map_err(|e| of_ring(py, e))?)
            }
        };
        Ok(PyPolynomial { inner })
    }

    /// The ideal these polynomials generate.
    ///
    /// Every polynomial must belong to this ring. Raises RingError
    /// otherwise.
    #[pyo3(text_signature = "($self, generators)")]
    fn ideal(&self, py: Python<'_>, generators: Vec<Py<PyPolynomial>>) -> PyResult<PyIdeal> {
        let inner = match &self.inner {
            AnyRing::Prime(ring) => {
                let mut collected = Vec::with_capacity(generators.len());
                for (index, generator) in generators.iter().enumerate() {
                    match &generator.get().inner {
                        AnyPolynomial::Prime(f) => collected.push(f.clone()),
                        AnyPolynomial::Rational(_) => return Err(foreign_polynomial(py, index)),
                    }
                }
                AnyIdeal::Prime(ring.ideal(collected).map_err(|e| of_ring(py, e))?)
            }
            AnyRing::Rational(ring) => {
                let mut collected = Vec::with_capacity(generators.len());
                for (index, generator) in generators.iter().enumerate() {
                    match &generator.get().inner {
                        AnyPolynomial::Rational(f) => collected.push(f.clone()),
                        AnyPolynomial::Prime(_) => return Err(foreign_polynomial(py, index)),
                    }
                }
                AnyIdeal::Rational(ring.ideal(collected).map_err(|e| of_ring(py, e))?)
            }
        };
        Ok(PyIdeal { inner })
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        match other.downcast::<PyRing>() {
            Ok(other) => self.inner == other.get().inner,
            Err(_) => false,
        }
    }

    fn __hash__(&self) -> isize {
        let mut hasher = DefaultHasher::new();
        match &self.inner {
            AnyRing::Prime(ring) => ring.hash(&mut hasher),
            AnyRing::Rational(ring) => ring.hash(&mut hasher),
        }
        hasher.finish() as isize
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            AnyRing::Prime(ring) => format!(
                "PolynomialRing.prime_field({}, {:?})",
                ring.modulus(),
                ring.variables()
            ),
            AnyRing::Rational(ring) => {
                format!("PolynomialRing.rationals({:?})", ring.variables())
            }
        }
    }
}

/// A polynomial of one ring.
///
/// Build one with `ring.parse(text)` or `ring.polynomial(terms)`. The value
/// is immutable, compares by content, and hashes by content.
#[pyclass(frozen, name = "Polynomial", module = "sylvester")]
struct PyPolynomial {
    inner: AnyPolynomial,
}

#[pymethods]
impl PyPolynomial {
    /// The terms, largest monomial first.
    ///
    /// Each item is `(coefficient, exponents)`. The coefficient is an int
    /// over `F_p`, in `[1, p)`, and a `fractions.Fraction` over `Q`. The
    /// exponents hold one entry per variable.
    #[pyo3(text_signature = "($self)")]
    fn terms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        match &self.inner {
            AnyPolynomial::Prime(f) => PyList::new(
                py,
                f.terms()
                    .map(|(coeff, exps)| (coeff.value(), exps.to_vec()))
                    .collect::<Vec<_>>(),
            ),
            AnyPolynomial::Rational(f) => PyList::new(
                py,
                f.terms()
                    .map(|(coeff, exps)| (coeff.clone(), exps.to_vec()))
                    .collect::<Vec<_>>(),
            ),
        }
    }

    /// The total degree, or None for the zero polynomial.
    #[pyo3(text_signature = "($self)")]
    fn degree(&self) -> Option<u32> {
        match &self.inner {
            AnyPolynomial::Prime(f) => f.degree(),
            AnyPolynomial::Rational(f) => f.degree(),
        }
    }

    /// Report whether the polynomial is zero.
    #[pyo3(text_signature = "($self)")]
    fn is_zero(&self) -> bool {
        match &self.inner {
            AnyPolynomial::Prime(f) => f.is_zero(),
            AnyPolynomial::Rational(f) => f.is_zero(),
        }
    }

    /// The ring the polynomial belongs to.
    #[getter]
    fn ring(&self) -> PyRing {
        PyRing {
            inner: match &self.inner {
                AnyPolynomial::Prime(f) => AnyRing::Prime(f.ring().clone()),
                AnyPolynomial::Rational(f) => AnyRing::Rational(f.ring().clone()),
            },
        }
    }

    fn __str__(&self) -> String {
        match &self.inner {
            AnyPolynomial::Prime(f) => f.to_string(),
            AnyPolynomial::Rational(f) => f.to_string(),
        }
    }

    fn __repr__(&self) -> String {
        format!("Polynomial({:?})", self.__str__())
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        match other.downcast::<PyPolynomial>() {
            Ok(other) => self.inner == other.get().inner,
            Err(_) => false,
        }
    }

    fn __hash__(&self) -> isize {
        hash_of(&self.inner)
    }
}

/// The ideal a list of polynomials generates.
///
/// Build one with `ring.ideal(generators)`.
#[pyclass(frozen, name = "Ideal", module = "sylvester")]
struct PyIdeal {
    inner: AnyIdeal,
}

#[pymethods]
impl PyIdeal {
    /// The generators, in the order the caller gave them.
    #[getter]
    fn generators(&self) -> Vec<PyPolynomial> {
        match &self.inner {
            AnyIdeal::Prime(ideal) => ideal
                .generators()
                .iter()
                .map(|f| PyPolynomial {
                    inner: AnyPolynomial::Prime(f.clone()),
                })
                .collect(),
            AnyIdeal::Rational(ideal) => ideal
                .generators()
                .iter()
                .map(|f| PyPolynomial {
                    inner: AnyPolynomial::Rational(f.clone()),
                })
                .collect(),
        }
    }

    /// The ring the ideal lives in.
    #[getter]
    fn ring(&self) -> PyRing {
        PyRing {
            inner: match &self.inner {
                AnyIdeal::Prime(ideal) => AnyRing::Prime(ideal.ring().clone()),
                AnyIdeal::Rational(ideal) => AnyRing::Rational(ideal.ring().clone()),
            },
        }
    }

    /// Report whether every generator is homogeneous.
    ///
    /// The test reads the generators, so it is sufficient and not
    /// necessary: a homogeneous ideal can be given by inhomogeneous
    /// generators. `GroebnerBasis.is_homogeneous` decides the ideal.
    #[pyo3(text_signature = "($self)")]
    fn has_homogeneous_generators(&self) -> bool {
        match &self.inner {
            AnyIdeal::Prime(ideal) => ideal.has_homogeneous_generators(),
            AnyIdeal::Rational(ideal) => ideal.has_homogeneous_generators(),
        }
    }

    /// The reduced Gröbner basis under grevlex.
    ///
    /// Every keyword defaults to None, which takes the library default.
    /// `backend` is `"f4"` or `"classic"`, `timeout` is seconds, and
    /// `memory_limit` is bytes. `threads` is the thread count of the run,
    /// and over `Q` it is the number of prime runs the driver starts at
    /// once. `stop` (`"unchanged"` or `"contains_input"`) and
    /// `extra_primes` describe the multimodular engine, so passing either
    /// on a prime-field ideal raises ValueError.
    ///
    /// The result carries no proof. Over `F_p` use
    /// `groebner_basis_certified` for a basis an independent verifier
    /// accepted. Over `Q` the value is the result of a heuristic modular
    /// computation, and `GroebnerBasis.lift` says what its stopping rule
    /// observed.
    #[pyo3(signature = (*, backend=None, timeout=None, memory_limit=None, threads=None, stop=None, extra_primes=None))]
    #[pyo3(
        text_signature = "($self, *, backend=None, timeout=None, memory_limit=None, threads=None, stop=None, extra_primes=None)"
    )]
    fn groebner_basis(
        slf: &Bound<'_, Self>,
        backend: Option<&str>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
        threads: Option<usize>,
        stop: Option<&str>,
        extra_primes: Option<usize>,
    ) -> PyResult<PyBasis> {
        let py = slf.py();
        let inner = match &slf.get().inner {
            AnyIdeal::Prime(ideal) => {
                let options =
                    compute_options(backend, timeout, memory_limit, threads, stop, extra_primes)?;
                let ideal = ideal.clone();
                AnyBasis::Prime(
                    py.allow_threads(|| ideal.groebner_basis(options))
                        .map_err(|e| of_compute(py, e))?,
                )
            }
            AnyIdeal::Rational(ideal) => {
                let options =
                    rational_options(backend, timeout, memory_limit, threads, stop, extra_primes)?;
                let ideal = ideal.clone();
                AnyBasis::Rational(
                    py.allow_threads(|| ideal.groebner_basis(options))
                        .map_err(|e| of_compute(py, e))?,
                )
            }
        };
        Ok(PyBasis { inner })
    }

    /// The reduced Gröbner basis and a report of the run.
    ///
    /// The basis is the one `groebner_basis` returns for the same
    /// keywords, and it carries no proof either.
    #[pyo3(signature = (*, backend=None, timeout=None, memory_limit=None, threads=None, stop=None, extra_primes=None))]
    #[pyo3(
        text_signature = "($self, *, backend=None, timeout=None, memory_limit=None, threads=None, stop=None, extra_primes=None)"
    )]
    fn groebner_basis_with_report(
        slf: &Bound<'_, Self>,
        backend: Option<&str>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
        threads: Option<usize>,
        stop: Option<&str>,
        extra_primes: Option<usize>,
    ) -> PyResult<(PyBasis, PyComputeReport)> {
        let py = slf.py();
        let (inner, report) = match &slf.get().inner {
            AnyIdeal::Prime(ideal) => {
                let options =
                    compute_options(backend, timeout, memory_limit, threads, stop, extra_primes)?;
                let ideal = ideal.clone();
                let (basis, report) = py
                    .allow_threads(|| ideal.groebner_basis_with_report(options))
                    .map_err(|e| of_compute(py, e))?;
                (AnyBasis::Prime(basis), report)
            }
            AnyIdeal::Rational(ideal) => {
                let options =
                    rational_options(backend, timeout, memory_limit, threads, stop, extra_primes)?;
                let ideal = ideal.clone();
                let (basis, report) = py
                    .allow_threads(|| ideal.groebner_basis_with_report(options))
                    .map_err(|e| of_compute(py, e))?;
                (AnyBasis::Rational(basis), report)
            }
        };
        Ok((PyBasis { inner }, PyComputeReport { inner: report }))
    }

    /// The reduced Gröbner basis and the certificate an independent
    /// verifier accepted.
    ///
    /// The value exists only after the verifier accepts the bytes, and the
    /// basis is decoded from them. The certificate format follows the
    /// backend: the classic backend writes `sylv-gb-cert-v1` and F4, the
    /// default, writes `sylv-gb-cert-v2`.
    ///
    /// There is no certified path over `Q`, so a rational ideal raises
    /// ValueError.
    #[pyo3(signature = (*, backend=None, timeout=None, memory_limit=None, threads=None))]
    #[pyo3(
        text_signature = "($self, *, backend=None, timeout=None, memory_limit=None, threads=None)"
    )]
    fn groebner_basis_certified(
        slf: &Bound<'_, Self>,
        backend: Option<&str>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
        threads: Option<usize>,
    ) -> PyResult<PyCertified> {
        let py = slf.py();
        match &slf.get().inner {
            AnyIdeal::Prime(ideal) => {
                let options = compute_options(backend, timeout, memory_limit, threads, None, None)?;
                let ideal = ideal.clone();
                let certified = py
                    .allow_threads(|| ideal.groebner_basis_certified(options))
                    .map_err(|e| of_certify(py, e))?;
                Ok(PyCertified { inner: certified })
            }
            AnyIdeal::Rational(_) => Err(PyValueError::new_err(
                "the rational numbers have no certified path: the lift is outside every certificate"
                    .to_string(),
            )),
        }
    }

    fn __repr__(&self) -> String {
        let count = match &self.inner {
            AnyIdeal::Prime(ideal) => ideal.generators().len(),
            AnyIdeal::Rational(ideal) => ideal.generators().len(),
        };
        format!("Ideal({count} generators)")
    }
}

/// A reduced Gröbner basis under grevlex.
///
/// The value is a sequence: `len`, indexing with negative indices and
/// slices, and iteration all work on it. There is no `__contains__`, so
/// `f in basis` iterates and compares: it asks whether `f` is one of the
/// listed polynomials. `basis.contains(f)` is the ideal membership test.
#[pyclass(frozen, name = "GroebnerBasis", module = "sylvester")]
struct PyBasis {
    inner: AnyBasis,
}

fn checked_prime_basis(
    py: Python<'_>,
    ring: &PolynomialRing<PrimeField>,
    polynomials: &[Py<PyPolynomial>],
    budget: Budget,
) -> PyResult<GroebnerBasis<PrimeField>> {
    let mut collected = Vec::with_capacity(polynomials.len());
    for (index, polynomial) in polynomials.iter().enumerate() {
        match &polynomial.get().inner {
            AnyPolynomial::Prime(value) => collected.push(value.clone()),
            AnyPolynomial::Rational(_) => return Err(foreign_polynomial(py, index)),
        }
    }
    let ring = ring.clone();
    py.allow_threads(|| GroebnerBasis::<PrimeField>::from_polynomials(&ring, collected, budget))
        .map_err(|error| of_basis(py, error))
}

fn checked_rational_basis(
    py: Python<'_>,
    ring: &PolynomialRing<Rationals>,
    polynomials: &[Py<PyPolynomial>],
    budget: Budget,
) -> PyResult<GroebnerBasis<Rationals>> {
    let mut collected = Vec::with_capacity(polynomials.len());
    for (index, polynomial) in polynomials.iter().enumerate() {
        match &polynomial.get().inner {
            AnyPolynomial::Rational(value) => collected.push(value.clone()),
            AnyPolynomial::Prime(_) => return Err(foreign_polynomial(py, index)),
        }
    }
    let ring = ring.clone();
    py.allow_threads(|| GroebnerBasis::<Rationals>::from_polynomials(&ring, collected, budget))
        .map_err(|error| of_basis(py, error))
}

#[pymethods]
impl PyBasis {
    /// Take a basis the caller supplies, after checking that it is one.
    ///
    /// The check runs under the budget: every polynomial belongs to the
    /// ring, none is zero, each is monic, the leading monomials run
    /// strictly descending, no monomial is divisible by the leading
    /// monomial of another element, and every S-polynomial reduces to
    /// zero. Raises BasisError with the first failure.
    ///
    /// A checked basis is a check and not a certificate.
    #[staticmethod]
    #[pyo3(signature = (ring, polynomials, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "(ring, polynomials, *, timeout=None, memory_limit=None)")]
    fn from_polynomials(
        py: Python<'_>,
        ring: &PyRing,
        polynomials: Vec<Py<PyPolynomial>>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyBasis> {
        let budget = budget_of(timeout, memory_limit)?;
        let inner = match &ring.inner {
            AnyRing::Prime(ring) => {
                AnyBasis::Prime(checked_prime_basis(py, ring, &polynomials, budget)?)
            }
            AnyRing::Rational(ring) => {
                AnyBasis::Rational(checked_rational_basis(py, ring, &polynomials, budget)?)
            }
        };
        Ok(PyBasis { inner })
    }

    /// The ring the basis lives in.
    #[getter]
    fn ring(&self) -> PyRing {
        PyRing {
            inner: match &self.inner {
                AnyBasis::Prime(basis) => AnyRing::Prime(basis.ring().clone()),
                AnyBasis::Rational(basis) => AnyRing::Rational(basis.ring().clone()),
            },
        }
    }

    /// Reduce `f` modulo the basis and return the remainder.
    ///
    /// The remainder holds no monomial divisible by the leading monomial of
    /// a basis element, so it is the unique normal form of `f` modulo the
    /// ideal. `timeout` is seconds and `memory_limit` is bytes. Raises
    /// RingError for a polynomial of another ring.
    #[pyo3(signature = (f, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, f, *, timeout=None, memory_limit=None)")]
    fn normal_form(
        &self,
        py: Python<'_>,
        f: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyPolynomial> {
        let budget = budget_of(timeout, memory_limit)?;
        let inner = match (&self.inner, &f.inner) {
            (AnyBasis::Prime(basis), AnyPolynomial::Prime(f)) => {
                let (basis, f) = (basis.clone(), f.clone());
                AnyPolynomial::Prime(
                    py.allow_threads(|| basis.normal_form(&f, budget))
                        .map_err(|e| of_normal_form(py, e))?,
                )
            }
            (AnyBasis::Rational(basis), AnyPolynomial::Rational(f)) => {
                let (basis, f) = (basis.clone(), f.clone());
                AnyPolynomial::Rational(
                    py.allow_threads(|| basis.normal_form(&f, budget))
                        .map_err(|e| of_normal_form(py, e))?,
                )
            }
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyPolynomial { inner })
    }

    /// Report whether `f` belongs to the ideal this basis generates.
    ///
    /// A polynomial belongs to that ideal exactly when its normal form is
    /// zero, so this carries the same budget and the same errors as
    /// `normal_form`.
    ///
    /// The answer is about the ideal this basis generates, which is not
    /// always the ideal the caller started from. `from_polynomials` takes
    /// the basis whole, so there the two ideals are the same. Over Q the
    /// driver is a heuristic and `lift()["established"]` says what holds.
    /// Under "unchanged" nothing relates the two ideals, so neither
    /// answer is established for the input ideal. Under "contains_input"
    /// the ideal of this basis contains the input ideal, so False holds
    /// for the input ideal as well, and True does not.
    #[pyo3(signature = (f, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, f, *, timeout=None, memory_limit=None)")]
    fn contains(
        &self,
        py: Python<'_>,
        f: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<bool> {
        let budget = budget_of(timeout, memory_limit)?;
        match (&self.inner, &f.inner) {
            (AnyBasis::Prime(basis), AnyPolynomial::Prime(f)) => {
                let (basis, f) = (basis.clone(), f.clone());
                py.allow_threads(|| basis.contains(&f, budget))
                    .map_err(|e| of_normal_form(py, e))
            }
            (AnyBasis::Rational(basis), AnyPolynomial::Rational(f)) => {
                let (basis, f) = (basis.clone(), f.clone());
                py.allow_threads(|| basis.contains(&f, budget))
                    .map_err(|e| of_normal_form(py, e))
            }
            _ => Err(foreign_argument(py)),
        }
    }

    /// The Hilbert series of the quotient by the leading monomial ideal
    /// of this basis.
    ///
    /// `HilbertSeries` states what the value says about the ideal this
    /// basis generates. That ideal is not always the ideal the caller
    /// started from: over Q, `lift()["established"]` says what relates
    /// the two. Under "unchanged" nothing does, and under
    /// "contains_input" the series is that of an ideal containing the
    /// input ideal.
    ///
    /// The recursion behind the numerator branches and this release
    /// offers no bound on its work, so the budget stops it at a recursion
    /// node.
    #[pyo3(signature = (*, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, *, timeout=None, memory_limit=None)")]
    fn hilbert_series(
        &self,
        py: Python<'_>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyHilbertSeries> {
        let budget = budget_of(timeout, memory_limit)?;
        let series = match &self.inner {
            AnyBasis::Prime(basis) => {
                let basis = basis.clone();
                py.allow_threads(|| basis.hilbert_series(budget))
            }
            AnyBasis::Rational(basis) => {
                let basis = basis.clone();
                py.allow_threads(|| basis.hilbert_series(budget))
            }
        };
        Ok(PyHilbertSeries {
            inner: series.map_err(|e| of_hilbert(py, e))?,
        })
    }

    /// The Krull dimension of the quotient by the ideal this basis
    /// generates, or None for the unit ideal.
    ///
    /// The value is the one `hilbert_series().dimension()` returns, with
    /// the same budget and the same errors. Over Q, under
    /// "contains_input" the ideal of this basis contains the input ideal,
    /// so the value is at most the Krull dimension of the quotient by the
    /// input ideal. Under "unchanged" nothing follows. None is the unit
    /// ideal, whose quotient ring is zero and whose dimension is -1 by
    /// the usual convention.
    #[pyo3(signature = (*, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, *, timeout=None, memory_limit=None)")]
    fn krull_dimension(
        &self,
        py: Python<'_>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Option<usize>> {
        Ok(self
            .hilbert_series(py, timeout, memory_limit)?
            .inner
            .dimension())
    }

    /// Report whether the ideal the basis generates is homogeneous.
    ///
    /// The test reads every element, and it decides the ideal: under a
    /// graded order an ideal is homogeneous exactly when its reduced
    /// Gröbner basis is.
    #[pyo3(text_signature = "($self)")]
    fn is_homogeneous(&self) -> bool {
        match &self.inner {
            AnyBasis::Prime(basis) => basis.is_homogeneous(),
            AnyBasis::Rational(basis) => basis.is_homogeneous(),
        }
    }

    /// The record of the multimodular run that produced the basis.
    ///
    /// The keys are `primes_consumed`, `primes_skipped`, `primes_folded`,
    /// `primes_discarded`, `confirming_primes`, `modulus_bits`, and
    /// `established` (`"unchanged"` or `"contains_input"`). The counters
    /// describe the run and establish nothing on their own.
    ///
    /// Returns None for a rational basis the checked constructor accepted,
    /// because no modular run produced one. Raises ValueError on a
    /// prime-field basis, where no such record exists at all.
    #[pyo3(text_signature = "($self)")]
    fn lift<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        match &self.inner {
            AnyBasis::Prime(_) => Err(PyValueError::new_err(
                "a prime-field basis has no modular lift: the engine computed it over F_p"
                    .to_string(),
            )),
            AnyBasis::Rational(basis) => basis.lift().map(|lift| lift_dict(py, lift)).transpose(),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__(slf: &Bound<'_, Self>, index: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        let basis = &slf.get().inner;
        let length = basis.len();
        if let Ok(slice) = index.downcast::<PySlice>() {
            let indices = slice.indices(length as isize)?;
            let mut taken = Vec::new();
            let mut position = indices.start;
            for _ in 0..indices.slicelength {
                taken.push(PyPolynomial {
                    inner: basis.get(position as usize),
                });
                position += indices.step;
            }
            return Ok(PyList::new(py, taken)?.into_any().unbind());
        }
        let mut position = index.extract::<isize>()?;
        if position < 0 {
            position += length as isize;
        }
        if position < 0 || position as usize >= length {
            return Err(PyIndexError::new_err("basis index out of range"));
        }
        Ok(Py::new(
            py,
            PyPolynomial {
                inner: basis.get(position as usize),
            },
        )?
        .into_any())
    }

    fn __iter__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyIterator>> {
        let py = slf.py();
        let basis = &slf.get().inner;
        let items: Vec<PyPolynomial> = (0..basis.len())
            .map(|index| PyPolynomial {
                inner: basis.get(index),
            })
            .collect();
        PyList::new(py, items)?.try_iter()
    }

    fn __repr__(&self) -> String {
        format!("GroebnerBasis({} polynomials)", self.inner.len())
    }
}

/// A basis and the certificate bytes an independent verifier accepted.
#[pyclass(frozen, name = "CertifiedGroebnerBasis", module = "sylvester")]
struct PyCertified {
    inner: CertifiedGroebnerBasis,
}

#[pymethods]
impl PyCertified {
    /// The basis the verifier accepted, decoded from the certificate.
    #[getter]
    fn basis(&self) -> PyBasis {
        PyBasis {
            inner: AnyBasis::Prime(self.inner.basis().clone()),
        }
    }

    /// The certificate bytes, which `sylvester.verify` reads back.
    #[getter]
    fn certificate<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, self.inner.certificate())
    }

    fn __repr__(&self) -> String {
        format!(
            "CertifiedGroebnerBasis({} polynomials, {} certificate bytes)",
            self.inner.basis().len(),
            self.inner.certificate().len()
        )
    }
}

/// The Hilbert series of the quotient by the leading monomial ideal.
///
/// The value is `N(t) / (1 - t)^n`, with `n` the number of variables. For a
/// homogeneous ideal it is the Hilbert series of the quotient by the ideal
/// itself. For any ideal the coefficient of `t^d` is the first difference
/// of the affine Hilbert function, and the dimension read off the series is
/// the Krull dimension of the quotient.
#[pyclass(frozen, name = "HilbertSeries", module = "sylvester")]
struct PyHilbertSeries {
    inner: HilbertSeries,
}

#[pymethods]
impl PyHilbertSeries {
    /// The numerator, dense and lowest degree first.
    ///
    /// The last coefficient is nonzero, and the empty list is the zero
    /// numerator, which is the unit ideal.
    #[getter]
    fn numerator(&self) -> Vec<BigInt> {
        self.inner.numerator().to_vec()
    }

    /// The power of `1 - t` in the denominator, which is the number of
    /// variables, before any cancellation.
    #[getter]
    fn denominator_power(&self) -> usize {
        self.inner.denominator_power()
    }

    /// The Krull dimension of the quotient ring, or None for the zero ring.
    #[pyo3(text_signature = "($self)")]
    fn dimension(&self) -> Option<usize> {
        self.inner.dimension()
    }

    /// The multiplicity of the quotient ring, or None for the zero ring.
    ///
    /// It is the value at `t = 1` of the numerator after every factor of
    /// `1 - t` is cancelled.
    #[pyo3(text_signature = "($self)")]
    fn multiplicity(&self) -> Option<BigInt> {
        self.inner.multiplicity()
    }

    /// The dimension over the field of the piece of degree `degree`.
    #[pyo3(text_signature = "($self, degree)")]
    fn coefficient(&self, degree: u32) -> BigInt {
        self.inner.coefficient(degree)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        format!("HilbertSeries({:?})", self.inner.to_string())
    }
}

fn lift_dict<'py>(py: Python<'py>, lift: &ModularLift) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("primes_consumed", lift.primes_consumed)?;
    dict.set_item("primes_skipped", lift.primes_skipped)?;
    dict.set_item("primes_folded", lift.primes_folded)?;
    dict.set_item("primes_discarded", lift.primes_discarded)?;
    dict.set_item("confirming_primes", lift.confirming_primes)?;
    dict.set_item("modulus_bits", lift.modulus_bits)?;
    dict.set_item(
        "established",
        match lift.established {
            Established::Unchanged => "unchanged",
            Established::ContainsInput => "contains_input",
        },
    )?;
    Ok(dict)
}

/// What one computation did, next to the basis it produced.
///
/// The report is for a benchmark harness or a test. Nothing in the package
/// reads it back.
#[pyclass(frozen, name = "ComputeReport", module = "sylvester")]
struct PyComputeReport {
    inner: ComputeReport,
}

#[pymethods]
impl PyComputeReport {
    /// The backend that ran: `"f4"` or `"classic"`.
    #[getter]
    fn backend(&self) -> &'static str {
        match self.inner.backend {
            Backend::F4 => "f4",
            Backend::Classic => "classic",
        }
    }

    /// What the F4 run counted, or None under the classic backend and over
    /// `Q`.
    ///
    /// A rational run has many prime runs, and summing their counters would
    /// name a run that never happened.
    #[getter]
    fn counters<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let Some(counters) = self.inner.counters else {
            return Ok(None);
        };
        let dict = PyDict::new(py);
        set_matrix_counters(&dict, &counters)?;
        set_pair_counters(&dict, &counters)?;
        Ok(Some(dict))
    }

    /// The wall time of the engine call, in seconds.
    #[getter]
    fn elapsed(&self) -> f64 {
        self.inner.elapsed.as_secs_f64()
    }

    /// The number of threads the run had.
    #[getter]
    fn threads_used(&self) -> usize {
        self.inner.threads_used
    }

    /// What the multimodular run did, or None over a prime field.
    ///
    /// The keys are the ones `GroebnerBasis.lift` returns.
    #[getter]
    fn modular<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        self.inner
            .modular
            .as_ref()
            .map(|lift| lift_dict(py, lift))
            .transpose()
    }

    /// The number of prime runs the driver ran at once, or None over a
    /// prime field.
    #[getter]
    fn modular_concurrency(&self) -> Option<usize> {
        self.inner.modular_concurrency
    }

    fn __repr__(&self) -> String {
        format!(
            "ComputeReport(backend={:?}, elapsed={:.6}, threads_used={})",
            self.backend(),
            self.elapsed(),
            self.threads_used()
        )
    }
}

fn set_matrix_counters(dict: &Bound<'_, PyDict>, counters: &F4Counters) -> PyResult<()> {
    dict.set_item("batches", counters.batches)?;
    dict.set_item("matrix_rows", counters.matrix_rows)?;
    dict.set_item("matrix_columns", counters.matrix_columns)?;
    dict.set_item("matrix_nonzeros", counters.matrix_nonzeros)?;
    dict.set_item("zero_rows", counters.zero_rows)?;
    dict.set_item("new_pivots", counters.new_pivots)?;
    dict.set_item("lane_restarts", counters.lane_restarts)?;
    dict.set_item("batch_retries", counters.batch_retries)?;
    dict.set_item("basis_monomials", counters.basis_monomials)?;
    Ok(())
}

fn set_pair_counters(dict: &Bound<'_, PyDict>, counters: &F4Counters) -> PyResult<()> {
    dict.set_item("pairs_generated", counters.pairs_generated)?;
    dict.set_item("pairs_discarded_product", counters.pairs_discarded_product)?;
    dict.set_item("pairs_discarded_b", counters.pairs_discarded_b)?;
    dict.set_item("pairs_discarded_m", counters.pairs_discarded_m)?;
    dict.set_item("pairs_discarded_f", counters.pairs_discarded_f)?;
    Ok(())
}

/// A basis an independent verifier accepted, with the input it was checked
/// against.
///
/// Neither certificate schema carries variable names, so the polynomials
/// print under the synthetic names `x1 .. xn`. They are positions, not the
/// names the caller used.
#[pyclass(frozen, name = "VerifiedGroebnerBasis", module = "sylvester")]
struct PyVerified {
    modulus: u64,
    nvars: usize,
    input: Vec<Polynomial<PrimeField>>,
    basis: Vec<Polynomial<PrimeField>>,
}

#[pymethods]
impl PyVerified {
    /// The prime the certificate names. The verifier proved it prime.
    #[getter]
    fn modulus(&self) -> u64 {
        self.modulus
    }

    /// The number of variables.
    #[getter]
    fn nvars(&self) -> usize {
        self.nvars
    }

    /// The input polynomials, in certificate order.
    #[getter]
    fn input(&self) -> Vec<PyPolynomial> {
        self.input
            .iter()
            .map(|f| PyPolynomial {
                inner: AnyPolynomial::Prime(f.clone()),
            })
            .collect()
    }

    /// The reduced basis, sorted strictly descending by leading monomial.
    #[getter]
    fn basis(&self) -> Vec<PyPolynomial> {
        self.basis
            .iter()
            .map(|f| PyPolynomial {
                inner: AnyPolynomial::Prime(f.clone()),
            })
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "VerifiedGroebnerBasis(modulus={}, nvars={}, {} polynomials)",
            self.modulus,
            self.nvars,
            self.basis.len()
        )
    }
}

/// Check certificate bytes with the independent verifier.
///
/// Acceptance proves that the input and the basis inside the certificate
/// generate the same ideal, and that the basis is the reduced Gröbner basis
/// of that ideal under grevlex. It says nothing about who wrote the bytes.
/// The verifier shares no code with the engines.
///
/// Raises CertificateInvalid when the bytes do not hold, and Timeout or
/// LimitExceeded when the verifier ran out of its own budget, which says
/// nothing about the certificate.
#[pyfunction(name = "verify")]
#[pyo3(text_signature = "(data)")]
fn verify_certificate(py: Python<'_>, data: Vec<u8>) -> PyResult<PyVerified> {
    let verified = py
        .allow_threads(|| verify::verify(&data))
        .map_err(|e| of_verify(py, e))?;
    let names: Vec<String> = (1..=verified.nvars())
        .map(|index| format!("x{index}"))
        .collect();
    let ring =
        PolynomialRing::prime_field(verified.modulus(), names).map_err(|e| of_ring(py, e))?;
    let read = |polynomials: &[verify::Poly]| -> PyResult<Vec<Polynomial<PrimeField>>> {
        polynomials
            .iter()
            .map(|f| {
                let terms: Vec<(u64, Vec<u16>)> = f
                    .terms()
                    .iter()
                    .map(|term| {
                        let exps = term
                            .mono()
                            .exps()
                            .iter()
                            .map(|&exp| u16::try_from(exp))
                            .collect::<Result<Vec<u16>, _>>()
                            .map_err(|_| {
                                raise(
                                    py,
                                    &exceptions(py).limit_exceeded,
                                    format!(
                                        "the certificate holds an exponent above {MAX_EXPONENT}"
                                    ),
                                )
                            })?;
                        Ok((term.coeff(), exps))
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                ring.polynomial(terms).map_err(|e| of_ring(py, e))
            })
            .collect()
    };
    Ok(PyVerified {
        modulus: verified.modulus(),
        nvars: verified.nvars(),
        input: read(verified.input())?,
        basis: read(verified.basis())?,
    })
}

/// Gröbner bases over prime fields and over the rational numbers, under the
/// grevlex order, with independent certificate verifiers.
#[pymodule(name = "_sylvester")]
fn sylvester_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let classes = build_exceptions(py)?;
    add_exception_classes(m, &classes)?;
    EXCEPTIONS
        .set(py, classes)
        .map_err(|_| PyRuntimeError::new_err("the exception classes are already built"))?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    add_value_classes(m)?;
    m.add_function(wrap_pyfunction!(verify_certificate, m)?)?;
    Ok(())
}

fn add_exception_classes(m: &Bound<'_, PyModule>, classes: &Exceptions) -> PyResult<()> {
    let py = m.py();
    for (name, class) in [
        ("SylvesterError", &classes.base),
        ("RingError", &classes.ring),
        ("ParseError", &classes.parse),
        ("BasisError", &classes.basis),
        ("CertificateInvalid", &classes.certificate_invalid),
        ("BudgetExhausted", &classes.budget_exhausted),
        ("Timeout", &classes.timeout),
        ("MemoryLimitExceeded", &classes.memory_limit),
        ("LimitExceeded", &classes.limit_exceeded),
        ("InternalDefect", &classes.internal_defect),
    ] {
        m.add(name, class.bind(py))?;
    }
    Ok(())
}

fn add_value_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRing>()?;
    m.add_class::<PyPolynomial>()?;
    m.add_class::<PyIdeal>()?;
    m.add_class::<PyBasis>()?;
    m.add_class::<PyCertified>()?;
    m.add_class::<PyHilbertSeries>()?;
    m.add_class::<PyComputeReport>()?;
    m.add_class::<PyVerified>()?;
    Ok(())
}
