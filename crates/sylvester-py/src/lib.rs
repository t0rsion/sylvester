//! Python bindings for the `sylvester` crate: reduced Gröbner bases under
//! grevlex, over a prime field or over the rational numbers.
//!
//! The Rust library carries the coefficient domain as a type parameter.
//! Python learns the domain at run time, so each pair of Rust
//! instantiations collapses to one Python class holding an enum over the
//! two, and the ring constructor picks the arm. `docs/rational-design.md`
//! fixes the base surface and the rule around the GIL.
//!
//! Every computation releases the GIL. Python values are converted while the
//! GIL is held, then immutable Rust values are borrowed or moved into a
//! joined worker. No borrowed Python value crosses into the released region.
//!
//! The calling Python main thread polls pending signals while the worker
//! runs. Python delivers signals to that thread; calls from other Python
//! threads still release the GIL but do not receive signal exceptions there.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::num::NonZeroUsize;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::mpsc::{RecvTimeoutError, sync_channel};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::{ToPrimitive, Zero};
use pyo3::exceptions::{PyException, PyIndexError, PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::sync::GILOnceCell;
use pyo3::types::{PyBytes, PyDict, PyIterator, PyList, PySlice, PyTuple, PyType};

use sylvester::verify::{self, VerifyError};
use sylvester::{
    ArithmeticError, Backend, BasisError, Budget, CancellationToken, CertifiedGroebnerBasis,
    CertifyError, ClaimedProvenance, Coefficient, ComputeError, ComputeOptions, ComputeReport,
    EnvelopeDomain, EnvelopeError, EqualityCheckError, Established, ExpressionError,
    FiniteQuotient, GroebnerBasis, HilbertError, HilbertSeries, Ideal, ModularLift,
    MultiplicationMatrix, NormalFormError, ParseError, Polynomial, PolynomialRing, PrimeField,
    QuotientError, RationalEqualityCheck, RationalOptions, RationalStop, Rationals, ResultEnvelope,
    RingError, UnivariatePolynomial,
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
    let leaves = build_leaf_exceptions(py, &base_type, &budget_type, &runtime_error)?;
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
    runtime_error: &Bound<'_, PyType>,
) -> PyResult<LeafExceptions> {
    let value_error = py.get_type::<PyValueError>();
    Ok(LeafExceptions {
        ring: new_exception(
            py,
            "RingError",
            &[base, &value_error],
            "A construction through the ring stopped: a modulus that is not \
             prime, a variable name the parser cannot read, an exponent \
             count that is not the variable count, or a polynomial of \
             another ring.",
        )?,
        parse: new_exception(
            py,
            "ParseError",
            &[base, &value_error],
            "The text is not a polynomial of the ring. The message names \
             the byte position and what the parser found there.",
        )?,
        basis: new_exception(
            py,
            "BasisError",
            &[base, &value_error],
            "A supplied list is not a reduced Gröbner basis. `index` names \
             the element that failed a shape check, and `left` and `right` \
             the first pair whose S-polynomial has a nonzero remainder.",
        )?,
        certificate_invalid: new_exception(
            py,
            "CertificateInvalid",
            &[base, &value_error],
            "The verifier rejected the certificate bytes. The bytes are \
             untrusted input, so this is a ValueError and says nothing \
             about the engines.",
        )?,
        timeout: new_exception(
            py,
            "Timeout",
            &[budget],
            "The deadline passed before the computation finished.",
        )?,
        memory_limit: new_exception(
            py,
            "MemoryLimitExceeded",
            &[budget],
            "The live data of the computation passed the memory limit.",
        )?,
        limit_exceeded: new_exception(
            py,
            "LimitExceeded",
            &[base, runtime_error],
            "The computation reached a structural limit of this release: a \
             degree or an exponent past the width of one exponent, a full \
             monomial table, an exhausted prime sequence, or a certificate \
             cap. `limit` carries the bound where the limit names one.",
        )?,
        internal_defect: new_exception(
            py,
            "InternalDefect",
            &[base, runtime_error],
            "The crate contradicted itself: it wrote a certificate its own \
             verifier rejected, or the emitter reported a fault. This is a \
             defect in sylvester, never bad input.",
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

fn of_expression(py: Python<'_>, error: ExpressionError) -> PyErr {
    let classes = exceptions(py);
    let message = error.to_string();
    match error {
        ExpressionError::Parse(inner) => of_parse(py, inner),
        ExpressionError::Ring(inner) => of_ring(py, inner),
        ExpressionError::Timeout => raise(py, &classes.timeout, message),
        ExpressionError::MemoryLimitExceeded => raise(py, &classes.memory_limit, message),
        ExpressionError::ExponentLimit { limit } => {
            raise_with(py, &classes.limit_exceeded, message, |value| {
                value.setattr("limit", limit)
            })
        }
        ExpressionError::NestingLimit { limit, .. } => {
            raise_with(py, &classes.limit_exceeded, message, |value| {
                value.setattr("limit", limit)
            })
        }
        ExpressionError::Compute(inner) => of_compute(py, inner),
        ExpressionError::Arithmetic(inner) => operation_error(py, inner),
    }
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

fn of_quotient(py: Python<'_>, error: QuotientError) -> PyErr {
    let message = error.to_string();
    match error {
        QuotientError::RingMismatch => raise(py, &exceptions(py).ring, message),
        QuotientError::ExponentLimit { limit } => {
            raise_with(py, &exceptions(py).limit_exceeded, message, |value| {
                value.setattr("limit", limit)
            })
        }
        QuotientError::Timeout => raise(py, &exceptions(py).timeout, message),
        QuotientError::MemoryLimitExceeded => raise(py, &exceptions(py).memory_limit, message),
        QuotientError::NotFinite => PyValueError::new_err(message),
    }
}

fn of_envelope(py: Python<'_>, error: EnvelopeError) -> PyErr {
    match error {
        EnvelopeError::RingMismatch | EnvelopeError::Ring(_) => {
            raise(py, &exceptions(py).ring, error.to_string())
        }
        EnvelopeError::Expression(inner) => of_expression(py, inner),
        EnvelopeError::Verification(inner) => of_verify(py, inner),
        EnvelopeError::EqualityCheck(inner) => of_equality(py, inner),
        EnvelopeError::Timeout => raise(py, &exceptions(py).timeout, error.to_string()),
        EnvelopeError::MemoryLimitExceeded => {
            raise(py, &exceptions(py).memory_limit, error.to_string())
        }
        EnvelopeError::CertificateMismatch | EnvelopeError::Format(_) => {
            raise(py, &exceptions(py).certificate_invalid, error.to_string())
        }
    }
}

fn of_equality(py: Python<'_>, error: EqualityCheckError) -> PyErr {
    let message = error.to_string();
    match error {
        EqualityCheckError::RingMismatch => raise(py, &exceptions(py).ring, message),
        EqualityCheckError::Candidate(inner) => of_basis(py, inner),
        EqualityCheckError::ReverseMembership { generator } => {
            raise_basis(py, message, Some(generator), None)
        }
        EqualityCheckError::ForwardMembership { basis } => {
            raise_basis(py, message, Some(basis), None)
        }
        EqualityCheckError::ExponentLimit { limit } => {
            raise_with(py, &exceptions(py).limit_exceeded, message, |value| {
                value.setattr("limit", limit)
            })
        }
        EqualityCheckError::Timeout => raise(py, &exceptions(py).timeout, message),
        EqualityCheckError::MemoryLimitExceeded => raise(py, &exceptions(py).memory_limit, message),
    }
}

fn of_certify(py: Python<'_>, error: CertifyError) -> PyErr {
    let classes = exceptions(py);
    let message = error.to_string();
    match error {
        CertifyError::Engine(inner) | CertifyError::WriterExhausted(inner) => of_compute(py, inner),
        CertifyError::VerifierExhausted(inner) => of_exhausted_verifier(py, inner),
        CertifyError::CapExceeded { limit, .. } => {
            raise_with(py, &classes.limit_exceeded, message, |value| {
                value.setattr("limit", limit)
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

fn cancellable_budget(
    timeout: Option<f64>,
    memory_limit: Option<usize>,
    cancellation: &CancellationToken,
) -> PyResult<Budget> {
    Ok(budget_of(timeout, memory_limit)?.cancellation(cancellation.clone()))
}

fn run_with_signals<T, E, F, C>(py: Python<'_>, cancel: C, work: F) -> Result<Result<T, E>, PyErr>
where
    T: Send,
    E: Send,
    F: FnOnce() -> Result<T, E> + Send,
    C: Fn(),
{
    let (sender, receiver) = sync_channel::<Result<Result<T, E>, ()>>(1);
    let receiver = Arc::new(Mutex::new(receiver));
    thread::scope(|scope| {
        let spawned = thread::Builder::new().spawn_scoped(scope, move || {
            let result = catch_unwind(AssertUnwindSafe(work));
            let message = match result {
                Ok(result) => Ok(result),
                Err(_) => Err(()),
            };
            let _ = sender.send(message);
        });
        if let Err(error) = spawned {
            return Err(PyRuntimeError::new_err(format!(
                "the worker could not start: {error}"
            )));
        }
        loop {
            let receiver_for_wait = Arc::clone(&receiver);
            let received = py.allow_threads(move || match receiver_for_wait.lock() {
                Ok(receiver) => receiver.recv_timeout(Duration::from_millis(20)),
                Err(_) => Err(RecvTimeoutError::Disconnected),
            });
            match received {
                Ok(Ok(result)) => return Ok(result),
                Ok(Err(())) => {
                    return Err(PyRuntimeError::new_err("the worker panicked"));
                }
                Err(RecvTimeoutError::Timeout) => {
                    if let Err(signal) = py.check_signals() {
                        cancel();
                        let receiver_for_cleanup = Arc::clone(&receiver);
                        let _ = py.allow_threads(move || {
                            receiver_for_cleanup
                                .lock()
                                .ok()
                                .and_then(|receiver| receiver.recv().ok())
                        });
                        return Err(signal);
                    }
                }
                Err(RecvTimeoutError::Disconnected) => {
                    return Err(PyRuntimeError::new_err("the worker ended without a result"));
                }
            }
        }
    })
}

/// Run one budgeted operation with a signal-aware worker.
fn run_budgeted<T, E, F>(
    py: Python<'_>,
    timeout: Option<f64>,
    memory_limit: Option<usize>,
    work: F,
) -> PyResult<Result<T, E>>
where
    T: Send,
    E: Send,
    F: FnOnce(Budget) -> Result<T, E> + Send,
{
    let cancellation = CancellationToken::new();
    let budget = cancellable_budget(timeout, memory_limit, &cancellation)?;
    let cancel = cancellation.clone();
    run_with_signals(py, move || cancel.cancel(), move || work(budget))
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
        None | Some("contains_input") => Ok(RationalStop::ContainsInput { extra }),
        Some("unchanged") => Ok(RationalStop::Unchanged { extra }),
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
    cancellation: Option<&CancellationToken>,
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
    let mut budget = budget_of(timeout, memory_limit)?;
    if let Some(token) = cancellation {
        budget = budget.cancellation((*token).clone());
    }
    let mut options = ComputeOptions::new().budget(budget);
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
    cancellation: Option<&CancellationToken>,
) -> PyResult<RationalOptions> {
    let mut budget = budget_of(timeout, memory_limit)?;
    if let Some(token) = cancellation {
        budget = budget.cancellation((*token).clone());
    }
    let mut options = ComputeOptions::new().budget(budget);
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

fn exponents_of(exponents: &[BigInt]) -> PyResult<Vec<u16>> {
    exponents
        .iter()
        .map(|exponent| {
            let Some(exponent) = exponent.to_u64() else {
                return Err(PyValueError::new_err(format!(
                    "an exponent is between 0 and {MAX_EXPONENT}"
                )));
            };
            u16::try_from(exponent).map_err(|_| {
                PyValueError::new_err(format!("an exponent is between 0 and {MAX_EXPONENT}"))
            })
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

fn ring_of(polynomial: &AnyPolynomial) -> AnyRing {
    match polynomial {
        AnyPolynomial::Prime(value) => AnyRing::Prime(value.ring().clone()),
        AnyPolynomial::Rational(value) => AnyRing::Rational(value.ring().clone()),
    }
}

fn polynomial_from_terms(
    template: &AnyPolynomial,
    terms: Vec<(Coefficient, Vec<u16>)>,
) -> Result<AnyPolynomial, RingError> {
    match template {
        AnyPolynomial::Prime(value) => Ok(AnyPolynomial::Prime(value.ring().polynomial(terms)?)),
        AnyPolynomial::Rational(value) => {
            Ok(AnyPolynomial::Rational(value.ring().polynomial(terms)?))
        }
    }
}

type OperationBudget = Budget;

fn operation_error(py: Python<'_>, error: ArithmeticError) -> PyErr {
    match error {
        ArithmeticError::RingMismatch => foreign_argument(py),
        ArithmeticError::CoefficientConversion(error) => of_ring(py, error),
        ArithmeticError::ExponentLimit { limit } => raise_with(
            py,
            &exceptions(py).limit_exceeded,
            error.to_string(),
            |value| value.setattr("limit", limit),
        ),
        ArithmeticError::Timeout => raise(py, &exceptions(py).timeout, error.to_string()),
        ArithmeticError::MemoryLimitExceeded => {
            raise(py, &exceptions(py).memory_limit, error.to_string())
        }
    }
}

#[derive(Clone, Copy)]
enum BinaryOperation {
    Add,
    Subtract,
    Multiply,
    ReverseSubtract,
}

#[derive(Clone, Copy)]
enum UnsupportedOperand {
    ReturnNotImplemented,
    RaiseTypeError,
}

fn apply_binary(
    left: &AnyPolynomial,
    right: &AnyPolynomial,
    operation: BinaryOperation,
    budget: OperationBudget,
) -> Result<AnyPolynomial, ArithmeticError> {
    match operation {
        BinaryOperation::Add => add_polynomials(left, right, budget),
        BinaryOperation::Subtract => subtract_polynomials(left, right, budget),
        BinaryOperation::Multiply => multiply_polynomials(left, right, budget),
        BinaryOperation::ReverseSubtract => subtract_polynomials(right, left, budget),
    }
}

fn binary_operation(
    py: Python<'_>,
    left: &AnyPolynomial,
    other: &Bound<'_, PyAny>,
    operation: BinaryOperation,
    unsupported: UnsupportedOperand,
    timeout: Option<f64>,
    memory_limit: Option<usize>,
) -> PyResult<Py<PyAny>> {
    let right = match other.downcast::<PyPolynomial>() {
        Ok(value) => value.get().inner.clone(),
        Err(_) => {
            let Some(scalar) = scalar_from(other)? else {
                return match unsupported {
                    UnsupportedOperand::ReturnNotImplemented => Ok(py.NotImplemented()),
                    UnsupportedOperand::RaiseTypeError => Err(PyTypeError::new_err(
                        "the operand is not a polynomial or an integer, Fraction, or (numerator, denominator) scalar",
                    )),
                };
            };
            scalar_polynomial(left, &scalar).map_err(|error| of_ring(py, error))?
        }
    };
    if ring_of(left) != ring_of(&right) {
        return Err(foreign_argument(py));
    }
    let result = run_budgeted(py, timeout, memory_limit, |budget| {
        apply_binary(left, &right, operation, budget)
    })?
    .map_err(|error| operation_error(py, error))?;
    polynomial_result(py, result)
}

fn power_result(
    py: Python<'_>,
    value: &AnyPolynomial,
    exponent: u64,
    timeout: Option<f64>,
    memory_limit: Option<usize>,
) -> PyResult<Py<PyAny>> {
    let result = run_budgeted(py, timeout, memory_limit, |budget| {
        power_polynomial(value, exponent, budget)
    })?
    .map_err(|error| operation_error(py, error))?;
    polynomial_result(py, result)
}

fn negated_polynomial(
    value: &AnyPolynomial,
    budget: OperationBudget,
) -> Result<AnyPolynomial, ArithmeticError> {
    match value {
        AnyPolynomial::Prime(value) => Ok(AnyPolynomial::Prime(value.try_neg(budget)?)),
        AnyPolynomial::Rational(value) => Ok(AnyPolynomial::Rational(value.try_neg(budget)?)),
    }
}

fn negation_result(
    py: Python<'_>,
    value: &AnyPolynomial,
    timeout: Option<f64>,
    memory_limit: Option<usize>,
) -> PyResult<Py<PyAny>> {
    let result = run_budgeted(py, timeout, memory_limit, |budget| {
        negated_polynomial(value, budget)
    })?
    .map_err(|error| operation_error(py, error))?;
    polynomial_result(py, result)
}

fn power_polynomial(
    value: &AnyPolynomial,
    exponent: u64,
    budget: OperationBudget,
) -> Result<AnyPolynomial, ArithmeticError> {
    let exponent = u32::try_from(exponent).map_err(|_| ArithmeticError::ExponentLimit {
        limit: u16::MAX as u32,
    })?;
    match value {
        AnyPolynomial::Prime(value) => Ok(AnyPolynomial::Prime(value.try_pow(exponent, budget)?)),
        AnyPolynomial::Rational(value) => {
            Ok(AnyPolynomial::Rational(value.try_pow(exponent, budget)?))
        }
    }
}

fn nonnegative_exponent_from(value: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    let exponent = value
        .extract::<BigInt>()
        .map_err(|_| PyValueError::new_err("the polynomial exponent is a nonnegative integer"))?;
    if exponent.sign() == Sign::Minus {
        return Err(PyValueError::new_err(
            "the polynomial exponent is a nonnegative integer",
        ));
    }
    Ok(exponent)
}

fn polynomial_exponent_from(py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<u64> {
    let exponent = nonnegative_exponent_from(value)?;
    exponent.try_into().map_err(|_| {
        raise_with(
            py,
            &exceptions(py).limit_exceeded,
            "the polynomial exponent is above 65535".to_string(),
            |error| error.setattr("limit", MAX_EXPONENT),
        )
    })
}

fn quotient_exponent_from(value: &Bound<'_, PyAny>) -> PyResult<usize> {
    let exponent = nonnegative_exponent_from(value)?;
    usize::try_from(exponent).map_err(|_| {
        PyValueError::new_err("the quotient exponent is a nonnegative integer that fits usize")
    })
}

fn index_from(value: &Bound<'_, PyAny>) -> PyResult<Option<usize>> {
    let index = value
        .extract::<BigInt>()
        .map_err(|_| PyTypeError::new_err("an index is an integer"))?;
    if index.sign() == Sign::Minus {
        return Ok(None);
    }
    Ok(usize::try_from(index).ok())
}

fn add_polynomials(
    left: &AnyPolynomial,
    right: &AnyPolynomial,
    budget: OperationBudget,
) -> Result<AnyPolynomial, ArithmeticError> {
    match (left, right) {
        (AnyPolynomial::Prime(left), AnyPolynomial::Prime(right)) => {
            Ok(AnyPolynomial::Prime(left.try_add(right, budget)?))
        }
        (AnyPolynomial::Rational(left), AnyPolynomial::Rational(right)) => {
            Ok(AnyPolynomial::Rational(left.try_add(right, budget)?))
        }
        _ => Err(ArithmeticError::RingMismatch),
    }
}

fn subtract_polynomials(
    left: &AnyPolynomial,
    right: &AnyPolynomial,
    budget: OperationBudget,
) -> Result<AnyPolynomial, ArithmeticError> {
    match (left, right) {
        (AnyPolynomial::Prime(left), AnyPolynomial::Prime(right)) => {
            Ok(AnyPolynomial::Prime(left.try_sub(right, budget)?))
        }
        (AnyPolynomial::Rational(left), AnyPolynomial::Rational(right)) => {
            Ok(AnyPolynomial::Rational(left.try_sub(right, budget)?))
        }
        _ => Err(ArithmeticError::RingMismatch),
    }
}

fn multiply_polynomials(
    left: &AnyPolynomial,
    right: &AnyPolynomial,
    budget: OperationBudget,
) -> Result<AnyPolynomial, ArithmeticError> {
    match (left, right) {
        (AnyPolynomial::Prime(left), AnyPolynomial::Prime(right)) => {
            Ok(AnyPolynomial::Prime(left.try_mul(right, budget)?))
        }
        (AnyPolynomial::Rational(left), AnyPolynomial::Rational(right)) => {
            Ok(AnyPolynomial::Rational(left.try_mul(right, budget)?))
        }
        _ => Err(ArithmeticError::RingMismatch),
    }
}

fn scalar_polynomial(
    template: &AnyPolynomial,
    scalar: &Coefficient,
) -> Result<AnyPolynomial, RingError> {
    let nvars = match template {
        AnyPolynomial::Prime(value) => value.ring().nvars(),
        AnyPolynomial::Rational(value) => value.ring().nvars(),
    };
    polynomial_from_terms(template, vec![(scalar.clone(), vec![0; nvars])])
}

fn scalar_from(value: &Bound<'_, PyAny>) -> PyResult<Option<Coefficient>> {
    if value.extract::<BigInt>().is_ok()
        || value.extract::<BigRational>().is_ok()
        || value.extract::<(BigInt, BigInt)>().is_ok()
    {
        return Ok(Some(coefficient_of(value)?));
    }
    Ok(None)
}

fn polynomial_result(py: Python<'_>, inner: AnyPolynomial) -> PyResult<Py<PyAny>> {
    Ok(Py::new(py, PyPolynomial { inner })?.into_any())
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

#[derive(Clone, Debug)]
enum AnyQuotient {
    Prime(FiniteQuotient<PrimeField>),
    Rational(FiniteQuotient<Rationals>),
}

#[derive(Clone, Debug)]
enum AnyMatrix {
    Prime(MultiplicationMatrix<PrimeField>),
    Rational(MultiplicationMatrix<Rationals>),
}

#[derive(Clone, Debug)]
enum AnyUnivariate {
    Prime(UnivariatePolynomial<PrimeField>),
    Rational(UnivariatePolynomial<Rationals>),
}

fn ring_of_quotient(quotient: &AnyQuotient) -> AnyRing {
    match quotient {
        AnyQuotient::Prime(value) => AnyRing::Prime(value.ring().clone()),
        AnyQuotient::Rational(value) => AnyRing::Rational(value.ring().clone()),
    }
}

fn basis_of_quotient(quotient: &AnyQuotient) -> AnyBasis {
    match quotient {
        AnyQuotient::Prime(value) => AnyBasis::Prime(value.source_basis().clone()),
        AnyQuotient::Rational(value) => AnyBasis::Rational(value.source_basis().clone()),
    }
}

fn ring_of_matrix(matrix: &AnyMatrix) -> AnyRing {
    match matrix {
        AnyMatrix::Prime(value) => AnyRing::Prime(value.ring().clone()),
        AnyMatrix::Rational(value) => AnyRing::Rational(value.ring().clone()),
    }
}

fn ring_of_univariate(polynomial: &AnyUnivariate) -> AnyRing {
    match polynomial {
        AnyUnivariate::Prime(value) => AnyRing::Prime(value.ring().clone()),
        AnyUnivariate::Rational(value) => AnyRing::Rational(value.ring().clone()),
    }
}

fn prime_coefficient_object<'py>(
    py: Python<'py>,
    coefficient: &sylvester::Felt,
) -> PyResult<Py<PyAny>> {
    let list = PyList::new(py, [coefficient.value()])?;
    Ok(list.get_item(0)?.unbind())
}

fn rational_coefficient_object<'py>(
    py: Python<'py>,
    coefficient: &BigRational,
) -> PyResult<Py<PyAny>> {
    let list = PyList::new(py, [coefficient])?;
    Ok(list.get_item(0)?.unbind())
}

fn prime_zero() -> &'static sylvester::Felt {
    static ZERO: std::sync::OnceLock<sylvester::Felt> = std::sync::OnceLock::new();
    ZERO.get_or_init(sylvester::Felt::default)
}

fn rational_zero() -> &'static BigRational {
    static ZERO: std::sync::OnceLock<BigRational> = std::sync::OnceLock::new();
    ZERO.get_or_init(BigRational::zero)
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

    /// The degree-one variable polynomials, in ring order.
    #[getter]
    fn gens(&self) -> Vec<PyPolynomial> {
        match &self.inner {
            AnyRing::Prime(ring) => ring
                .generators()
                .into_iter()
                .map(|value| PyPolynomial {
                    inner: AnyPolynomial::Prime(value),
                })
                .collect(),
            AnyRing::Rational(ring) => ring
                .generators()
                .into_iter()
                .map(|value| PyPolynomial {
                    inner: AnyPolynomial::Rational(value),
                })
                .collect(),
        }
    }

    /// The zero polynomial of this ring.
    #[getter]
    fn zero(&self) -> PyPolynomial {
        match &self.inner {
            AnyRing::Prime(ring) => PyPolynomial {
                inner: AnyPolynomial::Prime(ring.zero()),
            },
            AnyRing::Rational(ring) => PyPolynomial {
                inner: AnyPolynomial::Rational(ring.zero()),
            },
        }
    }

    /// The constant one polynomial of this ring.
    #[getter]
    fn one(&self) -> PyPolynomial {
        match &self.inner {
            AnyRing::Prime(ring) => PyPolynomial {
                inner: AnyPolynomial::Prime(ring.one()),
            },
            AnyRing::Rational(ring) => PyPolynomial {
                inner: AnyPolynomial::Rational(ring.one()),
            },
        }
    }

    /// Read a polynomial from text.
    ///
    /// The syntax is the one `str(polynomial)` writes: terms separated by
    /// `+` and `-`, each a coefficient and a product of powers, as in
    /// `x^2*y - 3*z + 1`. Over `Q` a coefficient may name a fraction
    /// (`1/2*x`). Raises ParseError on text the ring cannot read, and
    /// RingError on a coefficient the domain has no value for.
    #[pyo3(signature = (text, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, text, *, timeout=None, memory_limit=None)")]
    fn parse(
        &self,
        py: Python<'_>,
        text: &str,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyPolynomial> {
        let text = text.to_owned();
        let inner = match &self.inner {
            AnyRing::Prime(ring) => AnyPolynomial::Prime(
                run_budgeted(py, timeout, memory_limit, |budget| {
                    ring.parse_polynomial_with_budget(&text, budget)
                })?
                .map_err(|e| of_expression(py, e))?,
            ),
            AnyRing::Rational(ring) => AnyPolynomial::Rational(
                run_budgeted(py, timeout, memory_limit, |budget| {
                    ring.parse_polynomial_with_budget(&text, budget)
                })?
                .map_err(|e| of_expression(py, e))?,
            ),
        };
        Ok(PyPolynomial { inner })
    }

    /// Build a polynomial from `(coefficient, exponents)` pairs.
    ///
    /// A coefficient is an int, a `fractions.Fraction`, or a
    /// `(numerator, denominator)` pair. Exponents hold one entry per
    /// variable, each between 0 and 65535. Repeated monomials add up and a
    /// term that reduces to zero drops out. A negative or larger Python
    /// integer raises ValueError, including one too large for a machine
    /// integer.
    #[pyo3(text_signature = "($self, terms)")]
    fn polynomial(
        &self,
        py: Python<'_>,
        terms: Vec<(Py<PyAny>, Vec<BigInt>)>,
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

    /// Add another polynomial or a scalar.
    fn __add__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        binary_operation(
            py,
            &self.inner,
            other,
            BinaryOperation::Add,
            UnsupportedOperand::ReturnNotImplemented,
            None,
            None,
        )
    }

    /// Add another polynomial or a scalar from the right.
    fn __radd__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        binary_operation(
            py,
            &self.inner,
            other,
            BinaryOperation::Add,
            UnsupportedOperand::ReturnNotImplemented,
            None,
            None,
        )
    }

    /// Subtract another polynomial or a scalar.
    fn __sub__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        binary_operation(
            py,
            &self.inner,
            other,
            BinaryOperation::Subtract,
            UnsupportedOperand::ReturnNotImplemented,
            None,
            None,
        )
    }

    /// Subtract a polynomial from a scalar or another polynomial.
    fn __rsub__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        binary_operation(
            py,
            &self.inner,
            other,
            BinaryOperation::ReverseSubtract,
            UnsupportedOperand::ReturnNotImplemented,
            None,
            None,
        )
    }

    /// Multiply by another polynomial or a scalar.
    fn __mul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        binary_operation(
            py,
            &self.inner,
            other,
            BinaryOperation::Multiply,
            UnsupportedOperand::ReturnNotImplemented,
            None,
            None,
        )
    }

    /// Multiply a polynomial by this value from the right.
    fn __rmul__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        binary_operation(
            py,
            &self.inner,
            other,
            BinaryOperation::Multiply,
            UnsupportedOperand::ReturnNotImplemented,
            None,
            None,
        )
    }

    /// Negate the polynomial.
    fn __neg__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        negation_result(py, &self.inner, None, None)
    }

    /// Raise the polynomial to a nonnegative integer power.
    fn __pow__(
        &self,
        py: Python<'_>,
        exponent: &Bound<'_, PyAny>,
        modulo: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        if modulo.is_some() {
            return Ok(py.NotImplemented());
        }
        let exponent = polynomial_exponent_from(py, exponent)?;
        power_result(py, &self.inner, exponent, None, None)
    }

    /// Negate the polynomial under an optional budget.
    #[pyo3(signature = (*, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, *, timeout=None, memory_limit=None)")]
    fn neg(
        &self,
        py: Python<'_>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Py<PyAny>> {
        negation_result(py, &self.inner, timeout, memory_limit)
    }

    /// Add another polynomial or scalar under an optional budget.
    #[pyo3(signature = (other, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, other, *, timeout=None, memory_limit=None)")]
    fn add(
        &self,
        py: Python<'_>,
        other: &Bound<'_, PyAny>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Py<PyAny>> {
        binary_operation(
            py,
            &self.inner,
            other,
            BinaryOperation::Add,
            UnsupportedOperand::RaiseTypeError,
            timeout,
            memory_limit,
        )
    }

    /// Subtract another polynomial or scalar under an optional budget.
    #[pyo3(signature = (other, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, other, *, timeout=None, memory_limit=None)")]
    fn sub(
        &self,
        py: Python<'_>,
        other: &Bound<'_, PyAny>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Py<PyAny>> {
        binary_operation(
            py,
            &self.inner,
            other,
            BinaryOperation::Subtract,
            UnsupportedOperand::RaiseTypeError,
            timeout,
            memory_limit,
        )
    }

    /// Multiply by another polynomial or scalar under an optional budget.
    #[pyo3(signature = (other, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, other, *, timeout=None, memory_limit=None)")]
    fn mul(
        &self,
        py: Python<'_>,
        other: &Bound<'_, PyAny>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Py<PyAny>> {
        binary_operation(
            py,
            &self.inner,
            other,
            BinaryOperation::Multiply,
            UnsupportedOperand::RaiseTypeError,
            timeout,
            memory_limit,
        )
    }

    /// Raise the polynomial to a nonnegative integer power under a budget.
    #[pyo3(signature = (exponent, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, exponent, *, timeout=None, memory_limit=None)")]
    fn pow(
        &self,
        py: Python<'_>,
        exponent: &Bound<'_, PyAny>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Py<PyAny>> {
        let exponent = polynomial_exponent_from(py, exponent)?;
        power_result(py, &self.inner, exponent, timeout, memory_limit)
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
    /// `stop` defaults to `"contains_input"` over `Q`. `extra_primes`
    /// describes the multimodular engine, so passing either rational-only
    /// keyword on a prime-field ideal raises ValueError.
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
    // Each argument is a public Python keyword. A Rust options value would
    // replace that Python call shape rather than shorten it.
    #[allow(clippy::too_many_arguments)]
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
                let cancellation = CancellationToken::new();
                let options = compute_options(
                    backend,
                    timeout,
                    memory_limit,
                    threads,
                    stop,
                    extra_primes,
                    Some(&cancellation),
                )?;
                compute_prime_basis(py, ideal, options, cancellation)?
            }
            AnyIdeal::Rational(ideal) => {
                let cancellation = CancellationToken::new();
                let options = rational_options(
                    backend,
                    timeout,
                    memory_limit,
                    threads,
                    stop,
                    extra_primes,
                    Some(&cancellation),
                )?;
                compute_rational_basis(py, ideal, options, cancellation)?
            }
        };
        Ok(PyBasis { inner })
    }

    /// The reduced Gröbner basis and the report of the engine run.
    #[pyo3(signature = (*, backend=None, timeout=None, memory_limit=None, threads=None, stop=None, extra_primes=None))]
    #[pyo3(
        text_signature = "($self, *, backend=None, timeout=None, memory_limit=None, threads=None, stop=None, extra_primes=None)"
    )]
    #[allow(clippy::too_many_arguments)]
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
                let cancellation = CancellationToken::new();
                let options = compute_options(
                    backend,
                    timeout,
                    memory_limit,
                    threads,
                    stop,
                    extra_primes,
                    Some(&cancellation),
                )?;
                compute_prime_basis_with_report(py, ideal, options, cancellation)?
            }
            AnyIdeal::Rational(ideal) => {
                let cancellation = CancellationToken::new();
                let options = rational_options(
                    backend,
                    timeout,
                    memory_limit,
                    threads,
                    stop,
                    extra_primes,
                    Some(&cancellation),
                )?;
                compute_rational_basis_with_report(py, ideal, options, cancellation)?
            }
        };
        Ok((PyBasis { inner }, report))
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
                let cancellation = CancellationToken::new();
                let options = compute_options(
                    backend,
                    timeout,
                    memory_limit,
                    threads,
                    None,
                    None,
                    Some(&cancellation),
                )?;
                let certified = run_with_signals(
                    py,
                    || cancellation.cancel(),
                    move || ideal.groebner_basis_certified(options),
                )?
                .map_err(|e| of_certify(py, e))?;
                Ok(PyCertified { inner: certified })
            }
            AnyIdeal::Rational(_) => Err(PyValueError::new_err(
                "the rational numbers have no certified path: the lift is outside every certificate"
                    .to_string(),
            )),
        }
    }

    /// Check exact equality between a rational input ideal and a candidate.
    ///
    /// The check validates the candidate, proves both ideal inclusions with
    /// tracked rational cofactors, and retains the source data. It is a
    /// library check, not an independent certificate.
    #[pyo3(signature = (basis, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, basis, *, timeout=None, memory_limit=None)")]
    fn check_basis_equality(
        &self,
        py: Python<'_>,
        basis: &PyBasis,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyRationalEqualityCheck> {
        let (ideal, basis) = match (&self.inner, &basis.inner) {
            (AnyIdeal::Rational(ideal), AnyBasis::Rational(basis)) => (ideal, basis),
            (AnyIdeal::Rational(_), _) => return Err(foreign_argument(py)),
            (AnyIdeal::Prime(_), _) => {
                return Err(PyValueError::new_err(
                    "exact ideal equality checks are available over the rational numbers only",
                ));
            }
        };
        let checked = run_budgeted(py, timeout, memory_limit, move |budget| {
            ideal.check_basis_equality(basis, budget)
        })?
        .map_err(|error| of_equality(py, error))?;
        Ok(PyRationalEqualityCheck { inner: checked })
    }

    fn __repr__(&self) -> String {
        let count = match &self.inner {
            AnyIdeal::Prime(ideal) => ideal.generators().len(),
            AnyIdeal::Rational(ideal) => ideal.generators().len(),
        };
        format!("Ideal({count} generators)")
    }
}

fn compute_prime_basis(
    py: Python<'_>,
    ideal: &Ideal<PrimeField>,
    options: ComputeOptions,
    cancellation: CancellationToken,
) -> PyResult<AnyBasis> {
    let basis = run_with_signals(
        py,
        || cancellation.cancel(),
        move || ideal.groebner_basis(options),
    )?
    .map_err(|error| of_compute(py, error))?;
    Ok(AnyBasis::Prime(basis))
}

fn compute_rational_basis(
    py: Python<'_>,
    ideal: &Ideal<Rationals>,
    options: RationalOptions,
    cancellation: CancellationToken,
) -> PyResult<AnyBasis> {
    let basis = run_with_signals(
        py,
        || cancellation.cancel(),
        move || ideal.groebner_basis(options),
    )?
    .map_err(|error| of_compute(py, error))?;
    Ok(AnyBasis::Rational(basis))
}

fn compute_prime_basis_with_report(
    py: Python<'_>,
    ideal: &Ideal<PrimeField>,
    options: ComputeOptions,
    cancellation: CancellationToken,
) -> PyResult<(AnyBasis, PyComputeReport)> {
    let (basis, report) = run_with_signals(
        py,
        || cancellation.cancel(),
        move || ideal.groebner_basis_with_report(options),
    )?
    .map_err(|error| of_compute(py, error))?;
    Ok((AnyBasis::Prime(basis), PyComputeReport { inner: report }))
}

fn compute_rational_basis_with_report(
    py: Python<'_>,
    ideal: &Ideal<Rationals>,
    options: RationalOptions,
    cancellation: CancellationToken,
) -> PyResult<(AnyBasis, PyComputeReport)> {
    let (basis, report) = run_with_signals(
        py,
        || cancellation.cancel(),
        move || ideal.groebner_basis_with_report(options),
    )?
    .map_err(|error| of_compute(py, error))?;
    Ok((AnyBasis::Rational(basis), PyComputeReport { inner: report }))
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
        let inner = match &ring.inner {
            AnyRing::Prime(ring) => {
                let mut collected = Vec::with_capacity(polynomials.len());
                for (index, polynomial) in polynomials.iter().enumerate() {
                    match &polynomial.get().inner {
                        AnyPolynomial::Prime(f) => collected.push(f.clone()),
                        AnyPolynomial::Rational(_) => return Err(foreign_polynomial(py, index)),
                    }
                }
                AnyBasis::Prime(
                    run_budgeted(py, timeout, memory_limit, move |budget| {
                        GroebnerBasis::<PrimeField>::from_polynomials(ring, collected, budget)
                    })?
                    .map_err(|e| of_basis(py, e))?,
                )
            }
            AnyRing::Rational(ring) => {
                let mut collected = Vec::with_capacity(polynomials.len());
                for (index, polynomial) in polynomials.iter().enumerate() {
                    match &polynomial.get().inner {
                        AnyPolynomial::Rational(f) => collected.push(f.clone()),
                        AnyPolynomial::Prime(_) => return Err(foreign_polynomial(py, index)),
                    }
                }
                AnyBasis::Rational(
                    run_budgeted(py, timeout, memory_limit, move |budget| {
                        GroebnerBasis::<Rationals>::from_polynomials(ring, collected, budget)
                    })?
                    .map_err(|e| of_basis(py, e))?,
                )
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
        let inner = match (&self.inner, &f.inner) {
            (AnyBasis::Prime(basis), AnyPolynomial::Prime(f)) => AnyPolynomial::Prime(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    basis.normal_form(f, budget)
                })?
                .map_err(|e| of_normal_form(py, e))?,
            ),
            (AnyBasis::Rational(basis), AnyPolynomial::Rational(f)) => AnyPolynomial::Rational(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    basis.normal_form(f, budget)
                })?
                .map_err(|e| of_normal_form(py, e))?,
            ),
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyPolynomial { inner })
    }

    /// Report whether this basis defines a finite quotient algebra.
    fn is_zero_dimensional(&self) -> bool {
        match &self.inner {
            AnyBasis::Prime(basis) => basis.is_zero_dimensional(),
            AnyBasis::Rational(basis) => basis.is_zero_dimensional(),
        }
    }

    /// Build the finite quotient algebra defined by this basis.
    #[pyo3(signature = (*, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, *, timeout=None, memory_limit=None)")]
    fn finite_quotient(
        &self,
        py: Python<'_>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyFiniteQuotient> {
        let inner = match &self.inner {
            AnyBasis::Prime(basis) => AnyQuotient::Prime(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    basis.finite_quotient(budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            AnyBasis::Rational(basis) => AnyQuotient::Rational(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    basis.finite_quotient(budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
        };
        Ok(PyFiniteQuotient { inner })
    }

    /// Divide `f` by every basis element and return the quotients and remainder.
    #[pyo3(signature = (f, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, f, *, timeout=None, memory_limit=None)")]
    fn divide(
        &self,
        py: Python<'_>,
        f: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<(Vec<PyPolynomial>, PyPolynomial)> {
        match (&self.inner, &f.inner) {
            (AnyBasis::Prime(basis), AnyPolynomial::Prime(value)) => {
                let result = run_budgeted(py, timeout, memory_limit, move |budget| {
                    basis.divide(value, budget)
                })?
                .map_err(|error| of_normal_form(py, error))?;
                let quotients = result
                    .quotients
                    .into_iter()
                    .map(|quotient| PyPolynomial {
                        inner: AnyPolynomial::Prime(quotient),
                    })
                    .collect();
                Ok((
                    quotients,
                    PyPolynomial {
                        inner: AnyPolynomial::Prime(result.remainder),
                    },
                ))
            }
            (AnyBasis::Rational(basis), AnyPolynomial::Rational(value)) => {
                let result = run_budgeted(py, timeout, memory_limit, move |budget| {
                    basis.divide(value, budget)
                })?
                .map_err(|error| of_normal_form(py, error))?;
                let quotients = result
                    .quotients
                    .into_iter()
                    .map(|quotient| PyPolynomial {
                        inner: AnyPolynomial::Rational(quotient),
                    })
                    .collect();
                Ok((
                    quotients,
                    PyPolynomial {
                        inner: AnyPolynomial::Rational(result.remainder),
                    },
                ))
            }
            _ => Err(foreign_argument(py)),
        }
    }

    /// Divide `f` by this basis and return `(quotients, remainder)`.
    #[pyo3(signature = (f, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, f, *, timeout=None, memory_limit=None)")]
    fn divmod(
        &self,
        py: Python<'_>,
        f: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<(Vec<PyPolynomial>, PyPolynomial)> {
        self.divide(py, f, timeout, memory_limit)
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
        match (&self.inner, &f.inner) {
            (AnyBasis::Prime(basis), AnyPolynomial::Prime(f)) => {
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    basis.contains(f, budget)
                })?
                .map_err(|e| of_normal_form(py, e))
            }
            (AnyBasis::Rational(basis), AnyPolynomial::Rational(f)) => {
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    basis.contains(f, budget)
                })?
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
        let series = match &self.inner {
            AnyBasis::Prime(basis) => run_budgeted(py, timeout, memory_limit, move |budget| {
                basis.hilbert_series(budget)
            })?,
            AnyBasis::Rational(basis) => run_budgeted(py, timeout, memory_limit, move |budget| {
                basis.hilbert_series(budget)
            })?,
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

/// A successful exact rational ideal equality check with origin data.
#[pyclass(frozen, name = "RationalEqualityCheck", module = "sylvester")]
struct PyRationalEqualityCheck {
    inner: RationalEqualityCheck,
}

#[pymethods]
impl PyRationalEqualityCheck {
    /// The input ideal retained by the check.
    #[getter]
    fn input(&self) -> PyIdeal {
        PyIdeal {
            inner: AnyIdeal::Rational(self.inner.input().clone()),
        }
    }

    /// The candidate basis retained by the check.
    #[getter]
    fn basis(&self) -> PyBasis {
        PyBasis {
            inner: AnyBasis::Rational(self.inner.basis().clone()),
        }
    }

    /// The cofactor rows, in basis and input generator order.
    #[getter]
    fn origins(&self) -> Vec<Vec<PyPolynomial>> {
        self.inner
            .origins()
            .iter()
            .map(|row| {
                row.iter()
                    .map(|polynomial| PyPolynomial {
                        inner: AnyPolynomial::Rational(polynomial.clone()),
                    })
                    .collect()
            })
            .collect()
    }

    /// Counters from the exact check.
    #[getter]
    fn metrics<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let metrics = self.inner.metrics();
        let result = PyDict::new(py);
        result.set_item("raw_basis_elements", metrics.raw_basis_elements)?;
        result.set_item("raw_basis_terms", metrics.raw_basis_terms)?;
        result.set_item("final_working_bytes", metrics.final_working_bytes)?;
        Ok(result)
    }

    fn __repr__(&self) -> String {
        format!(
            "RationalEqualityCheck({} basis elements, {} origin rows)",
            self.inner.basis().len(),
            self.inner.origins().len()
        )
    }
}

fn prime_optional_coefficients<'py, I>(
    py: Python<'py>,
    coefficients: I,
) -> PyResult<Bound<'py, PyList>>
where
    I: IntoIterator<Item = Option<u64>>,
    I::IntoIter: ExactSizeIterator,
{
    PyList::new(py, coefficients.into_iter().map(|value| value.unwrap_or(0)))
}

fn rational_optional_coefficients<'py, I>(
    py: Python<'py>,
    coefficients: I,
) -> PyResult<Bound<'py, PyList>>
where
    I: IntoIterator<Item = Option<BigRational>>,
    I::IntoIter: ExactSizeIterator,
{
    PyList::new(
        py,
        coefficients
            .into_iter()
            .map(|value| value.unwrap_or_else(BigRational::zero)),
    )
}

/// A finite-dimensional quotient algebra defined by a Gröbner basis.
#[pyclass(frozen, name = "FiniteQuotient", module = "sylvester")]
struct PyFiniteQuotient {
    inner: AnyQuotient,
}

#[pymethods]
impl PyFiniteQuotient {
    /// The polynomial ring of the quotient.
    #[getter]
    fn ring(&self) -> PyRing {
        PyRing {
            inner: ring_of_quotient(&self.inner),
        }
    }

    /// The basis that defines this quotient, including rational lift data.
    #[getter]
    fn source_basis(&self) -> PyBasis {
        PyBasis {
            inner: basis_of_quotient(&self.inner),
        }
    }

    /// The standard monomials as exponent vectors in ring order.
    #[getter]
    fn basis(&self) -> Vec<Vec<u16>> {
        match &self.inner {
            AnyQuotient::Prime(value) => value.standard_monomials().to_vec(),
            AnyQuotient::Rational(value) => value.standard_monomials().to_vec(),
        }
    }

    /// The standard monomials as exponent vectors in ring order.
    #[getter]
    fn standard_monomials(&self) -> Vec<Vec<u16>> {
        self.basis()
    }

    /// The vector-space dimension of the quotient.
    #[getter]
    fn dimension(&self) -> usize {
        match &self.inner {
            AnyQuotient::Prime(value) => value.vector_space_dimension(),
            AnyQuotient::Rational(value) => value.vector_space_dimension(),
        }
    }

    /// The vector-space dimension of the quotient.
    fn vector_space_dimension(&self) -> usize {
        self.dimension()
    }

    /// Reduce a polynomial to its residue class.
    #[pyo3(signature = (polynomial, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, polynomial, *, timeout=None, memory_limit=None)")]
    fn reduce(
        &self,
        py: Python<'_>,
        polynomial: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyPolynomial> {
        let inner = match (&self.inner, &polynomial.inner) {
            (AnyQuotient::Prime(quotient), AnyPolynomial::Prime(value)) => AnyPolynomial::Prime(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.reduce(value, budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            (AnyQuotient::Rational(quotient), AnyPolynomial::Rational(value)) => {
                AnyPolynomial::Rational(
                    run_budgeted(py, timeout, memory_limit, move |budget| {
                        quotient.reduce(value, budget)
                    })?
                    .map_err(|error| of_quotient(py, error))?,
                )
            }
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyPolynomial { inner })
    }

    /// Return coordinates of a residue in the standard monomial basis.
    #[pyo3(signature = (polynomial, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, polynomial, *, timeout=None, memory_limit=None)")]
    fn coordinates<'py>(
        &self,
        py: Python<'py>,
        polynomial: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Bound<'py, PyList>> {
        match (&self.inner, &polynomial.inner) {
            (AnyQuotient::Prime(quotient), AnyPolynomial::Prime(value)) => {
                let coordinates = run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.coordinates(value, budget)
                })?
                .map_err(|error| of_quotient(py, error))?;
                prime_optional_coefficients(
                    py,
                    coordinates
                        .into_iter()
                        .map(|coefficient| coefficient.map(|value| value.value())),
                )
            }
            (AnyQuotient::Rational(quotient), AnyPolynomial::Rational(value)) => {
                let coordinates = run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.coordinates(value, budget)
                })?
                .map_err(|error| of_quotient(py, error))?;
                rational_optional_coefficients(py, coordinates)
            }
            _ => Err(foreign_argument(py)),
        }
    }

    /// Add two residue classes and reduce the result.
    #[pyo3(signature = (left, right, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, left, right, *, timeout=None, memory_limit=None)")]
    fn add(
        &self,
        py: Python<'_>,
        left: &PyPolynomial,
        right: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyPolynomial> {
        let inner = match (&self.inner, &left.inner, &right.inner) {
            (
                AnyQuotient::Prime(quotient),
                AnyPolynomial::Prime(left),
                AnyPolynomial::Prime(right),
            ) => AnyPolynomial::Prime(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.add(left, right, budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            (
                AnyQuotient::Rational(quotient),
                AnyPolynomial::Rational(left),
                AnyPolynomial::Rational(right),
            ) => AnyPolynomial::Rational(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.add(left, right, budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyPolynomial { inner })
    }

    /// Multiply two residue classes and reduce the result.
    #[pyo3(signature = (left, right, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, left, right, *, timeout=None, memory_limit=None)")]
    fn multiply(
        &self,
        py: Python<'_>,
        left: &PyPolynomial,
        right: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyPolynomial> {
        let inner = match (&self.inner, &left.inner, &right.inner) {
            (
                AnyQuotient::Prime(quotient),
                AnyPolynomial::Prime(left),
                AnyPolynomial::Prime(right),
            ) => AnyPolynomial::Prime(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.multiply(left, right, budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            (
                AnyQuotient::Rational(quotient),
                AnyPolynomial::Rational(left),
                AnyPolynomial::Rational(right),
            ) => AnyPolynomial::Rational(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.multiply(left, right, budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyPolynomial { inner })
    }

    /// Raise a residue class to a nonnegative integer power.
    #[pyo3(signature = (polynomial, exponent, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, polynomial, exponent, *, timeout=None, memory_limit=None)")]
    fn pow(
        &self,
        py: Python<'_>,
        polynomial: &PyPolynomial,
        exponent: &Bound<'_, PyAny>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyPolynomial> {
        let exponent = quotient_exponent_from(exponent)?;
        let inner = match (&self.inner, &polynomial.inner) {
            (AnyQuotient::Prime(quotient), AnyPolynomial::Prime(value)) => AnyPolynomial::Prime(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.pow(value, exponent, budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            (AnyQuotient::Rational(quotient), AnyPolynomial::Rational(value)) => {
                AnyPolynomial::Rational(
                    run_budgeted(py, timeout, memory_limit, move |budget| {
                        quotient.pow(value, exponent, budget)
                    })?
                    .map_err(|error| of_quotient(py, error))?,
                )
            }
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyPolynomial { inner })
    }

    /// Build the matrix of multiplication by a residue class.
    #[pyo3(signature = (polynomial, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, polynomial, *, timeout=None, memory_limit=None)")]
    fn multiplication_matrix(
        &self,
        py: Python<'_>,
        polynomial: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyMultiplicationMatrix> {
        let inner = match (&self.inner, &polynomial.inner) {
            (AnyQuotient::Prime(quotient), AnyPolynomial::Prime(value)) => AnyMatrix::Prime(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.multiplication_matrix(value, budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            (AnyQuotient::Rational(quotient), AnyPolynomial::Rational(value)) => {
                AnyMatrix::Rational(
                    run_budgeted(py, timeout, memory_limit, move |budget| {
                        quotient.multiplication_matrix(value, budget)
                    })?
                    .map_err(|error| of_quotient(py, error))?,
                )
            }
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyMultiplicationMatrix { inner })
    }

    /// Return the characteristic polynomial of multiplication by a residue class.
    #[pyo3(signature = (polynomial, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, polynomial, *, timeout=None, memory_limit=None)")]
    fn characteristic_polynomial(
        &self,
        py: Python<'_>,
        polynomial: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyUnivariatePolynomial> {
        let inner = match (&self.inner, &polynomial.inner) {
            (AnyQuotient::Prime(quotient), AnyPolynomial::Prime(value)) => AnyUnivariate::Prime(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.characteristic_polynomial(value, budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            (AnyQuotient::Rational(quotient), AnyPolynomial::Rational(value)) => {
                AnyUnivariate::Rational(
                    run_budgeted(py, timeout, memory_limit, move |budget| {
                        quotient.characteristic_polynomial(value, budget)
                    })?
                    .map_err(|error| of_quotient(py, error))?,
                )
            }
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyUnivariatePolynomial { inner })
    }

    /// Return the minimal polynomial of multiplication by a residue class.
    #[pyo3(signature = (polynomial, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, polynomial, *, timeout=None, memory_limit=None)")]
    fn minimal_polynomial(
        &self,
        py: Python<'_>,
        polynomial: &PyPolynomial,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyUnivariatePolynomial> {
        let inner = match (&self.inner, &polynomial.inner) {
            (AnyQuotient::Prime(quotient), AnyPolynomial::Prime(value)) => AnyUnivariate::Prime(
                run_budgeted(py, timeout, memory_limit, move |budget| {
                    quotient.minimal_polynomial(value, budget)
                })?
                .map_err(|error| of_quotient(py, error))?,
            ),
            (AnyQuotient::Rational(quotient), AnyPolynomial::Rational(value)) => {
                AnyUnivariate::Rational(
                    run_budgeted(py, timeout, memory_limit, move |budget| {
                        quotient.minimal_polynomial(value, budget)
                    })?
                    .map_err(|error| of_quotient(py, error))?,
                )
            }
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyUnivariatePolynomial { inner })
    }

    fn __repr__(&self) -> String {
        format!("FiniteQuotient(dimension={})", self.dimension())
    }
}

/// A square matrix of multiplication by a quotient residue class.
#[pyclass(frozen, name = "MultiplicationMatrix", module = "sylvester")]
struct PyMultiplicationMatrix {
    inner: AnyMatrix,
}

#[pymethods]
impl PyMultiplicationMatrix {
    /// The coefficient ring used by the matrix.
    #[getter]
    fn ring(&self) -> PyRing {
        PyRing {
            inner: ring_of_matrix(&self.inner),
        }
    }

    /// The number of rows and columns.
    #[getter]
    fn dimension(&self) -> usize {
        match &self.inner {
            AnyMatrix::Prime(value) => value.dimension(),
            AnyMatrix::Rational(value) => value.dimension(),
        }
    }

    /// The entries as rows of a Python list. Zero entries are ordinary zero
    /// coefficients in the ring's Python domain.
    #[getter]
    fn entries<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let dimension = self.dimension();
        match &self.inner {
            AnyMatrix::Prime(value) => PyList::new(
                py,
                (0..dimension)
                    .map(|row| {
                        (0..dimension)
                            .map(|column| {
                                value
                                    .entry(row, column)
                                    .map_or(0, |coefficient| coefficient.value())
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            ),
            AnyMatrix::Rational(value) => PyList::new(
                py,
                (0..dimension)
                    .map(|row| {
                        (0..dimension)
                            .map(|column| {
                                value
                                    .entry(row, column)
                                    .cloned()
                                    .unwrap_or_else(BigRational::zero)
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            ),
        }
    }

    /// The entries in row-major order.
    #[getter]
    fn flat_entries<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        match &self.inner {
            AnyMatrix::Prime(value) => prime_optional_coefficients(
                py,
                value
                    .entries()
                    .iter()
                    .map(|coefficient| coefficient.as_ref().map(|value| value.value())),
            ),
            AnyMatrix::Rational(value) => {
                rational_optional_coefficients(py, value.entries().iter().map(Clone::clone))
            }
        }
    }

    /// Read one matrix entry, or None for an out-of-range index.
    /// Negative and too large integer indexes are out of range.
    fn entry<'py>(
        &self,
        py: Python<'py>,
        row: &Bound<'_, PyAny>,
        column: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAny>>> {
        let Some(row) = index_from(row)? else {
            return Ok(None);
        };
        let Some(column) = index_from(column)? else {
            return Ok(None);
        };
        if row >= self.dimension() || column >= self.dimension() {
            return Ok(None);
        }
        match &self.inner {
            AnyMatrix::Prime(value) => Ok(Some(prime_coefficient_object(
                py,
                value.entry(row, column).unwrap_or_else(|| prime_zero()),
            )?)),
            AnyMatrix::Rational(value) => Ok(Some(rational_coefficient_object(
                py,
                value.entry(row, column).unwrap_or_else(|| rational_zero()),
            )?)),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "MultiplicationMatrix({}x{})",
            self.dimension(),
            self.dimension()
        )
    }
}

/// A univariate polynomial returned by a quotient spectrum computation.
#[pyclass(frozen, name = "UnivariatePolynomial", module = "sylvester")]
struct PyUnivariatePolynomial {
    inner: AnyUnivariate,
}

#[pymethods]
impl PyUnivariatePolynomial {
    /// The ring that supplies the coefficients.
    #[getter]
    fn ring(&self) -> PyRing {
        PyRing {
            inner: ring_of_univariate(&self.inner),
        }
    }

    /// Coefficients from the constant term upward. Zero entries are ordinary
    /// zero coefficients in the ring's Python domain.
    #[getter]
    fn coefficients<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        match &self.inner {
            AnyUnivariate::Prime(value) => prime_optional_coefficients(
                py,
                value
                    .coefficients()
                    .iter()
                    .map(|coefficient| coefficient.as_ref().map(|value| value.value())),
            ),
            AnyUnivariate::Rational(value) => {
                rational_optional_coefficients(py, value.coefficients().iter().map(Clone::clone))
            }
        }
    }

    /// Return the coefficient of `t**degree`, or None past the coefficient list.
    /// A negative or too large integer degree is out of range.
    fn coefficient<'py>(
        &self,
        py: Python<'py>,
        degree: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyAny>>> {
        let Some(degree) = index_from(degree)? else {
            return Ok(None);
        };
        match &self.inner {
            AnyUnivariate::Prime(value) => {
                if degree >= value.coefficients().len() {
                    return Ok(None);
                }
                Ok(Some(prime_coefficient_object(
                    py,
                    value.coefficient(degree).unwrap_or_else(|| prime_zero()),
                )?))
            }
            AnyUnivariate::Rational(value) => {
                if degree >= value.coefficients().len() {
                    return Ok(None);
                }
                Ok(Some(rational_coefficient_object(
                    py,
                    value.coefficient(degree).unwrap_or_else(|| rational_zero()),
                )?))
            }
        }
    }

    /// The largest nonzero degree, or None for the zero polynomial.
    #[getter]
    fn degree(&self) -> Option<usize> {
        match &self.inner {
            AnyUnivariate::Prime(value) => value.degree(),
            AnyUnivariate::Rational(value) => value.degree(),
        }
    }

    /// Report whether every coefficient is zero.
    #[getter]
    fn is_zero(&self) -> bool {
        match &self.inner {
            AnyUnivariate::Prime(value) => value.is_zero(),
            AnyUnivariate::Rational(value) => value.is_zero(),
        }
    }

    fn __str__(&self) -> String {
        match &self.inner {
            AnyUnivariate::Prime(value) => value.to_string(),
            AnyUnivariate::Rational(value) => value.to_string(),
        }
    }

    fn __repr__(&self) -> String {
        format!("UnivariatePolynomial({:?})", self.__str__())
    }
}

/// A saved input, basis, and untrusted provenance claim.
#[pyclass(frozen, name = "ResultEnvelope", module = "sylvester")]
struct PyResultEnvelope {
    inner: ResultEnvelope,
}

fn provenance_name(provenance: ClaimedProvenance) -> &'static str {
    match provenance {
        ClaimedProvenance::Unverified => "unverified",
        ClaimedProvenance::SuppliedBasis => "supplied_basis",
        ClaimedProvenance::Unchanged => "unchanged",
        ClaimedProvenance::ContainsInput => "contains_input",
        ClaimedProvenance::EqualsInput => "equals_input",
        ClaimedProvenance::Certified => "certified",
    }
}

#[pymethods]
impl PyResultEnvelope {
    /// Record a prime-field ideal and basis with an optional unverified certificate.
    #[staticmethod]
    #[pyo3(signature = (ideal, basis, certificate=None, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "(ideal, basis, certificate=None, *, timeout=None, memory_limit=None)")]
    fn from_prime(
        py: Python<'_>,
        ideal: &PyIdeal,
        basis: &PyBasis,
        certificate: Option<Vec<u8>>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Self> {
        let cancellation = CancellationToken::new();
        let budget = cancellable_budget(timeout, memory_limit, &cancellation)?;
        let inner = match (&ideal.inner, &basis.inner) {
            (AnyIdeal::Prime(ideal), AnyBasis::Prime(basis)) => run_with_signals(
                py,
                || cancellation.cancel(),
                move || ResultEnvelope::from_prime(ideal, basis, certificate.as_deref(), budget),
            )?
            .map_err(|error| of_envelope(py, error))?,
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyResultEnvelope { inner })
    }

    /// Record a certified prime-field basis after matching its certificate.
    #[staticmethod]
    #[pyo3(signature = (ideal, certified, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "(ideal, certified, *, timeout=None, memory_limit=None)")]
    fn from_certified(
        py: Python<'_>,
        ideal: &PyIdeal,
        certified: &PyCertified,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Self> {
        let cancellation = CancellationToken::new();
        let budget = cancellable_budget(timeout, memory_limit, &cancellation)?;
        let inner = match &ideal.inner {
            AnyIdeal::Prime(ideal) => {
                let certified = &certified.inner;
                run_with_signals(
                    py,
                    || cancellation.cancel(),
                    move || ResultEnvelope::from_certified(ideal, certified, budget),
                )?
                .map_err(|error| of_envelope(py, error))?
            }
            AnyIdeal::Rational(_) => {
                return Err(PyValueError::new_err(
                    "the rational numbers have no certified path".to_string(),
                ));
            }
        };
        Ok(PyResultEnvelope { inner })
    }

    /// Record a rational ideal and basis with its lift claim.
    #[staticmethod]
    #[pyo3(signature = (ideal, basis, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "(ideal, basis, *, timeout=None, memory_limit=None)")]
    fn from_rational(
        py: Python<'_>,
        ideal: &PyIdeal,
        basis: &PyBasis,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Self> {
        let cancellation = CancellationToken::new();
        let budget = cancellable_budget(timeout, memory_limit, &cancellation)?;
        let inner = match (&ideal.inner, &basis.inner) {
            (AnyIdeal::Rational(ideal), AnyBasis::Rational(basis)) => run_with_signals(
                py,
                || cancellation.cancel(),
                move || ResultEnvelope::from_rational(ideal, basis, budget),
            )?
            .map_err(|error| of_envelope(py, error))?,
            _ => return Err(foreign_argument(py)),
        };
        Ok(PyResultEnvelope { inner })
    }

    /// Record a successful exact rational equality check.
    #[staticmethod]
    #[pyo3(signature = (checked, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "(checked, *, timeout=None, memory_limit=None)")]
    fn from_checked_rational(
        py: Python<'_>,
        checked: &PyRationalEqualityCheck,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Self> {
        let checked = &checked.inner;
        let inner = run_budgeted(py, timeout, memory_limit, move |budget| {
            ResultEnvelope::from_checked_rational(checked, budget)
        })?
        .map_err(|error| of_envelope(py, error))?;
        Ok(PyResultEnvelope { inner })
    }

    /// Recheck a loaded rational record for exact ideal equality.
    #[pyo3(signature = (*, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, *, timeout=None, memory_limit=None)")]
    fn check_rational(
        &self,
        py: Python<'_>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<PyRationalEqualityCheck> {
        let record = &self.inner;
        let checked = run_budgeted(py, timeout, memory_limit, move |budget| {
            record.check_rational(budget)
        })?
        .map_err(|error| of_envelope(py, error))?;
        Ok(PyRationalEqualityCheck { inner: checked })
    }

    /// Load a record without trusting its provenance or certificate.
    #[staticmethod]
    #[pyo3(signature = (data, *, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "(data, *, timeout=None, memory_limit=None)")]
    fn from_json(
        py: Python<'_>,
        data: &Bound<'_, PyAny>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Self> {
        let cancellation = CancellationToken::new();
        let budget = cancellable_budget(timeout, memory_limit, &cancellation)?;
        let data = envelope_bytes(py, data, memory_limit)?;
        let inner = run_with_signals(
            py,
            || cancellation.cancel(),
            move || ResultEnvelope::from_json(&data, budget),
        )?
        .map_err(|error| of_envelope(py, error))?;
        Ok(PyResultEnvelope { inner })
    }

    /// Serialize the record to JSON bytes.
    #[pyo3(signature = (*, timeout=None, memory_limit=None))]
    #[pyo3(text_signature = "($self, *, timeout=None, memory_limit=None)")]
    fn to_json<'py>(
        &self,
        py: Python<'py>,
        timeout: Option<f64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Bound<'py, PyBytes>> {
        let cancellation = CancellationToken::new();
        let budget = cancellable_budget(timeout, memory_limit, &cancellation)?;
        let record = &self.inner;
        let data = run_with_signals(py, || cancellation.cancel(), move || record.to_json(budget))?
            .map_err(|error| of_envelope(py, error))?;
        Ok(PyBytes::new(py, &data))
    }

    /// The domain metadata as a dictionary.
    #[getter]
    fn domain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let out = PyDict::new(py);
        match self.inner.domain() {
            EnvelopeDomain::PrimeField { modulus } => {
                out.set_item("kind", "prime_field")?;
                out.set_item("modulus", modulus)?;
            }
            EnvelopeDomain::Rationals => {
                out.set_item("kind", "rationals")?;
            }
        }
        Ok(out)
    }

    /// The variable names in ring order.
    #[getter]
    fn variables(&self) -> Vec<String> {
        self.inner.variables().to_vec()
    }

    /// The input expressions stored in the record.
    #[getter]
    fn input(&self) -> Vec<String> {
        self.inner.input().to_vec()
    }

    /// The basis expressions stored in the record.
    #[getter]
    fn basis(&self) -> Vec<String> {
        self.inner.basis().to_vec()
    }

    /// The producer's untrusted provenance claim.
    #[getter]
    fn claimed_provenance(&self) -> &'static str {
        provenance_name(self.inner.claimed_provenance())
    }

    /// Certificate bytes stored in the record, if any.
    #[getter]
    fn certificate<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyBytes>> {
        self.inner
            .certificate()
            .map(|bytes| PyBytes::new(py, bytes))
    }

    /// Verify and match the stored prime certificate.
    #[pyo3(signature = (*, timeout=None, max_bytes=None, max_work_units=None, max_intermediate_bytes=None, max_live_bytes=None))]
    #[pyo3(
        text_signature = "($self, *, timeout=None, max_bytes=None, max_work_units=None, max_intermediate_bytes=None, max_live_bytes=None)"
    )]
    fn verify_prime(
        &self,
        py: Python<'_>,
        timeout: Option<f64>,
        max_bytes: Option<usize>,
        max_work_units: Option<u64>,
        max_intermediate_bytes: Option<usize>,
        max_live_bytes: Option<usize>,
    ) -> PyResult<PyVerified> {
        let mut caps = verification_limits(
            None,
            max_bytes,
            max_work_units,
            max_intermediate_bytes,
            max_live_bytes,
        )?;
        let cancellation = CancellationToken::new();
        caps.cancellation = Some(cancellation.shared_flag());
        let budget = cancellable_budget(timeout, None, &cancellation)?;
        let record = &self.inner;
        let verified = run_with_signals(
            py,
            || cancellation.cancel(),
            move || record.verify_prime(&caps, budget),
        )?
        .map_err(|error| of_envelope(py, error))?;
        verified_to_python(py, verified)
    }

    fn __repr__(&self) -> String {
        format!(
            "ResultEnvelope(domain={:?}, {} input, {} basis)",
            self.domain_name(),
            self.inner.input().len(),
            self.inner.basis().len()
        )
    }
}

impl PyResultEnvelope {
    fn domain_name(&self) -> &'static str {
        match self.inner.domain() {
            EnvelopeDomain::PrimeField { .. } => "prime_field",
            EnvelopeDomain::Rationals => "rationals",
        }
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

fn set_matrix_counters(dict: &Bound<'_, PyDict>, counters: &sylvester::F4Counters) -> PyResult<()> {
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

fn set_pair_counters(dict: &Bound<'_, PyDict>, counters: &sylvester::F4Counters) -> PyResult<()> {
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

/// Build verifier limits from Python keyword values.
fn verification_limits(
    timeout: Option<f64>,
    max_bytes: Option<usize>,
    max_work_units: Option<u64>,
    max_intermediate_bytes: Option<usize>,
    max_live_bytes: Option<usize>,
) -> PyResult<verify::Limits> {
    let mut limits = verify::Limits::default();
    if let Some(seconds) = timeout {
        limits.deadline = Instant::now().checked_add(duration_of("timeout", seconds)?);
    }
    if let Some(value) = max_bytes {
        limits.max_bytes = value;
    }
    if let Some(value) = max_work_units {
        limits.max_work_units = value;
    }
    if let Some(value) = max_intermediate_bytes {
        limits.max_intermediate_bytes = value;
    }
    if let Some(value) = max_live_bytes {
        limits.max_live_bytes = value;
    }
    Ok(limits)
}

fn certificate_bytes(
    py: Python<'_>,
    data: &Bound<'_, PyAny>,
    max_bytes: usize,
) -> PyResult<Vec<u8>> {
    let length = data.len()?;
    if length > max_bytes {
        return Err(raise_with(
            py,
            &exceptions(py).limit_exceeded,
            format!("the certificate is larger than the {max_bytes}-byte cap"),
            |value| value.setattr("limit", max_bytes),
        ));
    }
    let bytes = data.extract::<Vec<u8>>()?;
    if bytes.len() > max_bytes {
        return Err(raise_with(
            py,
            &exceptions(py).limit_exceeded,
            format!("the certificate is larger than the {max_bytes}-byte cap"),
            |value| value.setattr("limit", max_bytes),
        ));
    }
    Ok(bytes)
}

fn envelope_bytes(
    py: Python<'_>,
    data: &Bound<'_, PyAny>,
    memory_limit: Option<usize>,
) -> PyResult<Vec<u8>> {
    let length = data.len()?;
    if let Some(limit) = memory_limit {
        let estimate = length
            .saturating_mul(32)
            .saturating_add(std::mem::size_of::<ResultEnvelope>());
        if estimate > limit {
            return Err(raise(
                py,
                &exceptions(py).memory_limit,
                "the computation record passed its memory limit before decoding".to_string(),
            ));
        }
    }
    let bytes = data.extract::<Vec<u8>>()?;
    if let Some(limit) = memory_limit {
        let estimate = bytes
            .len()
            .saturating_mul(32)
            .saturating_add(std::mem::size_of::<ResultEnvelope>());
        if estimate > limit {
            return Err(raise(
                py,
                &exceptions(py).memory_limit,
                "the computation record passed its memory limit before decoding".to_string(),
            ));
        }
    }
    Ok(bytes)
}

fn verified_to_python(py: Python<'_>, verified: verify::VerifiedGb) -> PyResult<PyVerified> {
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

/// Check certificate bytes with the independent verifier.
///
/// `timeout` is seconds. The other keywords replace the matching default
/// resource caps. Exhaustion says nothing about the certificate.
#[pyfunction(name = "verify")]
#[pyo3(signature = (data, *, timeout=None, max_bytes=None, max_work_units=None, max_intermediate_bytes=None, max_live_bytes=None))]
#[pyo3(
    text_signature = "(data, *, timeout=None, max_bytes=None, max_work_units=None, max_intermediate_bytes=None, max_live_bytes=None)"
)]
fn verify_certificate(
    py: Python<'_>,
    data: &Bound<'_, PyAny>,
    timeout: Option<f64>,
    max_bytes: Option<usize>,
    max_work_units: Option<u64>,
    max_intermediate_bytes: Option<usize>,
    max_live_bytes: Option<usize>,
) -> PyResult<PyVerified> {
    let mut limits = verification_limits(
        timeout,
        max_bytes,
        max_work_units,
        max_intermediate_bytes,
        max_live_bytes,
    )?;
    let data = certificate_bytes(py, data, limits.max_bytes)?;
    let cancellation = CancellationToken::new();
    limits.cancellation = Some(cancellation.shared_flag());
    let verified = run_with_signals(
        py,
        || cancellation.cancel(),
        move || verify::verify_with_limits(&data, &limits),
    )?
    .map_err(|error| of_verify(py, error))?;
    verified_to_python(py, verified)
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
    add_ring_classes(m)?;
    add_quotient_classes(m)?;
    add_result_classes(m)?;
    Ok(())
}

fn add_ring_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRing>()?;
    m.add_class::<PyPolynomial>()?;
    m.add_class::<PyIdeal>()?;
    m.add_class::<PyBasis>()?;
    Ok(())
}

fn add_quotient_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRationalEqualityCheck>()?;
    m.add_class::<PyFiniteQuotient>()?;
    m.add_class::<PyMultiplicationMatrix>()?;
    m.add_class::<PyUnivariatePolynomial>()?;
    Ok(())
}

fn add_result_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyResultEnvelope>()?;
    m.add_class::<PyCertified>()?;
    m.add_class::<PyHilbertSeries>()?;
    m.add_class::<PyComputeReport>()?;
    m.add_class::<PyVerified>()?;
    Ok(())
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
