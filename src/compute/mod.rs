//! Gröbner basis engines and the options that drive them.

mod classic;
mod interreduce;
mod matrix;
mod signature;

use std::fmt;
use std::time::{Duration, Instant};

use crate::cert::{self, Budget};
use crate::certificate::{CertifiedGroebnerBasis, CertifyError};
use crate::ideal::GroebnerBasis;
use crate::poly::Polynomial;
use crate::ring::PolynomialRing;
use crate::verify::{self, Limits};

/// The engine that computes a basis.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Backend {
    /// Batch the critical pairs by degree and reduce them in sparse
    /// Macaulay matrices with F4-style elimination, under the F5 syzygy
    /// criterion.
    #[default]
    Matrix,
    /// Process the critical pairs one at a time in signature order,
    /// as classic F5 does.
    Classic,
}

/// The budget and the backend one computation runs under.
///
/// ```
/// use std::time::Duration;
/// use sylvester::{Backend, ComputeOptions};
///
/// let options = ComputeOptions::new()
///     .backend(Backend::Classic)
///     .timeout(Duration::from_secs(30))
///     .memory_limit(512 << 20);
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ComputeOptions {
    backend: Backend,
    timeout: Option<Duration>,
    memory_limit: Option<usize>,
}

impl ComputeOptions {
    /// The default options: the matrix backend, no deadline, no memory
    /// limit.
    pub fn new() -> Self {
        ComputeOptions::default()
    }

    /// Pick the backend.
    pub fn backend(mut self, backend: Backend) -> Self {
        self.backend = backend;
        self
    }

    /// Stop the computation after this much time.
    ///
    /// The deadline starts when the computation starts. An exhausted
    /// deadline is [`ComputeError::Timeout`]. A `timeout` too large to add
    /// to the start instant leaves the computation with no deadline,
    /// instead of stopping it before it starts.
    pub fn timeout(mut self, timeout: Duration) -> Self {
        self.timeout = Some(timeout);
        self
    }

    /// Stop the computation once the live engine data passes this many
    /// bytes.
    ///
    /// The count covers the basis, the pair queue, and the tracked
    /// cofactors. On the certified path it also covers the certificate
    /// emitter and sizes the verifier caps. It is not a cap on the process.
    /// An exhausted budget is [`ComputeError::MemoryLimitExceeded`].
    pub fn memory_limit(mut self, bytes: usize) -> Self {
        self.memory_limit = Some(bytes);
        self
    }

    fn deadline(&self) -> Option<Instant> {
        self.timeout
            .and_then(|timeout| Instant::now().checked_add(timeout))
    }
}

/// The width of one exponent, and so the largest total degree an input
/// generator or a critical pair may reach.
///
/// A monomial holds each exponent in a `u16`. Both backends reject a
/// generator above this degree at entry and a pair above it at creation.
/// Those two gates bound every product the polynomial path builds, by the
/// argument stated at the entry gate in [`classic`]. A signature or a
/// cofactor product carries no such bound and reports the limit where it
/// is formed.
pub(crate) const DEGREE_LIMIT: u32 = u16::MAX as u32;

/// Term operations between two deadline polls inside one reduction.
///
/// `Instant::now()` costs more than one subtraction, so a reduction polls
/// the clock once per this many term operations, off the inner arithmetic.
/// The stride bounds one long reduction's overrun to milliseconds on
/// recorded workloads.
pub(super) const REDUCE_DEADLINE_STRIDE: usize = 4096;

/// Return [`ComputeError::Timeout`] once the deadline has passed.
#[inline]
pub(super) fn poll_deadline(deadline: Option<Instant>) -> Result<(), ComputeError> {
    match deadline {
        Some(deadline) if Instant::now() >= deadline => Err(ComputeError::Timeout),
        _ => Ok(()),
    }
}

/// The bytes one polynomial term occupies over a ring of `nvars` variables:
/// the term itself, plus what its exponent vector spills to the heap.
pub(super) fn per_term_bytes(nvars: usize) -> usize {
    size_of::<crate::poly::Term>() + crate::poly::heap_exps_bytes(nvars)
}

/// Why a computation stops before it has a basis.
///
/// [`ComputeError::Timeout`] and [`ComputeError::MemoryLimitExceeded`]
/// report an exhausted budget; neither says anything about the ideal, and
/// the computation stopped rather than failed. [`ComputeError::DegreeLimit`]
/// reports an input past what this release supports: an input generator of
/// total degree above 65535, a critical pair whose leading monomials have a
/// least common multiple above that degree, or a signature or cofactor
/// product past the width one exponent holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComputeError {
    /// The deadline passed.
    Timeout,
    /// The live engine data passed the memory limit.
    MemoryLimitExceeded,
    /// A monomial the computation needs is past the degree an exponent's
    /// width supports: an input generator, a critical pair's least common
    /// multiple, or a signature or cofactor product.
    DegreeLimit {
        /// The largest degree a monomial may reach.
        limit: u32,
    },
}

impl fmt::Display for ComputeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ComputeError::Timeout => f.write_str("the computation passed its deadline"),
            ComputeError::MemoryLimitExceeded => {
                f.write_str("the computation passed its memory limit")
            }
            ComputeError::DegreeLimit { limit } => write!(
                f,
                "the computation reached a monomial of degree above {limit}, which an exponent's width does not support"
            ),
        }
    }
}

impl std::error::Error for ComputeError {}

impl From<crate::poly::ExponentOverflow> for ComputeError {
    fn from(_: crate::poly::ExponentOverflow) -> Self {
        ComputeError::DegreeLimit {
            limit: DEGREE_LIMIT,
        }
    }
}

/// Reject a generator of total degree above [`DEGREE_LIMIT`].
///
/// Grevlex is graded, so a polynomial's leading monomial carries its
/// largest degree. With every generator at or below the limit, and every
/// critical pair's least common multiple checked at creation, every
/// product the polynomial path builds stays inside the width one exponent
/// holds: a tail term's degree is at most its lead's, an S-pair multiplier
/// m satisfies deg(m * t) <= deg(lcm), and a reduction multiplier q
/// satisfies deg(q * t) <= deg(lm(p)). Signature and cofactor products
/// carry no such bound; they check the width where they are formed and
/// report the same limit.
pub(crate) fn check_input_degrees(generators: &[Polynomial]) -> Result<(), ComputeError> {
    let past_limit = generators
        .iter()
        .any(|f| f.degree().is_some_and(|deg| deg > DEGREE_LIMIT));
    if past_limit {
        return Err(ComputeError::DegreeLimit {
            limit: DEGREE_LIMIT,
        });
    }
    Ok(())
}

pub(crate) fn groebner_basis(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &ComputeOptions,
) -> Result<Vec<Polynomial>, ComputeError> {
    let deadline = options.deadline();
    let budget = options.memory_limit;
    match options.backend {
        Backend::Matrix => matrix::solve_checked(ring, generators, deadline, budget),
        Backend::Classic => classic::solve_checked(ring, generators, deadline, budget),
    }
}

/// Run the classic backend, write its certificate, and verify it.
///
/// `crate::cert` only assembles the certificate bytes; it never imports
/// [`crate::verify`]. This function is the one place that reads the
/// emitter's output back through the independent verifier, so a defect in
/// either module cannot hide behind a shared assumption.
///
/// `deadline` and `memory_limit` are the budget of the whole run. They
/// cover the engine, the emitter, and the verifier. An exhausted budget is
/// [`CertifyError::Engine`] before the verifier starts and
/// [`CertifyError::VerifierExhausted`] inside it. Neither says anything
/// about the basis.
pub(crate) fn groebner_basis_certified(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &ComputeOptions,
) -> Result<CertifiedGroebnerBasis, CertifyError> {
    let deadline = options.deadline();
    let memory_limit = options.memory_limit;

    let (basis, origins) = classic::solve_tracked(ring, generators, deadline, memory_limit)?;

    let mut budget = Budget::new(deadline, memory_limit, ring.nvars());
    let certificate = cert::assemble(ring, generators, basis, origins, &mut budget)?;
    cert::check_deadline(deadline)?;

    let limits = verifier_limits(
        certificate.len(),
        budget.held(),
        ring.nvars(),
        deadline,
        memory_limit,
    );
    let verified = verify::verify_with_limits(&certificate, &limits)?;

    // The verifier proves the certificate consistent with the input the
    // certificate itself carries. It cannot know which ideal the caller
    // asked about, so the caller checks that here.
    if verified.modulus() != ring.modulus() || verified.nvars() != ring.nvars() {
        return Err(CertifyError::InputMismatch);
    }
    let decoded_input: Vec<Polynomial> = verified
        .input()
        .iter()
        .map(|poly| decode(ring, poly))
        .collect();
    if decoded_input != generators {
        return Err(CertifyError::InputMismatch);
    }

    let polynomials = verified
        .basis()
        .iter()
        .map(|poly| decode(ring, poly))
        .collect();
    Ok(CertifiedGroebnerBasis::new(
        GroebnerBasis::new(ring.clone(), polynomials),
        certificate,
    ))
}

/// Map the budget of a run to the caps of one verification.
///
/// `certificate_len` is the length of the certificate bytes, and `held` is
/// what the emitter's budget holds after charging them, so `held` is at
/// least `certificate_len`. What the memory limit leaves once `held` is
/// paid for pays for the decoded certificate and for the results of the
/// identity checks; `max_bytes` is that allowance, which by construction
/// covers the certificate the caller is about to hand the verifier. A
/// check holds the decoded certificate and at most two results at once, so
/// each term cap gets a third of what is left. The map is deliberately
/// conservative: a limit too small for the whole run stops the verifier
/// with [`CertifyError::VerifierExhausted`], which is not a rejection.
///
/// Without a memory limit the caps are the verifier defaults, raised to the
/// length of the certificate where a default sits below it. The defaults
/// hold for bytes from an untrusted peer. These bytes come from the emitter
/// in this process, and one term of them costs at least one byte, so the
/// length bounds every count the caps name.
fn verifier_limits(
    certificate_len: usize,
    held: usize,
    nvars: usize,
    deadline: Option<Instant>,
    memory_limit: Option<usize>,
) -> Limits {
    let defaults = Limits::default();
    let Some(limit) = memory_limit else {
        let len = certificate_len;
        return Limits {
            max_bytes: defaults.max_bytes.max(len),
            max_polys: defaults.max_polys.max(len),
            max_terms_per_poly: defaults.max_terms_per_poly.max(len),
            max_total_terms: defaults.max_total_terms.max(len),
            max_entries: defaults.max_entries.max(len),
            max_intermediate_bytes: defaults.max_intermediate_bytes.max(len),
            deadline,
        };
    };

    let left = limit.saturating_sub(held) / 3;
    let term_bytes = size_of::<verify::Term>() + nvars * size_of::<verify::Exp>();
    // The verifier arithmetic cap counts bytes, not terms, so it takes
    // what is left directly. It never falls below one term.
    let terms = (left / term_bytes).max(1);
    let polys = (left / size_of::<verify::Poly>()).max(1);
    Limits {
        max_bytes: limit.saturating_sub(held) + certificate_len,
        max_polys: polys,
        max_terms_per_poly: terms,
        max_total_terms: terms,
        max_entries: polys,
        max_intermediate_bytes: left.max(term_bytes),
        deadline,
    }
}

/// Read one verified polynomial back into the ring.
///
/// The verifier accepted bytes the emitter wrote from ring polynomials, so
/// every coefficient and every exponent came from this ring.
fn decode(ring: &PolynomialRing, poly: &verify::Poly) -> Polynomial {
    let terms = poly.terms().iter().map(|term| {
        let exps: Vec<u16> = term
            .mono()
            .exps()
            .iter()
            // the emitter wrote exponents that came from u16 values.
            .map(|&exp| u16::try_from(exp).expect("a certificate exponent fits a u16"))
            .collect();
        (term.coeff() as i64, exps)
    });
    ring.polynomial(terms)
        // the emitter wrote one exponent per variable of this ring.
        .expect("a certificate polynomial belongs to the ring that produced it")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_timeout_too_large_to_add_to_now_leaves_no_deadline() {
        let options = ComputeOptions::new().timeout(Duration::MAX);
        assert_eq!(options.deadline(), None);
    }

    #[test]
    fn an_ordinary_timeout_still_sets_a_deadline() {
        let options = ComputeOptions::new().timeout(Duration::from_secs(60));
        assert!(options.deadline().is_some());
    }

    #[test]
    fn the_verifier_limits_come_from_what_is_left_of_the_budget() {
        let bytes = vec![0u8; 1000];
        let unlimited = verifier_limits(bytes.len(), bytes.len(), 2, None, None);
        assert_eq!(unlimited.max_bytes, Limits::default().max_bytes);
        assert_eq!(unlimited.max_total_terms, Limits::default().max_total_terms);

        // A certificate longer than a default cap raises that cap. One term
        // costs at least one byte, so the length bounds every count.
        let long = vec![0u8; Limits::default().max_bytes + 1];
        let raised = verifier_limits(long.len(), long.len(), 2, None, None);
        assert_eq!(raised.max_bytes, long.len());
        assert_eq!(raised.max_total_terms, long.len());
        assert_eq!(raised.max_polys, long.len());

        // Nothing but the certificate itself is charged to the budget, so
        // the whole limit is left over for the byte cap and the term caps.
        let bounded = verifier_limits(bytes.len(), bytes.len(), 2, None, Some(1 << 20));
        assert_eq!(bounded.max_bytes, 1 << 20);
        assert!(bounded.max_total_terms < Limits::default().max_total_terms);
        assert!(bounded.max_total_terms > 0);

        // A limit the certificate alone fills leaves the smallest caps, and
        // no cap is zero.
        let spent = verifier_limits(bytes.len(), bytes.len(), 2, None, Some(bytes.len()));
        assert_eq!(spent.max_bytes, bytes.len());
        assert_eq!(spent.max_polys, 1);
        assert_eq!(spent.max_total_terms, 1);
        assert_eq!(
            spent.max_intermediate_bytes,
            size_of::<verify::Term>() + 2 * size_of::<verify::Exp>()
        );

        // Bytes the emitter already held before it wrote the certificate
        // count against the limit too, so the term caps shrink with them
        // even though the certificate itself did not grow.
        let held = bytes.len() + 500;
        let with_other_data = verifier_limits(bytes.len(), held, 2, None, Some(1 << 20));
        assert!(with_other_data.max_total_terms < bounded.max_total_terms);
    }
}
