//! Gröbner basis engines and the options that drive them.

mod classic;
pub(crate) mod f4;
mod interreduce;
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
    /// Batch critical pairs by degree, build sparse matrices over interned
    /// monomials, and use Gebauer-Moller pair management. This is the
    /// default.
    #[default]
    F4,
    /// Process the critical pairs one at a time in signature order,
    /// as classic F5 does.
    Classic,
}

/// The backend, resource budget, and thread count of one computation.
///
/// ```
/// use std::time::Duration;
/// use sylvester::{Backend, ComputeOptions};
///
/// let options = ComputeOptions::new()
///     .backend(Backend::Classic)
///     .timeout(Duration::from_secs(30))
///     .memory_limit(512 << 20)
///     .threads(1);
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ComputeOptions {
    backend: Backend,
    timeout: Option<Duration>,
    memory_limit: Option<usize>,
    threads: Option<usize>,
}

impl ComputeOptions {
    /// The default options: F4, no deadline, no memory limit, and the
    /// rayon global pool's thread count.
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
    /// to the start instant leaves the computation with no deadline.
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

    /// Compute on this many threads.
    ///
    /// Zero selects the rayon global pool, as the default does. One keeps
    /// the whole computation on the calling thread. A larger count builds
    /// a pool for this F4 run. The classic backend and every certified run
    /// stay on one thread. The count changes no basis or certificate
    /// bytes.
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = (threads > 0).then_some(threads);
        self
    }

    fn deadline(&self) -> Option<Instant> {
        self.timeout
            .and_then(|timeout| Instant::now().checked_add(timeout))
    }
}

/// What one computation did, next to the basis it produced.
///
/// [`crate::Ideal::groebner_basis_with_report`] returns this value for
/// benchmarks and diagnostics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ComputeReport {
    /// The backend that ran.
    pub backend: Backend,
    /// F4 counters, or `None` for the classic backend.
    pub counters: Option<F4Counters>,
    /// Wall time of the engine call, excluding certificate work.
    pub elapsed: Duration,
    /// The number of threads available to the run.
    pub threads_used: usize,
}

/// Counters reported by one F4 run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct F4Counters {
    /// Batches reduced, including final interreduction.
    pub batches: u64,
    /// Rows across all matrices.
    pub matrix_rows: u64,
    /// Columns across all matrices.
    pub matrix_columns: u64,
    /// Nonzero entries before reduction.
    pub matrix_nonzeros: u64,
    /// Rows that reduced to zero.
    pub zero_rows: u64,
    /// Pivots installed by reduction.
    pub new_pivots: u64,
    /// Restarts at a wider exponent packing.
    pub lane_restarts: u32,
    /// Batches retried with fewer pairs to meet the memory limit.
    pub batch_retries: u64,
    /// Monomials in the final basis table.
    pub basis_monomials: usize,
    /// Candidate critical pairs formed.
    pub pairs_generated: u64,
    /// Candidates discarded by the product criterion.
    pub pairs_discarded_product: u64,
    /// Queued pairs deleted by criterion B.
    pub pairs_discarded_b: u64,
    /// Candidates discarded by criterion M.
    pub pairs_discarded_m: u64,
    /// Candidates discarded by criterion F.
    pub pairs_discarded_f: u64,
}

impl From<f4::RunCounters> for F4Counters {
    fn from(counters: f4::RunCounters) -> Self {
        F4Counters {
            batches: counters.batches,
            matrix_rows: counters.matrix_rows,
            matrix_columns: counters.matrix_columns,
            matrix_nonzeros: counters.matrix_nonzeros,
            zero_rows: counters.zero_rows,
            new_pivots: counters.new_pivots,
            lane_restarts: counters.lane_restarts,
            batch_retries: counters.batch_retries,
            basis_monomials: counters.basis_monomials,
            pairs_generated: counters.pairs.generated,
            pairs_discarded_product: counters.pairs.product,
            pairs_discarded_b: counters.pairs.criterion_b,
            pairs_discarded_m: counters.pairs.criterion_m,
            pairs_discarded_f: counters.pairs.criterion_f,
        }
    }
}

/// The width of one exponent, and so the largest total degree an input
/// generator or a critical pair may reach.
///
/// A monomial holds each exponent in a `u16`. Both backends reject a
/// generator above this degree at entry and a pair above it at creation.
/// Those two gates bound every product the polynomial path builds, by the
/// argument at [`check_input_degrees`]. A signature or a cofactor product
/// carries no such bound and reports the limit where it is formed.
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
/// report an exhausted budget. Neither says anything about the ideal.
/// [`ComputeError::DegreeLimit`] reports a monomial whose total degree is
/// past the width of one exponent. [`ComputeError::ExponentLimit`] reports
/// an F4 monomial with one exponent past that width.
/// [`ComputeError::TableFull`] reports an F4 monomial table that cannot
/// intern another monomial.
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
    /// An F4 monomial needs one exponent above `limit`.
    ExponentLimit {
        /// The largest value one exponent holds.
        limit: u32,
    },
    /// An F4 monomial table reached its maximum number of entries.
    TableFull,
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
            ComputeError::ExponentLimit { limit } => write!(
                f,
                "a monomial reached an exponent above {limit}, which is the largest one exponent holds"
            ),
            ComputeError::TableFull => {
                f.write_str("a monomial table reached the largest number of monomials it holds")
            }
        }
    }
}

impl std::error::Error for ComputeError {}

impl From<crate::poly::ExponentOverflow> for ComputeError {
    fn from(_: crate::poly::ExponentOverflow) -> Self {
        // Polynomial-path overflow is DegreeLimit: a product past the
        // width of one exponent. F4 reports a single exponent past that
        // width as ExponentLimit.
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
        Backend::F4 => f4::groebner_basis(ring, generators, options),
        Backend::Classic => classic::solve_checked(ring, generators, deadline, budget),
    }
}

/// Run the selected backend and report what it did.
pub(crate) fn groebner_basis_with_report(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &ComputeOptions,
) -> Result<(Vec<Polynomial>, ComputeReport), ComputeError> {
    let deadline = options.deadline();
    let budget = options.memory_limit;
    let start = Instant::now();
    let (polynomials, counters, threads) = match options.backend {
        Backend::F4 => {
            let (basis, counters) = f4::solve(ring, generators, options)?;
            (basis, Some(F4Counters::from(counters)), counters.threads)
        }
        Backend::Classic => (
            classic::solve_checked(ring, generators, deadline, budget)?,
            None,
            Some(1),
        ),
    };
    let report = ComputeReport {
        backend: options.backend,
        counters,
        elapsed: start.elapsed(),
        threads_used: threads.unwrap_or_else(rayon::current_num_threads),
    };
    Ok((polynomials, report))
}

/// Run the selected backend, write its certificate, and verify it.
///
/// This function is the one place that reads the emitter's output back
/// through the independent verifier, so a defect in either module cannot
/// hide behind a shared assumption.
///
/// The deadline and the memory limit are the budget of the whole run. They
/// cover the engine, the emitter, and the verifier. An exhausted budget is
/// [`CertifyError::Engine`] before the verifier starts and
/// [`CertifyError::VerifierExhausted`] inside it. Neither says anything
/// about the basis.
pub(crate) fn groebner_basis_certified(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &ComputeOptions,
) -> Result<CertifiedGroebnerBasis, CertifyError> {
    let run_limits = f4::Limits::of(options);
    let deadline = run_limits.deadline;
    let memory_limit = run_limits.memory;

    let (certificate, held) = match options.backend {
        Backend::Classic => classic_certificate(ring, generators, deadline, memory_limit)?,
        Backend::F4 => f4_certificate(ring, generators, options, &run_limits)?,
    };
    cert::check_deadline(deadline)?;

    let limits = verifier_limits(
        certificate.len(),
        held,
        ring.nvars(),
        deadline,
        memory_limit,
    );
    let verified = verify::verify_with_limits(&certificate, &limits)?;

    // The verifier accepts a certificate as consistent with the input the
    // certificate itself carries. It cannot know which ideal the caller
    // asked about, so this function checks that here.
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

/// Run classic and write its v1 cofactor certificate.
fn classic_certificate(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    deadline: Option<Instant>,
    memory_limit: Option<usize>,
) -> Result<(Vec<u8>, usize), CertifyError> {
    let (basis, origins) = classic::solve_tracked(ring, generators, deadline, memory_limit)?;
    let mut budget = Budget::new(deadline, memory_limit, ring.nvars());
    let certificate = cert::assemble(ring, generators, basis, origins, &mut budget)?;
    Ok((certificate, budget.held()))
}

/// Run F4 with a recorder and write its v2 trace certificate.
fn f4_certificate(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &ComputeOptions,
    limits: &f4::Limits,
) -> Result<(Vec<u8>, usize), CertifyError> {
    let mut recorder = cert::v2::Recorder::new(ring, generators, limits.deadline, limits.memory);
    let (basis, _) = f4::solve_recorded(ring, generators, options, limits, &mut recorder)?;
    cert::v2::assemble(ring, generators, &basis, recorder)
}

/// Map the budget of a run to the caps of one verification.
///
/// `certificate_len` is the length of the certificate bytes. `held` is
/// what the emitter's budget holds after charging them, so `held` is at
/// least `certificate_len`. After `held` is charged, the remainder of the
/// memory limit pays for the decoded certificate and the identity-check
/// results. `max_bytes` is that remainder plus `certificate_len`. A check
/// holds the decoded certificate and at most two results at once, so each
/// term cap gets a third of the remainder. A limit too small for the
/// whole run stops the verifier with [`CertifyError::VerifierExhausted`].
/// That is not a rejection.
///
/// Without a memory limit the caps are the verifier defaults, raised to
/// the certificate length where a default sits below it. The defaults
/// hold for bytes from an untrusted peer. These bytes come from the
/// emitter in this process, and one term of them costs at least one
/// byte, so the length bounds every count the caps name.
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
            max_pool_monomials: defaults.max_pool_monomials.max(len),
            max_pool_entries: defaults.max_pool_entries.max(len),
            max_input_polys: defaults.max_input_polys.max(len),
            max_nodes: defaults.max_nodes.max(len),
            max_comb_steps: defaults.max_comb_steps.max(len),
            max_trace_steps: defaults.max_trace_steps.max(len),
            max_basis: defaults.max_basis.max(len),
            max_pairs: defaults.max_pairs.max(len),
            max_division_steps: defaults.max_division_steps.max(len),
            max_total_division_steps: defaults.max_total_division_steps.max(len),
            max_live_bytes: defaults.max_live_bytes.max(len),
            max_work_units: u64::MAX,
            deadline,
        };
    };

    let left = limit.saturating_sub(held) / 3;
    let term_bytes = size_of::<verify::Term>() + nvars * size_of::<verify::Exp>();
    let terms = (left / term_bytes).max(1);
    let polys = (left / size_of::<verify::Poly>()).max(1);
    let slots = (left / 64).max(1);
    Limits {
        max_bytes: limit.saturating_sub(held) + certificate_len,
        max_polys: polys,
        max_terms_per_poly: terms,
        max_total_terms: terms,
        max_entries: polys,
        // `max_intermediate_bytes` counts bytes, not terms, so it takes
        // `left` directly and never falls below one term.
        max_intermediate_bytes: left.max(term_bytes),
        max_pool_monomials: slots,
        max_pool_entries: slots,
        max_input_polys: polys,
        max_nodes: slots,
        max_comb_steps: slots,
        max_trace_steps: slots,
        max_basis: polys,
        max_pairs: slots,
        max_division_steps: slots,
        max_total_division_steps: slots,
        max_live_bytes: left.max(term_bytes),
        max_work_units: u64::MAX,
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
            // The emitter wrote exponents that came from u16 values.
            .map(|&exp| u16::try_from(exp).expect("a certificate exponent fits a u16"))
            .collect();
        (term.coeff() as i64, exps)
    });
    ring.polynomial(terms)
        // The emitter wrote one exponent per variable of this ring.
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
    fn a_thread_count_of_zero_means_the_default_pool() {
        assert_eq!(ComputeOptions::new().threads(0), ComputeOptions::new());
        assert_ne!(ComputeOptions::new().threads(1), ComputeOptions::new());
    }

    #[test]
    fn the_verifier_limits_come_from_what_is_left_of_the_budget() {
        let bytes = vec![0u8; 1000];
        check_unlimited_limits(&bytes);
        check_raised_limits();
        let bounded = check_bounded_limits(&bytes);
        check_spent_limits(&bytes);
        check_held_data_limits(&bytes, &bounded);
    }

    fn check_unlimited_limits(bytes: &[u8]) {
        let unlimited = verifier_limits(bytes.len(), bytes.len(), 2, None, None);
        assert_eq!(unlimited.max_bytes, Limits::default().max_bytes);
        assert_eq!(unlimited.max_total_terms, Limits::default().max_total_terms);
    }

    fn check_raised_limits() {
        let long = vec![0u8; Limits::default().max_bytes + 1];
        let raised = verifier_limits(long.len(), long.len(), 2, None, None);
        assert_eq!(raised.max_bytes, long.len());
        assert_eq!(raised.max_total_terms, long.len());
        assert_eq!(raised.max_polys, long.len());
    }

    fn check_bounded_limits(bytes: &[u8]) -> Limits {
        let bounded = verifier_limits(bytes.len(), bytes.len(), 2, None, Some(1 << 20));
        assert_eq!(bounded.max_bytes, 1 << 20);
        assert!(bounded.max_total_terms < Limits::default().max_total_terms);
        assert!(bounded.max_total_terms > 0);
        bounded
    }

    fn check_spent_limits(bytes: &[u8]) {
        let spent = verifier_limits(bytes.len(), bytes.len(), 2, None, Some(bytes.len()));
        assert_eq!(spent.max_bytes, bytes.len());
        assert_eq!(spent.max_polys, 1);
        assert_eq!(spent.max_total_terms, 1);
        assert_eq!(
            spent.max_intermediate_bytes,
            size_of::<verify::Term>() + 2 * size_of::<verify::Exp>()
        );
    }

    fn check_held_data_limits(bytes: &[u8], bounded: &Limits) {
        let held = bytes.len() + 500;
        let with_other_data = verifier_limits(bytes.len(), held, 2, None, Some(1 << 20));
        assert!(with_other_data.max_total_terms < bounded.max_total_terms);
    }
}
