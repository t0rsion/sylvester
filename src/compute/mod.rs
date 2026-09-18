//! Gröbner basis engines and the options that drive them.

mod classic;
pub(crate) mod f4;
mod interreduce;
mod modular;
mod signature;

use std::fmt;
use std::num::NonZeroUsize;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

use crate::cert::{self, WriterBudget};
use crate::certificate::{CertifiedGroebnerBasis, CertifyError};
use crate::ideal::GroebnerBasis;
use crate::poly::Polynomial;
use crate::ring::{ModularLift, PolynomialRing, Rationals};
use crate::verify::{self, Limits as VerifyLimits};

/// The engine that computes a basis.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Backend {
    /// Batch the critical pairs by degree, build one sparse matrix per
    /// batch over interned monomials, and reduce it. This is F4 with
    /// Gebauer-Moller pair management, and it is the default.
    #[default]
    F4,
    /// Process the critical pairs one at a time in signature order,
    /// as classic F5 does.
    Classic,
}

/// A one-way cancellation token shared by one or more computations.
///
/// Clones observe the same flag. [`CancellationToken::cancel`] sets the flag
/// and no operation resets it. Equality compares token identity, so
/// cancellation does not change it.
#[derive(Clone, Debug)]
pub struct CancellationToken(Arc<AtomicBool>);

impl Default for CancellationToken {
    fn default() -> Self {
        CancellationToken(Arc::new(AtomicBool::new(false)))
    }
}

impl CancellationToken {
    /// Make an unset cancellation token.
    pub fn new() -> Self {
        CancellationToken::default()
    }

    /// Request cancellation.
    pub fn cancel(&self) {
        self.0.store(true, Ordering::Relaxed);
    }

    /// Report whether cancellation was requested.
    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }

    /// Return the atomic flag shared with this token.
    ///
    /// Pass the clone to [`crate::verify::Limits::cancellation`] when one
    /// cancellation request must stop an independent verification. Only set
    /// the returned flag to `true`; the token never resets it.
    pub fn shared_flag(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.0)
    }

    /// Return the flag shared with this token inside the crate.
    pub(crate) fn flag(&self) -> Arc<AtomicBool> {
        self.shared_flag()
    }
}

impl PartialEq for CancellationToken {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl Eq for CancellationToken {}

/// The deadline and the memory limit of one call.
///
/// The default budget stops nothing. A budget holds a duration and not an
/// instant: the deadline starts when the call that receives it starts. It
/// can also hold a shared cancellation token. [`ComputeOptions`] holds one.
///
/// ```
/// use std::time::Duration;
/// use sylvester::{Budget, ComputeOptions};
///
/// let budget = Budget::new()
///     .timeout(Duration::from_secs(30))
///     .memory_limit(512 << 20);
/// let options = ComputeOptions::new().budget(budget);
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Budget {
    timeout: Option<Duration>,
    memory_limit: Option<usize>,
    cancellation: Option<CancellationToken>,
}

impl Budget {
    /// The budget that stops nothing: no deadline, no memory limit, and no
    /// cancellation token.
    pub fn new() -> Self {
        Budget::default()
    }

    /// Stop the call after this much time.
    ///
    /// The deadline starts when the call starts. A `timeout` too large to
    /// add to that instant leaves the call with no deadline, instead of
    /// stopping it before it starts.
    pub fn timeout(mut self, timeout: Duration) -> Self {
        self.timeout = Some(timeout);
        self
    }

    /// Stop the call once the data it tracks passes this many bytes.
    ///
    /// The count is an estimate of the live data of the call. It is not a
    /// cap on the process.
    pub fn memory_limit(mut self, bytes: usize) -> Self {
        self.memory_limit = Some(bytes);
        self
    }

    /// Stop the call when `token` is cancelled.
    ///
    /// Cloning the token lets another thread request cancellation while the
    /// call runs. A cancelled call reports [`ComputeError::Timeout`].
    pub fn cancellation(mut self, token: CancellationToken) -> Self {
        self.cancellation = Some(token);
        self
    }

    /// The instant the deadline stands at, taken now.
    fn deadline(&self) -> Option<Instant> {
        self.timeout
            .and_then(|timeout| Instant::now().checked_add(timeout))
    }
}

/// The backend, the budget, and the thread count of one computation.
///
/// The default is the F4 backend, no deadline, no memory limit, no
/// cancellation token, and the thread count of the rayon global pool. Each
/// setter takes the value by value and returns the options, so a call chain
/// builds them.
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
    budget: Budget,
    threads: Option<usize>,
}

impl ComputeOptions {
    /// The default options: the F4 backend, no deadline, no memory limit, no
    /// cancellation token, and the thread count of the rayon global pool.
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
        self.budget = self.budget.timeout(timeout);
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
        self.budget = self.budget.memory_limit(bytes);
        self
    }

    /// Take the deadline, memory limit, and cancellation token from one budget.
    ///
    /// This overwrites what [`ComputeOptions::timeout`] and
    /// [`ComputeOptions::memory_limit`] wrote, because those two setters
    /// write into the budget the options hold.
    pub fn budget(mut self, budget: Budget) -> Self {
        self.budget = budget;
        self
    }

    /// Compute on this many threads.
    ///
    /// The default is the rayon global pool, which reads
    /// `RAYON_NUM_THREADS`. A count of 0 means the same as the default.
    /// A count of 1 keeps the whole computation on the calling thread.
    /// A count above 1 runs it in a rayon pool of that size, which the
    /// call builds and drops. Building the pool costs more than a small
    /// computation itself takes.
    ///
    /// The F4 backend reduces the rows of a batch in parallel once the
    /// batch passes an internal work threshold. The classic backend, and
    /// every run that writes a certificate, computes on one thread
    /// whatever the count says. The thread count changes no byte of a
    /// basis or of a certificate.
    ///
    /// The count is best effort. A machine that refuses the pool leaves
    /// the run on the calling thread, and
    /// [`ComputeReport::threads_used`] then reports 1.
    ///
    /// Over `Q` the count is the number of prime runs the driver starts at
    /// once, and each of those runs computes on one thread. The driver
    /// starts at most 256 of them, and
    /// [`ComputeReport::modular_concurrency`] reports the largest number
    /// it started.
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = (threads > 0).then_some(threads);
        self
    }
}

/// When the multimodular rational engine stops.
///
/// A stop observes what its variant names and nothing more. The default is
/// [`RationalStop::Unchanged`] with `extra` equal to 2.
///
/// `extra` is a [`NonZeroUsize`] because a rule that observes zero further
/// primes observes nothing.
///
/// ```
/// use std::num::NonZeroUsize;
/// use sylvester::{RationalOptions, RationalStop};
///
/// let extra = NonZeroUsize::new(3).expect("3 is not zero");
/// let options = RationalOptions::new().stop(RationalStop::ContainsInput { extra });
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RationalStop {
    /// Stop when `extra` further primes of the prevalent category leave
    /// the lifted basis unchanged.
    ///
    /// This is a heuristic. It establishes nothing about the ideal.
    Unchanged {
        /// The number of further primes that must leave the lift
        /// unchanged.
        extra: NonZeroUsize,
    },
    /// Stop as [`RationalStop::Unchanged`] does, and then check the basis
    /// over `Q`.
    ///
    /// Every generator must reduce to zero modulo the basis, and so must
    /// every S-polynomial of the basis. The basis is then the reduced
    /// Gröbner basis of an ideal that contains the input ideal. Equality is
    /// not established: the basis `{1}` passes both checks. The checks are
    /// exact arithmetic over `Q` and can cost more than the modular
    /// computation.
    ContainsInput {
        /// The number of further primes that must leave the lift
        /// unchanged.
        extra: NonZeroUsize,
    },
}

impl Default for RationalStop {
    fn default() -> Self {
        RationalStop::Unchanged {
            extra: NonZeroUsize::new(2).expect("2 is not zero"),
        }
    }
}

/// The options of one computation over `Q`.
///
/// It holds a [`ComputeOptions`], which names the backend of each prime
/// run, the budget of the whole run, and the number of prime runs the
/// driver starts at once. It adds the stopping rule.
///
/// ```
/// use sylvester::{ComputeOptions, RationalOptions};
///
/// let options = RationalOptions::new().compute(ComputeOptions::new().threads(4));
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RationalOptions {
    compute: ComputeOptions,
    stop: RationalStop,
}

impl RationalOptions {
    /// The default options: the default [`ComputeOptions`] and
    /// [`RationalStop::Unchanged`] with `extra` equal to 2.
    pub fn new() -> Self {
        RationalOptions::default()
    }

    /// Take the backend, the budget, and the thread count from these
    /// options.
    ///
    /// The thread count is the number of prime runs the driver starts at
    /// once. Each prime run computes on one thread, so the thread count of
    /// the whole call is the one named here.
    pub fn compute(mut self, options: ComputeOptions) -> Self {
        self.compute = options;
        self
    }

    /// Pick the stopping rule.
    pub fn stop(mut self, stop: RationalStop) -> Self {
        self.stop = stop;
        self
    }
}

/// What one computation did, next to the basis it produced.
///
/// [`Ideal::groebner_basis_with_report`] returns it. The report is for a
/// benchmark harness or a test. Nothing in the crate reads it back.
///
/// [`Ideal::groebner_basis_with_report`]: crate::Ideal::groebner_basis_with_report
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ComputeReport {
    /// The backend that ran.
    pub backend: Backend,
    /// What the F4 run counted, or `None` under the classic backend,
    /// which counts nothing.
    pub counters: Option<F4Counters>,
    /// The wall time of the engine call.
    ///
    /// It covers the engine, and under the F4 backend the conversion of
    /// the basis back to [`Polynomial`]. It covers neither the
    /// certificate writer nor the verifier.
    pub elapsed: Duration,
    /// The number of threads the run had.
    ///
    /// It is 1 for the classic backend, for a run that writes a
    /// certificate, and for `ComputeOptions::threads(1)`. Otherwise it is
    /// the size of the pool the call built, or the size of the rayon
    /// global pool when the options name no count. A machine that refuses
    /// the pool leaves the run on the calling thread, and the value is
    /// then 1. The F4 kernel uses the threads only for a batch above its
    /// internal work threshold, so a run can have threads it never uses.
    pub threads_used: usize,
    /// What the multimodular run did, or `None` over a prime field, where
    /// no modular run happened.
    pub modular: Option<ModularLift>,
    /// The largest number of prime runs the driver started at once, or
    /// `None` over a prime field.
    ///
    /// It is what the driver started, not what [`ComputeOptions::threads`]
    /// asked for. The driver starts at most 256 runs at once, a machine
    /// that refuses a thread leaves that run on the calling thread, and a
    /// run that needs no prime, such as the one over the zero ideal,
    /// reports 0.
    pub modular_concurrency: Option<usize>,
}

/// The counts one F4 run reported.
///
/// A count is a report of what happened. It is never an input to a
/// decision, and a later release may change what a run counts. The
/// classic backend counts nothing, so it reports no value of this type.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct F4Counters {
    /// The batches the run reduced, the final interreduction included.
    pub batches: u64,
    /// The rows of every batch.
    pub matrix_rows: u64,
    /// The columns of every batch.
    pub matrix_columns: u64,
    /// The nonzero entries of every batch, before reduction.
    pub matrix_nonzeros: u64,
    /// The rows that reduced to zero.
    pub zero_rows: u64,
    /// The pivots the reduction installed.
    pub new_pivots: u64,
    /// The times the run started again at a wider exponent packing.
    pub lane_restarts: u32,
    /// The batches the run built again with half the pairs, because the
    /// batch it built first did not fit the memory limit.
    pub batch_retries: u64,
    /// The monomials the basis holds at the end of the run.
    pub basis_monomials: usize,
    /// The candidate critical pairs the run formed.
    pub pairs_generated: u64,
    /// The candidates dropped because the two leading monomials are
    /// coprime.
    pub pairs_discarded_product: u64,
    /// The queued pairs deleted by criterion B.
    pub pairs_discarded_b: u64,
    /// The candidates dropped by criterion M.
    pub pairs_discarded_m: u64,
    /// The candidates dropped by criterion F.
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

/// The largest value one exponent holds, which is the bound both backend
/// rules come from.
///
/// A monomial holds each exponent in a `u16`. The classic backend bounds
/// the total degree of a critical pair by it, and reports
/// [`ComputeError::DegreeLimit`]: a pair whose leading monomials have a
/// least common multiple at or below this degree can only produce
/// monomials whose own exponents and degree stay inside the width,
/// because every monomial the engines build is a term-multiple bounded by
/// an existing term or by an S-pair least common multiple. The F4 backend
/// packs each exponent on its own and bounds one exponent by it, and
/// reports [`ComputeError::ExponentLimit`].
pub(crate) const DEGREE_LIMIT: u32 = u16::MAX as u32;

fn check_input_degrees(generators: &[Polynomial]) -> Result<(), RunError> {
    if generators
        .iter()
        .any(|poly| poly.degree().is_some_and(|degree| degree > DEGREE_LIMIT))
    {
        return Err(RunError::Compute(ComputeError::DegreeLimit {
            limit: DEGREE_LIMIT,
        }));
    }
    Ok(())
}

/// Why a computation stops before it has a basis.
///
/// [`ComputeError::Timeout`] and [`ComputeError::MemoryLimitExceeded`]
/// report an exhausted budget; neither says anything about the ideal, and
/// the computation stopped rather than failed. [`ComputeError::DegreeLimit`]
/// and [`ComputeError::ExponentLimit`] report an input past what this
/// release supports, one per backend rule.
/// [`ComputeError::TableFull`] reports an F4 run that interned more
/// monomials than one table holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComputeError {
    /// The deadline passed.
    Timeout,
    /// The live engine data passed the memory limit.
    MemoryLimitExceeded,
    /// A critical pair reached a least common multiple of total degree
    /// above `limit`.
    ///
    /// The classic backend reports it. It bounds the degree of a pair,
    /// because the monomials it builds stay inside that bound.
    DegreeLimit {
        /// The largest total degree a critical pair may reach, which is
        /// what one exponent's width supports.
        limit: u32,
    },
    /// A monomial the run needs holds an exponent above `limit`.
    ///
    /// The F4 backend reports it. It packs each exponent on its own, so
    /// one exponent, not the total degree, is what it bounds.
    ExponentLimit {
        /// The largest value one exponent holds.
        limit: u32,
    },
    /// An F4 monomial table reached the largest number of monomials it
    /// holds, which is `u32::MAX - 1`.
    TableFull,
    /// The multimodular rational engine ran out of primes.
    ///
    /// The prime sequence walks down from 2^31 - 1 and ends below 3. A run
    /// that reaches the end has consumed every prime the sequence holds.
    PrimesExhausted,
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
                "a critical pair reached total degree above {limit}, which an exponent's width does not support"
            ),
            ComputeError::ExponentLimit { limit } => write!(
                f,
                "a monomial reached an exponent above {limit}, which is the largest one exponent holds"
            ),
            ComputeError::TableFull => {
                f.write_str("a monomial table reached the largest number of monomials it holds")
            }
            ComputeError::PrimesExhausted => {
                f.write_str("the rational computation ran out of primes")
            }
        }
    }
}

impl std::error::Error for ComputeError {}

/// The absolute limits of one call: the deadline, the memory limit, the
/// thread count, and the cancellation flag.
///
/// The deadline is an instant, taken once. Every part of one call holds to
/// that one instant: the engine, the certificate writer, and the verifier.
/// A part that took the timeout of [`ComputeOptions`] for itself would
/// start the whole duration again.
///
/// Only the crate builds this value. A caller gives a [`Budget`], which
/// holds a duration.
#[derive(Clone, Debug, Default)]
pub(crate) struct ComputeLimits {
    /// The instant the call stops at, or `None` for no deadline.
    pub(crate) deadline: Option<Instant>,
    /// The bytes the call may hold, or `None` for no limit.
    pub(crate) memory: Option<usize>,
    /// The thread count the caller asked for, or `None` for the rayon
    /// global pool.
    pub(crate) threads: Option<usize>,
    /// The flag that stops the call where its deadline is read, or `None`
    /// for a call nothing cancels.
    ///
    /// A public budget stores its flag here. A modular wave replaces this
    /// with a fresh private flag and stores the public flag in
    /// `external_cancel`.
    pub(crate) cancel: Option<Arc<AtomicBool>>,
    /// The caller's flag when `cancel` belongs to a modular wave.
    ///
    /// The two flags are read together. The private wave flag is never
    /// reset, so it cannot cancel a later wave.
    pub(crate) external_cancel: Option<Arc<AtomicBool>>,
}

/// The number of items one loop takes between two reads of the clock.
///
/// [`ComputeLimits::stop_every`] is what a loop calls per item.
pub(crate) const TICK: usize = 64;

impl ComputeLimits {
    /// The limits `options` names, with the deadline and cancellation flag
    /// taken now.
    pub(crate) fn of(options: &ComputeOptions) -> Self {
        ComputeLimits {
            threads: options.threads,
            ..ComputeLimits::of_budget(&options.budget)
        }
    }

    /// The limits `budget` names, with the deadline and cancellation flag
    /// taken now, and the thread count of the rayon global pool.
    ///
    /// A call that takes a [`Budget`] and no options, as
    /// [`GroebnerBasis::normal_form`] does, builds its limits here.
    ///
    /// [`GroebnerBasis::normal_form`]: crate::GroebnerBasis::normal_form
    pub(crate) fn of_budget(budget: &Budget) -> Self {
        ComputeLimits {
            deadline: budget.deadline(),
            memory: budget.memory_limit,
            threads: None,
            cancel: budget.cancellation.as_ref().map(CancellationToken::flag),
            external_cancel: None,
        }
    }

    /// Report why the call stops now, or `None` to carry on.
    ///
    /// This reads the clock once every [`TICK`] items. `index` counts the
    /// items of the loop that calls it, from 0, so the first item of
    /// every loop is a full check. Reading the clock costs more than one
    /// pass of the loops that call this per item.
    pub(crate) fn stop_every(&self, index: usize) -> Option<RunError> {
        if !index.is_multiple_of(TICK) {
            return None;
        }
        self.stop()
    }

    /// Report why the call stops now, or `None` to carry on.
    ///
    /// Cancellation flags are read before the clock, because reading the
    /// clock costs more.
    pub(crate) fn stop(&self) -> Option<RunError> {
        let cancelled = self
            .cancel
            .as_ref()
            .is_some_and(|flag| flag.load(Ordering::Relaxed))
            || self
                .external_cancel
                .as_ref()
                .is_some_and(|flag| flag.load(Ordering::Relaxed));
        if cancelled {
            return Some(RunError::Cancelled);
        }
        match self.deadline {
            Some(at) if Instant::now() >= at => Some(RunError::Compute(ComputeError::Timeout)),
            _ => None,
        }
    }
}

/// Why one backend run stops, including the cancellation the public errors
/// do not name.
///
/// [`ComputeError`] is public and exhaustive, so it gains no cancellation
/// variant: a caller who never asked for cancellation would have to match
/// one. The crate keeps the difference here. A caller that sets a token
/// consumes [`RunError::Cancelled`] through [`RunError::reported`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RunError {
    /// The run stopped for a reason the public error names.
    Compute(ComputeError),
    /// The cancellation flag of the run was set.
    Cancelled,
}

impl RunError {
    /// The public report of a stopped run.
    ///
    /// A cancelled run reports [`ComputeError::Timeout`]: it stopped
    /// before it had a basis, and that is what the public error says.
    pub(crate) fn reported(self) -> ComputeError {
        match self {
            RunError::Compute(error) => error,
            RunError::Cancelled => ComputeError::Timeout,
        }
    }
}

impl From<RunError> for CertifyError {
    fn from(error: RunError) -> Self {
        CertifyError::Engine(error.reported())
    }
}

impl From<crate::poly::ExponentOverflow> for RunError {
    fn from(_: crate::poly::ExponentOverflow) -> Self {
        RunError::Compute(ComputeError::DegreeLimit {
            limit: DEGREE_LIMIT,
        })
    }
}

/// Run the backend the options name, under limits this call takes.
pub(crate) fn groebner_basis(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &ComputeOptions,
) -> Result<Vec<Polynomial>, ComputeError> {
    let limits = ComputeLimits::of(options);
    groebner_basis_with_limits(ring, generators, options, &limits).map_err(RunError::reported)
}

/// Run the backend the options name, under limits the caller took.
///
/// The basis is the one [`groebner_basis`] returns for the same options.
/// A caller that runs the backend many times takes one [`ComputeLimits`]
/// and passes it to every run, so the deadline covers the whole of what it
/// does and not each run on its own.
pub(crate) fn groebner_basis_with_limits(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &ComputeOptions,
    limits: &ComputeLimits,
) -> Result<Vec<Polynomial>, RunError> {
    match options.backend {
        Backend::F4 => f4::groebner_basis(ring, generators, limits),
        Backend::Classic => classic::solve_checked(ring, generators, limits),
    }
}

/// Run the backend the options name and report what it did.
///
/// The basis is the one [`groebner_basis`] returns for the same options.
/// The report adds the counters of the run and the wall time of the
/// engine call.
pub(crate) fn groebner_basis_with_report(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &ComputeOptions,
) -> Result<(Vec<Polynomial>, ComputeReport), ComputeError> {
    let limits = ComputeLimits::of(options);
    let start = Instant::now();
    let (polynomials, counters, threads) = match options.backend {
        Backend::F4 => {
            let (basis, counters) =
                f4::solve(ring, generators, &limits).map_err(RunError::reported)?;
            (basis, Some(F4Counters::from(counters)), counters.threads)
        }
        Backend::Classic => (
            classic::solve_checked(ring, generators, &limits).map_err(RunError::reported)?,
            None,
            Some(1),
        ),
    };
    let report = ComputeReport {
        backend: options.backend,
        counters,
        elapsed: start.elapsed(),
        threads_used: threads.unwrap_or_else(rayon::current_num_threads),
        modular: None,
        modular_concurrency: None,
    };
    Ok((polynomials, report))
}

/// Compute the reduced Gröbner basis over `Q`.
pub(crate) fn rational_groebner_basis(
    ring: &PolynomialRing<Rationals>,
    generators: &[Polynomial<Rationals>],
    options: &RationalOptions,
) -> Result<(Vec<Polynomial<Rationals>>, ModularLift), ComputeError> {
    let outcome = modular::groebner_basis(ring, generators, options)?;
    Ok((outcome.basis, outcome.lift))
}

/// Compute the reduced Gröbner basis over `Q` and report what the run did.
///
/// The report carries no F4 counters. They describe one prime run, a
/// rational run has many, and summing them would name a run that never
/// happened.
pub(crate) fn rational_groebner_basis_with_report(
    ring: &PolynomialRing<Rationals>,
    generators: &[Polynomial<Rationals>],
    options: &RationalOptions,
) -> Result<(Vec<Polynomial<Rationals>>, ModularLift, ComputeReport), ComputeError> {
    let start = Instant::now();
    let outcome = modular::groebner_basis(ring, generators, options)?;
    let report = ComputeReport {
        backend: options.compute.backend,
        counters: None,
        elapsed: start.elapsed(),
        threads_used: 1,
        modular: Some(outcome.lift),
        modular_concurrency: Some(outcome.concurrency),
    };
    Ok((outcome.basis, outcome.lift, report))
}

/// Run the backend the options name, write its certificate, and verify
/// it.
///
/// The certificate format follows the backend: the classic backend writes
/// `sylv-gb-cert-v1` from the cofactors it tracks, and the F4 backend
/// writes `sylv-gb-cert-v2` from the trace it records.
///
/// `crate::cert` only assembles the certificate bytes; it never imports
/// [`crate::verify`]. This is the one place that reads the writer's
/// output back through the independent verifier, so a defect in either
/// module cannot hide behind a shared assumption.
///
/// The deadline and the memory limit of `options` are the budget of the
/// whole run. They are taken once and threaded through the engine, the
/// certificate writer, and the verifier caps, so the parts of one run
/// share a budget instead of each taking the whole of it. An exhausted
/// budget is [`CertifyError::Engine`] before the verifier starts and
/// [`CertifyError::VerifierExhausted`] inside it. Neither says anything
/// about the basis.
pub(crate) fn groebner_basis_certified(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &ComputeOptions,
) -> Result<CertifiedGroebnerBasis, CertifyError> {
    let limits = ComputeLimits::of(options);

    let (certificate, held) = match options.backend {
        Backend::Classic => classic_certificate(ring, generators, &limits)?,
        Backend::F4 => f4_certificate(ring, generators, &limits)?,
    };
    if let Some(stop) = limits.stop() {
        return Err(CertifyError::WriterExhausted(stop.reported()));
    }

    let verifier = verifier_limits(certificate.len(), held, ring.nvars(), &limits);
    let verified = verify::verify_with_limits(&certificate, &verifier)?;

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
        GroebnerBasis::new(ring.clone(), polynomials, ()),
        certificate,
    ))
}

/// Run the classic backend and write its `sylv-gb-cert-v1` certificate.
///
/// The backend tracks the origin cofactors through reduction, and the
/// emitter writes them. The second value is what the emitter's budget
/// holds, which sizes the verifier's caps.
fn classic_certificate(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    limits: &ComputeLimits,
) -> Result<(Vec<u8>, usize), CertifyError> {
    let (basis, origins) = classic::solve_tracked(ring, generators, limits)?;
    let mut budget = WriterBudget::new(limits, ring.nvars());
    let certificate = cert::assemble(ring, generators, basis, origins, &mut budget)?;
    Ok((certificate, budget.held()))
}

/// Run the F4 backend and write its `sylv-gb-cert-v2` certificate.
///
/// The recorder holds the trace of the run, so the run stays on one
/// thread whatever the thread count says.
///
/// `limits` is the budget of the whole certified run. The engine, the
/// recorder, and the writer hold to that one deadline and that one memory
/// limit: the engine's estimate counts the recorder's nodes, and the
/// writer carries on with the recorder's budget.
fn f4_certificate(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    limits: &ComputeLimits,
) -> Result<(Vec<u8>, usize), CertifyError> {
    let mut recorder = cert::v2::Recorder::new(ring, generators, limits);
    let (basis, _) = f4::solve_recorded(ring, generators, limits, &mut recorder)?;
    cert::v2::assemble(ring, generators, &basis, recorder)
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
/// length bounds every count the caps name. The work cap stays open on
/// both paths: it guards against hostile bytes, and here the deadline
/// bounds the time.
pub(crate) fn verifier_limits(
    certificate_len: usize,
    held: usize,
    nvars: usize,
    limits: &ComputeLimits,
) -> VerifyLimits {
    let deadline = limits.deadline;
    let cancellation = limits
        .external_cancel
        .clone()
        .or_else(|| limits.cancel.clone());
    let defaults = VerifyLimits::default();
    let Some(limit) = limits.memory else {
        let len = certificate_len;
        return VerifyLimits {
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
            // The work cap is a guard against hostile bytes. These bytes
            // come from this process, and the deadline bounds the time.
            max_work_units: u64::MAX,
            deadline,
            cancellation,
        };
    };

    let left = limit.saturating_sub(held) / 3;
    let term_bytes = size_of::<verify::Term>() + nvars * size_of::<verify::Exp>();
    let terms = (left / term_bytes).max(1);
    let polys = (left / size_of::<verify::Poly>()).max(1);
    let slots = (left / 64).max(1);
    VerifyLimits {
        max_bytes: limit.saturating_sub(held) + certificate_len,
        max_polys: polys,
        max_terms_per_poly: terms,
        max_total_terms: terms,
        max_entries: polys,
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
        cancellation,
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

    /// Limits of `bytes` memory and no deadline.
    fn bounded_limits(bytes: usize) -> ComputeLimits {
        ComputeLimits {
            memory: Some(bytes),
            ..ComputeLimits::default()
        }
    }

    #[test]
    fn a_timeout_too_large_to_add_to_now_leaves_no_deadline() {
        let options = ComputeOptions::new().timeout(Duration::MAX);
        assert_eq!(ComputeLimits::of(&options).deadline, None);
    }

    #[test]
    fn a_thread_count_of_zero_means_the_default_pool() {
        assert_eq!(ComputeOptions::new().threads(0), ComputeOptions::new());
        assert_ne!(ComputeOptions::new().threads(1), ComputeOptions::new());
    }

    #[test]
    fn an_ordinary_timeout_still_sets_a_deadline() {
        let options = ComputeOptions::new().timeout(Duration::from_secs(60));
        assert!(ComputeLimits::of(&options).deadline.is_some());
    }

    #[test]
    fn the_verifier_limits_come_from_what_is_left_of_the_budget() {
        let bytes = vec![0u8; 1000];
        assert_unlimited_verifier_limits(&bytes);
        assert_long_certificate_raises_caps();
        let bounded = assert_bounded_verifier_limits(&bytes);
        assert_spent_verifier_limits(&bytes);
        assert_held_bytes_reduce_term_caps(&bytes, &bounded);
    }

    fn assert_unlimited_verifier_limits(bytes: &[u8]) {
        let unlimited = verifier_limits(bytes.len(), bytes.len(), 2, &ComputeLimits::default());
        assert_eq!(unlimited.max_bytes, VerifyLimits::default().max_bytes);
        assert_eq!(
            unlimited.max_total_terms,
            VerifyLimits::default().max_total_terms
        );
    }

    fn assert_long_certificate_raises_caps() {
        let long = vec![0u8; VerifyLimits::default().max_bytes + 1];
        let raised = verifier_limits(long.len(), long.len(), 2, &ComputeLimits::default());
        assert_eq!(raised.max_bytes, long.len());
        assert_eq!(raised.max_total_terms, long.len());
        assert_eq!(raised.max_polys, long.len());
    }

    fn assert_bounded_verifier_limits(bytes: &[u8]) -> VerifyLimits {
        let bounded = verifier_limits(bytes.len(), bytes.len(), 2, &bounded_limits(1 << 20));
        assert_eq!(bounded.max_bytes, 1 << 20);
        assert!(bounded.max_total_terms < VerifyLimits::default().max_total_terms);
        assert!(bounded.max_total_terms > 0);
        bounded
    }

    fn assert_spent_verifier_limits(bytes: &[u8]) {
        let spent = verifier_limits(bytes.len(), bytes.len(), 2, &bounded_limits(bytes.len()));
        assert_eq!(spent.max_bytes, bytes.len());
        assert_eq!(spent.max_polys, 1);
        assert_eq!(spent.max_total_terms, 1);
        assert_eq!(
            spent.max_intermediate_bytes,
            size_of::<verify::Term>() + 2 * size_of::<verify::Exp>()
        );
    }

    fn assert_held_bytes_reduce_term_caps(bytes: &[u8], bounded: &VerifyLimits) {
        let held = bytes.len() + 500;
        let with_other_data = verifier_limits(bytes.len(), held, 2, &bounded_limits(1 << 20));
        assert!(with_other_data.max_total_terms < bounded.max_total_terms);
    }

    /// A ring of `nvars` variables named `x0 .. x(nvars - 1)`.
    fn wide_ring(nvars: usize) -> PolynomialRing {
        let names: Vec<String> = (0..nvars).map(|index| format!("x{index}")).collect();
        PolynomialRing::prime_field(1073741827, names).expect("the modulus is prime")
    }

    /// cyclic-n in n variables, as `benchmarks/gb-comparison/gen.py`
    /// defines it.
    fn cyclic(ring: &PolynomialRing, n: usize) -> Vec<Polynomial> {
        let mut system = Vec::new();
        for d in 1..n {
            let mut terms = Vec::new();
            for i in 0..n {
                let mut exps = vec![0u16; n];
                for j in 0..d {
                    exps[(i + j) % n] += 1;
                }
                terms.push((1i64, exps));
            }
            system.push(ring.polynomial(terms).expect("the exponents fit"));
        }
        let last = vec![(1i64, vec![1u16; n]), (-1, vec![0u16; n])];
        system.push(ring.polynomial(last).expect("the exponents fit"));
        system
    }

    /// katsura-n in n + 1 variables, the POSSO definition.
    fn katsura(ring: &PolynomialRing, n: usize) -> Vec<Polynomial> {
        let nv = n + 1;
        let mut system = Vec::new();
        let mut first = Vec::new();
        for i in 0..=n {
            let mut exps = vec![0u16; nv];
            exps[i] = 1;
            first.push((if i == 0 { 1i64 } else { 2 }, exps));
        }
        first.push((-1, vec![0u16; nv]));
        system.push(ring.polynomial(first).expect("the exponents fit"));

        let bound = n as i64;
        for m in 0..n as i64 {
            let mut terms: Vec<(i64, Vec<u16>)> = Vec::new();
            for i in -bound..=bound {
                let j = m - i;
                if j.abs() > bound {
                    continue;
                }
                let mut exps = vec![0u16; nv];
                exps[i.unsigned_abs() as usize] += 1;
                exps[j.unsigned_abs() as usize] += 1;
                terms.push((1, exps));
            }
            let mut exps = vec![0u16; nv];
            exps[m as usize] = 1;
            terms.push((-1, exps));
            system.push(ring.polynomial(terms).expect("the exponents fit"));
        }
        system
    }

    fn cancelled_run(ring: &PolynomialRing, system: &[Polynomial], options: &ComputeOptions) {
        let flag = Arc::new(AtomicBool::new(true));
        let limits = ComputeLimits {
            cancel: Some(flag),
            ..ComputeLimits::of(options)
        };
        assert_eq!(
            groebner_basis_with_limits(ring, system, options, &limits),
            Err(RunError::Cancelled)
        );
    }

    #[test]
    fn periodic_checks_read_cancellation_at_the_next_tick() {
        let flag = Arc::new(AtomicBool::new(false));
        let limits = ComputeLimits {
            cancel: Some(Arc::clone(&flag)),
            ..ComputeLimits::default()
        };
        assert_eq!(limits.stop_every(0), None);
        flag.store(true, Ordering::Relaxed);
        assert_eq!(limits.stop_every(1), None);
        assert_eq!(limits.stop_every(TICK), Some(RunError::Cancelled));
    }

    #[test]
    fn a_cancelled_f4_run_stops_at_the_flag() {
        let ring = wide_ring(8);
        let system = katsura(&ring, 7);
        cancelled_run(&ring, &system, &ComputeOptions::new().threads(1));
    }

    #[test]
    fn a_cancelled_parallel_f4_run_stops_at_the_flag() {
        let ring = wide_ring(8);
        let system = katsura(&ring, 7);
        cancelled_run(&ring, &system, &ComputeOptions::new().threads(8));
    }

    #[test]
    fn a_cancelled_classic_run_stops_at_the_flag() {
        let ring = wide_ring(6);
        let system = katsura(&ring, 5);
        cancelled_run(
            &ring,
            &system,
            &ComputeOptions::new().backend(Backend::Classic),
        );
    }

    #[test]
    fn the_limits_of_a_call_give_the_basis_the_options_give() {
        let ring = wide_ring(5);
        let system = cyclic(&ring, 5);
        for backend in [Backend::F4, Backend::Classic] {
            let options = ComputeOptions::new()
                .backend(backend)
                .threads(1)
                .timeout(Duration::from_secs(60))
                .memory_limit(1 << 30);
            let limits = ComputeLimits::of(&options);
            assert_eq!(
                groebner_basis_with_limits(&ring, &system, &options, &limits)
                    .map_err(RunError::reported),
                groebner_basis(&ring, &system, &options),
                "{backend:?}"
            );
        }
    }

    /// The certified path takes one instant and the engine, the writer,
    /// and the verifier hold to it. A writer that took the timeout again
    /// would run for the whole of it after the engine had already spent
    /// part of it.
    #[test]
    fn one_deadline_covers_the_engine_the_writer_and_the_verifier() {
        let ring = wide_ring(6);
        let system = cyclic(&ring, 6);
        let options = ComputeOptions::new().backend(Backend::F4).threads(1);

        let start = Instant::now();
        groebner_basis(&ring, &system, &options).expect("no budget stops it");
        let engine = start.elapsed();

        // Halves of the engine time, so the first deadline the engine
        // fits inside leaves the writer, which costs several times more,
        // without the time to finish.
        for halves in [3, 4, 6] {
            let timeout = engine * halves / 2;
            let start = Instant::now();
            let stopped =
                groebner_basis_certified(&ring, &system, &options.clone().timeout(timeout));
            let elapsed = start.elapsed();
            match stopped {
                // The engine spent the deadline itself on this run. A
                // longer one leaves it to the writer.
                Err(CertifyError::Engine(ComputeError::Timeout)) => continue,
                Err(CertifyError::WriterExhausted(ComputeError::Timeout))
                | Err(CertifyError::VerifierExhausted(_)) => {
                    assert!(
                        elapsed < timeout + engine / 2,
                        "the whole call stops at one deadline: {elapsed:?} against {timeout:?}"
                    );
                    return;
                }
                other => panic!("the deadline must stop the writer, not {other:?}"),
            }
        }
        panic!("no deadline left the engine time and the writer none");
    }
}
