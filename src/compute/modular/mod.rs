//! The multimodular rational engine.
//!
//! One run clears the denominators of the generators, computes a basis
//! modulo many primes with the backend the options name, combines the runs
//! by the Chinese remainder theorem, and lifts each residue to a rational
//! number. The result is a heuristic: no isolated verifier checks the
//! vote, the combination, or the lift, and `RationalStop` says what a stop
//! observes.

mod check;
mod crt;
mod primes;
mod reconstruct;

use std::collections::BTreeMap;
use std::mem::size_of;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread;

use num_traits::Zero;

use crate::compute::{
    self, Backend, ComputeError, ComputeLimits, ComputeOptions, RationalOptions, RationalStop,
    RunError,
};
use crate::poly::{Monomial, Polynomial, Term, heap_exps_bytes};
use crate::ring::field::residue;
use crate::ring::rational::{ClearedPolynomial, limb_bytes};
use crate::ring::{Established, ModularLift, PolynomialRing, PrimeField, Rationals};

use crt::Accumulator;
use primes::{Prime, Sequence};

/// The most prime runs the driver starts at once.
///
/// A larger [`ComputeOptions::threads`] is clamped to this. Each
/// concurrent run takes an operating system thread and a share of the
/// memory limit, and the share shrinks as the count grows.
/// [`ComputeReport::modular_concurrency`] reports the value the driver
/// used, not the value the caller asked for.
///
/// [`ComputeOptions::threads`]: crate::ComputeOptions::threads
/// [`ComputeReport::modular_concurrency`]: crate::ComputeReport::modular_concurrency
const MAX_CONCURRENCY: usize = 256;

/// The basis one rational run produced, with the record of the run.
pub(crate) struct Outcome {
    /// The lifted basis, monic and sorted strictly descending by leading
    /// monomial.
    pub(crate) basis: Vec<Polynomial<Rationals>>,
    /// The counters of the run and what it establishes.
    pub(crate) lift: ModularLift,
    /// The largest number of prime runs the driver started at once, which
    /// is 0 when it started none.
    pub(crate) concurrency: usize,
}

/// Compute the reduced Gröbner basis over `Q`.
pub(crate) fn groebner_basis(
    ring: &PolynomialRing<Rationals>,
    generators: &[Polynomial<Rationals>],
    options: &RationalOptions,
) -> Result<Outcome, ComputeError> {
    let input = Input::new(ring, generators, options).map_err(RunError::reported)?;
    let mut state = State::new();
    state.run(&input).map_err(RunError::reported)
}

/// What one rational run reads and never changes.
struct Input<'a> {
    ring: &'a PolynomialRing<Rationals>,
    /// The nonzero generators, in the order the caller gave them.
    generators: Vec<Polynomial<Rationals>>,
    /// The generators with their denominators cleared and their content
    /// divided out.
    cleared: Vec<ClearedPolynomial>,
    /// The bytes `cleared` holds, which the ledger charges for the whole
    /// run.
    cleared_bytes: usize,
    /// The bytes `generators` holds. The driver keeps the rational
    /// generators for T1, so the ledger charges them next to `cleared`.
    generators_bytes: usize,
    backend: Backend,
    stop: RationalStop,
    /// The absolute deadline and the memory limit of the whole run.
    limits: ComputeLimits,
    /// The number of prime runs the driver starts at once.
    concurrency: usize,
}

impl<'a> Input<'a> {
    /// Read the generators and clear their denominators.
    ///
    /// The deadline covers this too, so it is read once per generator.
    fn new(
        ring: &'a PolynomialRing<Rationals>,
        generators: &[Polynomial<Rationals>],
        options: &RationalOptions,
    ) -> Result<Self, RunError> {
        let limits = ComputeLimits {
            threads: Some(1),
            ..ComputeLimits::of(&options.compute)
        };
        let generators: Vec<Polynomial<Rationals>> = generators
            .iter()
            .filter(|f| !f.is_zero())
            .cloned()
            .collect();
        let mut cleared = Vec::with_capacity(generators.len());
        let mut cleared_bytes = 0usize;
        let mut generators_bytes = 0usize;
        for generator in &generators {
            if let Some(stop) = limits.stop() {
                return Err(stop);
            }
            let polynomial = ClearedPolynomial::of(generator, &limits)?;
            cleared_bytes = cleared_bytes.saturating_add(polynomial.heap_bytes(ring.nvars()));
            generators_bytes = generators_bytes.saturating_add(generator.heap_bytes());
            cleared.push(polynomial);
            if let Some(limit) = limits.memory
                && cleared_bytes.saturating_add(generators_bytes) > limit
            {
                return Err(RunError::Compute(ComputeError::MemoryLimitExceeded));
            }
        }
        let concurrency = options
            .compute
            .threads
            .unwrap_or_else(rayon::current_num_threads)
            .clamp(1, MAX_CONCURRENCY);
        Ok(Input {
            ring,
            generators,
            cleared,
            cleared_bytes,
            generators_bytes,
            backend: options.compute.backend,
            stop: options.stop,
            limits,
            concurrency,
        })
    }
}

/// One completed prime run the driver consumed.
struct Run {
    prime: u64,
    /// The primes the driver skipped before this one.
    skipped_before: usize,
    /// The leading monomials of the basis, which name the category.
    category: Vec<Monomial>,
    basis: Vec<Polynomial<PrimeField>>,
    /// The bytes the basis holds.
    bytes: usize,
}

/// What the driver decided after it consumed one run.
enum Progress {
    /// Consume the next prime.
    Continue,
    /// The run reported its memory share exhausted. Run the same prime
    /// alone with the whole residual.
    RetryAlone,
    /// The stopping rule is met and the basis is the result.
    Stop(Vec<Polynomial<Rationals>>, Established),
    /// The stopping rule is met and the exact tests decide.
    Test(Vec<Polynomial<Rationals>>),
}

/// The state machine of section 3.7 of `docs/rational-design.md`.
struct State {
    sequence: Sequence,
    /// The number of primes of the sequence the driver consumed.
    consumed: usize,
    retained: Vec<Run>,
    /// One count per category, keyed by its leading monomials.
    counts: BTreeMap<Vec<Monomial>, usize>,
    prevalent: Option<Vec<Monomial>>,
    accumulator: Option<Accumulator>,
    candidate: Option<Vec<Polynomial<Rationals>>>,
    /// The runs that left the candidate unchanged.
    confirmations: usize,
    /// The held candidate failed the exact tests.
    rejected: bool,
    /// The next wave runs one prime with the whole residual.
    retry_alone: bool,
    /// The largest number of prime runs one wave started at once. It
    /// counts what the driver started, not what the caller asked for: a
    /// thread the machine refuses leaves the run on the calling thread.
    concurrency: usize,
}

impl State {
    fn new() -> Self {
        State {
            sequence: Sequence::new(),
            consumed: 0,
            retained: Vec::new(),
            counts: BTreeMap::new(),
            prevalent: None,
            accumulator: None,
            candidate: None,
            confirmations: 0,
            rejected: false,
            retry_alone: false,
            concurrency: 0,
        }
    }

    fn run(&mut self, input: &Input) -> Result<Outcome, RunError> {
        if input.cleared.is_empty() {
            return Ok(Outcome {
                basis: Vec::new(),
                lift: self.lift(Established::Unchanged),
                concurrency: self.concurrency,
            });
        }
        loop {
            match self.wave(input)? {
                Progress::Continue => {}
                Progress::RetryAlone => self.retry_alone = true,
                Progress::Stop(basis, established) => {
                    return Ok(Outcome {
                        basis,
                        lift: self.lift(established),
                        concurrency: self.concurrency,
                    });
                }
                Progress::Test(basis) => {
                    // Every prime run is joined by now, so the tests hold
                    // the whole residual, less the copy of the candidate
                    // they run on, which is live next to the held one. The
                    // flag is fresh, because a flag is never reset.
                    let copy = basis_bytes(&basis);
                    let limits = ComputeLimits {
                        memory: self.residual(input).map(|bytes| bytes.saturating_sub(copy)),
                        cancel: Some(Arc::new(AtomicBool::new(false))),
                        ..input.limits.clone()
                    };
                    if check::contains_input(&input.generators, &basis, &limits)? {
                        return Ok(Outcome {
                            basis,
                            lift: self.lift(Established::ContainsInput),
                            concurrency: self.concurrency,
                        });
                    }
                    self.rejected = true;
                }
            }
        }
    }

    /// Run one wave of primes and consume the results in sequence order.
    ///
    /// The wave holds `concurrency` primes, or one after a run reported
    /// its memory share exhausted. A result the driver does not consume is
    /// discarded whole, and its prime stays unconsumed, so the primes the
    /// driver consumes are the same at every thread count.
    fn wave(&mut self, input: &Input) -> Result<Progress, RunError> {
        if let Some(stop) = input.limits.stop() {
            return Err(stop);
        }
        let residual = self.residual(input);
        if residual == Some(0) {
            return Err(RunError::Compute(ComputeError::MemoryLimitExceeded));
        }

        let primes = self.next_primes(input)?;
        if primes.is_empty() {
            return Err(RunError::Compute(ComputeError::PrimesExhausted));
        }

        let flag = Arc::new(AtomicBool::new(false));
        let run_limits = ComputeLimits {
            memory: residual.map(|bytes| (bytes / primes.len()).max(1)),
            cancel: Some(flag.clone()),
            ..input.limits.clone()
        };
        let (outcome, started) = self.run_primes(input, &primes, &run_limits, &flag);
        self.concurrency = self.concurrency.max(started);
        outcome
    }

    fn run_primes(
        &mut self,
        input: &Input,
        primes: &[Prime],
        run_limits: &ComputeLimits,
        flag: &Arc<AtomicBool>,
    ) -> (Result<Progress, RunError>, usize) {
        let mut outcome: Result<Progress, RunError> = Ok(Progress::Continue);
        let mut spawned = 0usize;
        let mut on_this_thread = false;
        thread::scope(|scope| {
            let mut handles = Vec::with_capacity(primes.len());
            for prime in primes {
                let value = prime.value;
                let limits = run_limits.clone();
                let started = thread::Builder::new()
                    .spawn_scoped(scope, move || run_one(input, value, &limits));
                spawned += usize::from(started.is_ok());
                handles.push(started.ok());
            }
            for (offset, handle) in handles.into_iter().enumerate() {
                if !matches!(outcome, Ok(Progress::Continue)) {
                    flag.store(true, Ordering::Relaxed);
                    continue;
                }
                let result = join_run(
                    handle,
                    input,
                    primes[offset].value,
                    run_limits,
                    &mut on_this_thread,
                );
                outcome = self.advance_prime(input, primes[offset], result);
                if matches!(outcome, Ok(Progress::Continue)) {
                    self.retry_alone = false;
                }
            }
        });
        let started = spawned.saturating_add(usize::from(on_this_thread));
        (outcome, started)
    }

    fn next_primes(&mut self, input: &Input) -> Result<Vec<Prime>, RunError> {
        let width = if self.retry_alone {
            1
        } else {
            input.concurrency
        };
        let mut primes = Vec::with_capacity(width);
        for offset in 0..width {
            let next = self.sequence.at(
                self.consumed + offset,
                |value| usable(value, &input.cleared),
                &input.limits,
            )?;
            let Some(prime) = next else {
                break;
            };
            primes.push(prime);
        }
        Ok(primes)
    }

    fn advance_prime(
        &mut self,
        input: &Input,
        prime: Prime,
        result: Result<Vec<Polynomial<PrimeField>>, RunError>,
    ) -> Result<Progress, RunError> {
        match result {
            Ok(basis) => self.consume(input, prime, basis),
            Err(RunError::Compute(ComputeError::MemoryLimitExceeded)) if !self.retry_alone => {
                Ok(Progress::RetryAlone)
            }
            Err(error) => Err(error),
        }
    }

    /// Fold one completed run into the state and report what follows.
    fn consume(
        &mut self,
        input: &Input,
        prime: Prime,
        basis: Vec<Polynomial<PrimeField>>,
    ) -> Result<Progress, RunError> {
        let category = self.record_run(prime, basis);
        if !self.fold_run(input, prime.value, &category)? {
            return Ok(Progress::Continue);
        }
        let accumulator = self
            .accumulator
            .as_ref()
            .expect("a folded run holds an accumulator");
        if accumulator.folded() < 2 {
            return Ok(Progress::Continue);
        }

        let Some(lifted) = self.lift_candidate(input)? else {
            self.reset_candidate();
            return Ok(Progress::Continue);
        };
        self.update_candidate(lifted);
        if self.confirmations < confirmations(input.stop) {
            return Ok(Progress::Continue);
        }
        self.finish_candidate(input)
    }

    fn record_run(&mut self, prime: Prime, basis: Vec<Polynomial<PrimeField>>) -> Vec<Monomial> {
        let category = basis
            .iter()
            .map(|element| {
                element
                    .lm()
                    .cloned()
                    .expect("a basis element is not the zero polynomial")
            })
            .collect::<Vec<_>>();
        let bytes = basis.iter().fold(0usize, |bytes, element| {
            bytes.saturating_add(element.heap_bytes())
        });
        self.retained.push(Run {
            prime: prime.value,
            skipped_before: prime.skipped_before,
            category: category.clone(),
            basis,
            bytes,
        });
        self.consumed += 1;
        *self.counts.entry(category.clone()).or_insert(0) += 1;
        category
    }

    fn fold_run(
        &mut self,
        input: &Input,
        prime: u64,
        category: &[Monomial],
    ) -> Result<bool, RunError> {
        let prevalent = self.prevalent();
        if self.prevalent.as_ref() != Some(&prevalent) {
            self.prevalent = Some(prevalent);
            self.rebuild(input)?;
            self.reset_candidate();
            return Ok(true);
        }
        if category != prevalent {
            return Ok(false);
        }
        let index = self.retained.len() - 1;
        let accumulator = self
            .accumulator
            .as_mut()
            .expect("a prevalent category holds an accumulator");
        accumulator.fold(&self.retained[index].basis, prime, &input.limits)?;
        Ok(true)
    }

    fn reset_candidate(&mut self) {
        self.candidate = None;
        self.rejected = false;
        self.confirmations = 0;
    }

    fn update_candidate(&mut self, lifted: Vec<Polynomial<Rationals>>) {
        if self.candidate.as_ref() == Some(&lifted) {
            self.confirmations += 1;
            return;
        }
        self.candidate = Some(lifted);
        self.rejected = false;
        self.confirmations = 0;
    }

    fn finish_candidate(&mut self, input: &Input) -> Result<Progress, RunError> {
        let held = self
            .candidate
            .as_ref()
            .expect("the counter counts a held candidate");
        // The clone is a second copy of the candidate, live next to the
        // held one until the run ends.
        let copy = basis_bytes(held);
        self.charge(input, copy)?;
        let candidate = held.clone();
        match input.stop {
            RationalStop::Unchanged { .. } => Ok(Progress::Stop(candidate, Established::Unchanged)),
            // The tests ran on this candidate already and rejected it.
            // Running them again would cost the same and observe the same.
            RationalStop::ContainsInput { .. } if self.rejected => Ok(Progress::Continue),
            RationalStop::ContainsInput { .. } => Ok(Progress::Test(candidate)),
        }
    }

    /// The category with the most runs.
    ///
    /// A tie goes to the lexicographically smallest list of leading
    /// monomials, which makes the choice deterministic and nothing more.
    fn prevalent(&self) -> Vec<Monomial> {
        let mut best: Option<(&Vec<Monomial>, usize)> = None;
        for (category, &count) in &self.counts {
            match best {
                Some((_, held)) if held >= count => {}
                _ => best = Some((category, count)),
            }
        }
        best.expect("a consumed run holds a category").0.clone()
    }

    /// Build the accumulator again from the retained runs of the prevalent
    /// category, in prime order.
    fn rebuild(&mut self, input: &Input) -> Result<(), RunError> {
        let prevalent = self
            .prevalent
            .clone()
            .expect("a rebuild follows a prevalent category");
        let mut accumulator = Accumulator::new(&prevalent);
        let nvars = input.ring.nvars();
        for run in &self.retained {
            if run.category == prevalent {
                accumulator.fold(&run.basis, run.prime, &input.limits)?;
                // The held accumulator is not dropped until the rebuild
                // ends, so the scratch one is charged next to it.
                self.charge(input, accumulator.heap_bytes(nvars))?;
            }
        }
        self.accumulator = Some(accumulator);
        Ok(())
    }

    /// Lift every residue of the accumulator, or report that one of them
    /// does not lift.
    fn lift_candidate(
        &self,
        input: &Input,
    ) -> Result<Option<Vec<Polynomial<Rationals>>>, RunError> {
        let accumulator = self
            .accumulator
            .as_ref()
            .expect("a lift follows a folded run");
        let modulus = accumulator.modulus();
        let nvars = input.ring.nvars();
        let mut basis = Vec::with_capacity(accumulator.len());
        // The held candidate is still live, so the basis this builds is
        // charged as it grows.
        let mut scratch = 0usize;
        for index in 0..accumulator.len() {
            let mut terms: Vec<Term<Rationals>> = Vec::new();
            for (mono, residue) in accumulator.residues(index) {
                let Some(coeff) = reconstruct::reconstruct(residue, modulus, &input.limits)? else {
                    return Ok(None);
                };
                if coeff.is_zero() {
                    continue;
                }
                scratch = scratch.saturating_add(term_bytes(&coeff, nvars));
                self.charge(input, scratch)?;
                terms.push(Term {
                    coeff,
                    mono: mono.clone(),
                });
            }
            basis.push(Polynomial::from_sorted_terms(input.ring.clone(), terms));
        }
        Ok(Some(basis))
    }

    /// The bytes the run may still hold, or `None` under no memory limit.
    fn residual(&self, input: &Input) -> Option<usize> {
        input
            .limits
            .memory
            .map(|limit| limit.saturating_sub(self.held(input)))
    }

    /// The bytes the driver holds.
    ///
    /// The count covers everything the driver keeps for the whole run: the
    /// rational generators, the cleared integer generators, every retained
    /// modular basis with its category, the category table, the primes the
    /// sequence found, the accumulator, and the held candidate.
    fn held(&self, input: &Input) -> usize {
        let nvars = input.ring.nvars();
        let start = input.cleared_bytes.saturating_add(input.generators_bytes);
        let retained = self.retained.iter().fold(start, |bytes, run| {
            bytes
                .saturating_add(run.bytes)
                .saturating_add(category_bytes(&run.category, nvars))
        });
        let counts = self.counts.keys().fold(0usize, |bytes, category| {
            bytes
                .saturating_add(size_of::<(Vec<Monomial>, usize)>())
                .saturating_add(category_bytes(category, nvars))
        });
        let accumulator = self
            .accumulator
            .as_ref()
            .map_or(0, |accumulator| accumulator.heap_bytes(nvars));
        let candidate = self.candidate.as_deref().map_or(0, basis_bytes);
        retained
            .saturating_add(counts)
            .saturating_add(self.sequence.heap_bytes())
            .saturating_add(accumulator)
            .saturating_add(candidate)
    }

    /// Report why the run stops now, or `Ok(())` to carry on.
    ///
    /// `extra` is what the caller holds on top of the ledger: a scratch
    /// value the driver is building next to the one it already keeps.
    fn charge(&self, input: &Input, extra: usize) -> Result<(), RunError> {
        if let Some(limit) = input.limits.memory
            && self.held(input).saturating_add(extra) > limit
        {
            return Err(RunError::Compute(ComputeError::MemoryLimitExceeded));
        }
        Ok(())
    }

    /// The record of the run so far.
    fn lift(&self, established: Established) -> ModularLift {
        let folded = self.accumulator.as_ref().map_or(0, Accumulator::folded);
        ModularLift {
            primes_consumed: self.consumed,
            primes_skipped: self.retained.last().map_or(0, |run| run.skipped_before),
            primes_folded: folded,
            primes_discarded: self.consumed - folded,
            confirming_primes: self.confirmations,
            modulus_bits: self
                .accumulator
                .as_ref()
                .map_or(0, |accumulator| accumulator.modulus().bits()),
            established,
        }
    }
}

fn confirmations(stop: RationalStop) -> usize {
    match stop {
        RationalStop::Unchanged { extra } | RationalStop::ContainsInput { extra } => extra.get(),
    }
}

fn join_run(
    handle: Option<thread::ScopedJoinHandle<'_, Result<Vec<Polynomial<PrimeField>>, RunError>>>,
    input: &Input,
    prime: u64,
    limits: &ComputeLimits,
    on_this_thread: &mut bool,
) -> Result<Vec<Polynomial<PrimeField>>, RunError> {
    match handle {
        Some(handle) => handle
            .join()
            .unwrap_or_else(|payload| std::panic::resume_unwind(payload)),
        None => {
            *on_this_thread = true;
            run_one(input, prime, limits)
        }
    }
}

/// The bytes one category holds.
fn category_bytes(category: &[Monomial], nvars: usize) -> usize {
    category
        .len()
        .saturating_mul(size_of::<Monomial>().saturating_add(heap_exps_bytes(nvars)))
}

/// The bytes one lifted basis holds.
fn basis_bytes(basis: &[Polynomial<Rationals>]) -> usize {
    basis.iter().fold(0usize, |bytes, element| {
        bytes.saturating_add(element.heap_bytes())
    })
}

/// The bytes one lifted term holds.
fn term_bytes(coeff: &num_rational::BigRational, nvars: usize) -> usize {
    size_of::<Term<Rationals>>()
        .saturating_add(heap_exps_bytes(nvars))
        .saturating_add(limb_bytes(coeff.numer()))
        .saturating_add(limb_bytes(coeff.denom()))
}

/// Compute the basis of the images of the generators modulo `prime`.
fn run_one(
    input: &Input,
    prime: u64,
    limits: &ComputeLimits,
) -> Result<Vec<Polynomial<PrimeField>>, RunError> {
    let ring = PolynomialRing::prime_field(prime, input.ring.variables())
        .expect("a prime of the sequence is a ring modulus");
    let mut generators: Vec<Polynomial<PrimeField>> = Vec::with_capacity(input.cleared.len());
    for cleared in &input.cleared {
        if let Some(stop) = limits.stop() {
            return Err(stop);
        }
        generators.push(cleared.image(&ring));
    }
    let options = ComputeOptions::new().backend(input.backend);
    compute::groebner_basis_with_limits(&ring, &generators, &options, limits)
}

/// Report whether the driver runs `prime`.
///
/// A prime that divides the leading coefficient of a cleared generator is
/// skipped before any run starts, because the image of that generator
/// leads on another monomial.
fn usable(prime: u64, cleared: &[ClearedPolynomial]) -> bool {
    cleared.iter().all(|f| match f.leading_coefficient() {
        Some(coeff) => residue(coeff, prime) != 0,
        None => true,
    })
}
