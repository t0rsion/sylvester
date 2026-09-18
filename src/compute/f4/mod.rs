//! The F4 engine (design section 3).
//!
//! [`groebner_basis`] is the whole engine: it interns the generators, runs
//! the batch loop, interreduces, and converts the result back to
//! [`Polynomial`]. Everything inside the run works on interned monomials
//! ([`monomial`]), on column indices ([`matrix`]), and on one field
//! context ([`field`]).
//!
//! The lane width comes from the input: a run starts at the narrowest
//! packing the generators' exponents fit and starts again at the next
//! width on [`F4Error::LaneOverflow`] (design section 2.3). The monomial
//! order is grevlex, sealed.

// The interfaces of design section 9 are frozen, so every module carries
// the whole API the later chunks compile against, not only what the run
// loop calls today.
#![allow(dead_code)]

pub(crate) mod basis;
pub(crate) mod field;
pub(crate) mod kernel;
pub(crate) mod matrix;
pub(crate) mod monomial;
pub(crate) mod pairs;
pub(crate) mod symbolic;
pub(crate) mod trace;

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use basis::{Basis, BasisPoly, InputPoly};
use field::{FieldOps, Small31};
use matrix::{Batch, Workspace};
use monomial::{Lanes, Lanes8, Lanes16, Lanes32, MAX_EXPONENT, MonomialId, MonomialTable};
use pairs::{PairCounters, PairSet, SelectOptions};
use symbolic::{Strategy, Symbolic};
use trace::{BatchRow, BatchRows, Insertion, NoTrace, Returned, RowId, Trace};

use super::{ComputeError, ComputeLimits, DEGREE_LIMIT, RunError};
use crate::poly::{Exps, Monomial, Polynomial, Term};
use crate::ring::PolynomialRing;
use crate::ring::field::Felt;

/// Why an F4 run stops before it has a basis.
///
/// [`F4Error::LaneOverflow`] is internal. The run loop catches it, discards
/// the run state, and starts again at the next lane width (design section
/// 2.3). [`F4Error::ExponentOverflow`] surfaces as
/// [`ComputeError::ExponentLimit`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum F4Error {
    /// The deadline passed.
    Timeout,
    /// The live engine data passed the memory limit.
    MemoryLimitExceeded,
    /// A result passed the current lane width's bound but not
    /// [`MAX_EXPONENT`]. The run restarts at the next width.
    LaneOverflow,
    /// A result passed [`MAX_EXPONENT`]. No width holds it.
    ExponentOverflow,
    /// A table reached `u32::MAX - 1` monomials.
    TableFull,
    /// The cancellation flag of the run was set.
    Cancelled,
}

/// The number of items one loop takes between two reads of the clock.
///
/// Reading the clock costs more than one pass of the loops that call
/// [`Deadline::tick`], so a loop counts its items and reads the clock
/// every [`TICK`] of them.
const TICK: u32 = 64;

/// The deadline and cancellation flags of one run, checked at a fixed
/// cadence.
///
/// [`Deadline::tick`] is what a loop calls per item. It reads the clock
/// once every [`TICK`] calls. [`Deadline::check`] reads it at once and is
/// what the run loop calls between batches. A cancelled run stops at the
/// same cadence, so either flag stops it within [`TICK`] items of being set.
pub(crate) struct Deadline {
    at: Option<Instant>,
    cancel: Option<Arc<AtomicBool>>,
    external_cancel: Option<Arc<AtomicBool>>,
    ticks: u32,
}

impl Deadline {
    /// The deadline and cancellation flags of `limits`.
    pub(crate) fn of(limits: &ComputeLimits) -> Self {
        Deadline {
            at: limits.deadline,
            cancel: limits.cancel.clone(),
            external_cancel: limits.external_cancel.clone(),
            ticks: 0,
        }
    }

    /// A deadline that never passes.
    #[cfg(test)]
    pub(crate) fn none() -> Self {
        Deadline {
            at: None,
            cancel: None,
            external_cancel: None,
            ticks: 0,
        }
    }

    /// Count one item, and check the clock every [`TICK`] items.
    #[inline]
    pub(crate) fn tick(&mut self) -> Result<(), F4Error> {
        if self.at.is_none() && self.cancel.is_none() && self.external_cancel.is_none() {
            return Ok(());
        }
        self.ticks += 1;
        if self.ticks < TICK {
            return Ok(());
        }
        self.ticks = 0;
        self.check()
    }

    /// A deadline at the same instant, with the same cancellation flags.
    ///
    /// The parallel phase of the kernel builds one per worker, because a
    /// worker counts its own items. A worker reads both flags, so cancelling
    /// a run stops it inside a parallel batch too.
    pub(crate) fn fork(&self) -> Self {
        Deadline {
            at: self.at,
            cancel: self.cancel.clone(),
            external_cancel: self.external_cancel.clone(),
            ticks: 0,
        }
    }

    /// Check the flags and the clock now.
    pub(crate) fn check(&mut self) -> Result<(), F4Error> {
        let cancelled = self
            .cancel
            .as_ref()
            .is_some_and(|flag| flag.load(Ordering::Relaxed))
            || self
                .external_cancel
                .as_ref()
                .is_some_and(|flag| flag.load(Ordering::Relaxed));
        if cancelled {
            return Err(F4Error::Cancelled);
        }
        match self.at {
            Some(at) if Instant::now() >= at => Err(F4Error::Timeout),
            _ => Ok(()),
        }
    }
}

/// The internal settings of one F4 run.
///
/// None of these is public (design open decision 7). `max_batch_pairs`
/// defaults to no cap, and the strategy defaults are the M1 baseline of
/// design section 3.5.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct F4Options {
    /// The largest number of pairs one batch takes, or `None` for no cap.
    pub(crate) max_batch_pairs: Option<usize>,
    /// The reducer and upper-row choices of symbolic preprocessing.
    pub(crate) strategy: Strategy,
    /// The thread count the caller asked for, or `None` for the rayon
    /// global pool. A count of 1 keeps every row of a batch on the calling
    /// thread.
    pub(crate) threads: Option<usize>,
}

/// What one run did, for the benchmark harness and the tests.
///
/// The counters are a report, never an input to a decision. Design section
/// 8.6 lists them.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct RunCounters {
    /// The batches the run reduced, including the final interreduction.
    pub(crate) batches: u64,
    /// The rows of every batch.
    pub(crate) matrix_rows: u64,
    /// The columns of every batch.
    pub(crate) matrix_columns: u64,
    /// The nonzero entries of every batch, before reduction.
    pub(crate) matrix_nonzeros: u64,
    /// The rows that reduced to zero.
    pub(crate) zero_rows: u64,
    /// The pivots the reduction installed.
    pub(crate) new_pivots: u64,
    /// The lane restarts of design section 2.3.
    pub(crate) lane_restarts: u32,
    /// The monomials of the basis table at the end of the run.
    pub(crate) basis_monomials: usize,
    /// The batches the run built, found too large for the budget, and
    /// built again with half the pairs.
    pub(crate) batch_retries: u64,
    /// The counters of the pair set.
    pub(crate) pairs: PairCounters,
    /// The pool the run computed in, or `None` for the rayon global pool.
    pub(crate) threads: Option<usize>,
}

/// Compute the reduced Gröbner basis with the F4 engine.
///
/// The result is monic and sorted strictly descending by leading monomial.
/// `limits` carries the deadline, the memory limit, the thread count, and
/// the cancellation flag.
pub(crate) fn groebner_basis(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    limits: &ComputeLimits,
) -> Result<Vec<Polynomial>, RunError> {
    solve(ring, generators, limits).map(|(basis, _)| basis)
}

/// Compute the basis and report the run counters.
pub(crate) fn solve(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    limits: &ComputeLimits,
) -> Result<(Vec<Polynomial>, RunCounters), RunError> {
    solve_recorded(ring, generators, limits, &mut NoTrace)
}

/// Compute the basis and report every operation to `trace`.
///
/// A run that records stays on the calling thread and builds no pool,
/// because the parallel phase of the kernel reports nothing. The reports
/// are the ones `docs/certificate-v2.md` section 9.4 fixes.
///
/// A thread count above 1 runs the whole computation inside a rayon pool
/// of that size, so the kernel's parallel phase gets the threads the
/// caller asked for and no others. Building the pool is best effort: a
/// machine that refuses it leaves the run on the calling thread, and
/// `RunCounters::threads` then reports 1. Neither the pool nor the thread
/// count changes the basis (design section 5).
pub(crate) fn solve_recorded<T: Trace + Send>(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    limits: &ComputeLimits,
    trace: &mut T,
) -> Result<(Vec<Polynomial>, RunCounters), RunError> {
    let mut f4 = F4Options {
        threads: limits.threads,
        ..F4Options::default()
    };
    if !T::RECORDS
        && let Some(count) = limits.threads
        && count > 1
    {
        match rayon::ThreadPoolBuilder::new().num_threads(count).build() {
            Ok(pool) => {
                let (basis, mut counters) = pool
                    .install(|| dispatch(ring, generators, &f4, limits, trace))
                    .map_err(report)?;
                counters.threads = Some(count);
                return Ok((basis, counters));
            }
            Err(_) => f4.threads = Some(1),
        }
    }
    let (basis, mut counters) = dispatch(ring, generators, &f4, limits, trace).map_err(report)?;
    counters.threads = if T::RECORDS { Some(1) } else { f4.threads };
    Ok((basis, counters))
}

/// Map an engine stop onto the report the crate passes on.
fn report(error: F4Error) -> RunError {
    let compute = match error {
        F4Error::Timeout => ComputeError::Timeout,
        F4Error::MemoryLimitExceeded => ComputeError::MemoryLimitExceeded,
        F4Error::ExponentOverflow => ComputeError::ExponentLimit {
            limit: DEGREE_LIMIT,
        },
        F4Error::TableFull => ComputeError::TableFull,
        // The run loop catches this one and starts again at the next
        // width, so it reaches here only from the widest width, where no
        // width is left to hold the value.
        F4Error::LaneOverflow => ComputeError::ExponentLimit {
            limit: DEGREE_LIMIT,
        },
        F4Error::Cancelled => return RunError::Cancelled,
    };
    RunError::Compute(compute)
}

/// One lane packing of the monomial tables.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Width {
    /// [`Lanes8`], exponents up to 127.
    W8,
    /// [`Lanes16`], exponents up to 32767.
    W16,
    /// [`Lanes32`], exponents up to 65535.
    W32,
}

impl Width {
    /// The narrowest width whose bound `max_exponent` satisfies.
    fn for_exponent(max_exponent: u32) -> Self {
        if max_exponent <= Lanes8::BOUND {
            Width::W8
        } else if max_exponent <= Lanes16::BOUND {
            Width::W16
        } else {
            Width::W32
        }
    }

    /// The next wider packing, or `None` for the widest.
    fn wider(self) -> Option<Self> {
        match self {
            Width::W8 => Some(Width::W16),
            Width::W16 => Some(Width::W32),
            Width::W32 => None,
        }
    }
}

/// Run at the width the input asks for, and restart wider on overflow.
///
/// A restart discards every run structure and starts again from the
/// generators, because promoting an append-only table in place would
/// rebuild every key, every hash, and every hash-table entry (design
/// section 2.3). At most two restarts happen.
fn dispatch<T: Trace>(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &F4Options,
    limits: &ComputeLimits,
    trace: &mut T,
) -> Result<(Vec<Polynomial>, RunCounters), F4Error> {
    // A public exponent is a u16, so it is at most MAX_EXPONENT and the
    // widest packing holds every input this crate builds.
    let max_exponent = max_exponent(generators);
    debug_assert!(max_exponent <= MAX_EXPONENT, "an exponent came from a u16");
    let mut width = Width::for_exponent(max_exponent);
    let mut restarts = 0;
    loop {
        let attempt = match width {
            Width::W8 => start::<Lanes8, T>(ring, generators, options, limits, max_exponent, trace),
            Width::W16 => {
                start::<Lanes16, T>(ring, generators, options, limits, max_exponent, trace)
            }
            Width::W32 => {
                start::<Lanes32, T>(ring, generators, options, limits, max_exponent, trace)
            }
        };
        match attempt {
            Err(F4Error::LaneOverflow) => match width.wider() {
                Some(next) => {
                    width = next;
                    restarts += 1;
                    // The attempt that overflowed wrote nodes the run
                    // discards with the rest of its state (contract
                    // section 9.2).
                    trace.restart();
                }
                // The widest packing holds every exponent the public
                // contract allows, so a lane overflow there says the value
                // is past the contract.
                None => return Err(F4Error::ExponentOverflow),
            },
            other => {
                return other.map(|(basis, mut counters)| {
                    counters.lane_restarts = restarts;
                    (basis, counters)
                });
            }
        }
    }
}

/// Build the field context and run at one lane width.
fn start<L: Lanes, T: Trace>(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    options: &F4Options,
    limits: &ComputeLimits,
    max_exponent: u32,
    trace: &mut T,
) -> Result<(Vec<Polynomial>, RunCounters), F4Error> {
    // The ring caps the modulus at 2^31 - 1, so the kernel takes it.
    let p = ring.modulus() as u32;
    let field = Small31::new(p);
    run::<L, _, _>(
        ring,
        generators,
        &field,
        options,
        limits,
        max_exponent,
        trace,
    )
}

/// One run at one lane width, with one field kernel (design section 3.2).
fn run<L: Lanes, F: FieldOps<Coeff = u32> + Sync, T: Trace>(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    field: &F,
    options: &F4Options,
    limits: &ComputeLimits,
    max_exponent: u32,
    trace: &mut T,
) -> Result<(Vec<Polynomial>, RunCounters), F4Error> {
    let nvars = ring.nvars();
    let mut clock = Deadline::of(limits);
    clock.check()?;

    if generators.is_empty() {
        if T::RECORDS {
            trace.returned(&[]);
        }
        return Ok((Vec::new(), RunCounters::default()));
    }

    let mut state = F4State::<L>::new(
        generators,
        field,
        RunSetup {
            nvars,
            max_exponent,
            limits,
        },
        trace,
        &mut clock,
    )?;
    state.process_batches(field, options, limits, trace, &mut clock)?;
    state.finish(ring, field, options, limits, trace, &mut clock)
}

struct F4State<L: Lanes> {
    basis: Basis<L>,
    update_ws: MonomialTable<L>,
    sym_table: MonomialTable<L>,
    pairs: PairSet,
    batch: Batch,
    ws: Workspace,
    counters: RunCounters,
    monos: Vec<MonomialId>,
    coeffs: Vec<u32>,
    report: Report,
    seed_unit: bool,
    nvars: usize,
}

struct RunSetup<'a> {
    nvars: usize,
    max_exponent: u32,
    limits: &'a ComputeLimits,
}

impl<L: Lanes> F4State<L> {
    fn new<F: FieldOps<Coeff = u32>, T: Trace>(
        generators: &[Polynomial],
        field: &F,
        setup: RunSetup<'_>,
        trace: &T,
        clock: &mut Deadline,
    ) -> Result<Self, F4Error> {
        let mut table = MonomialTable::<L>::with_max_exponent(setup.nvars, setup.max_exponent)
            .with_limit(setup.limits.memory);
        let inputs = intern_generators(generators, &mut table)?;
        let mut state = F4State {
            basis: Basis::new(table),
            update_ws: MonomialTable::<L>::with_max_exponent(setup.nvars, setup.max_exponent)
                .with_limit(setup.limits.memory),
            sym_table: MonomialTable::<L>::with_max_exponent(setup.nvars, setup.max_exponent)
                .with_limit(setup.limits.memory),
            pairs: PairSet::new(),
            batch: Batch::default(),
            ws: Workspace::default(),
            counters: RunCounters::default(),
            monos: Vec::new(),
            coeffs: Vec::new(),
            report: Report::default(),
            seed_unit: false,
            nvars: setup.nvars,
        };
        state
            .basis
            .seed(field, inputs, &mut state.pairs, &mut state.update_ws, clock)?;
        state.seed_unit = state.basis.is_unit();
        state.check_memory(setup.limits.memory, trace)?;
        Ok(state)
    }

    fn check_memory<T: Trace>(&self, memory: Option<usize>, trace: &T) -> Result<(), F4Error> {
        check_memory(
            held_bytes(
                &self.basis,
                &self.pairs,
                &self.update_ws,
                &self.sym_table,
                &self.batch,
                &self.ws,
                trace,
            ),
            memory,
        )
    }

    fn process_batches<F, T>(
        &mut self,
        field: &F,
        options: &F4Options,
        limits: &ComputeLimits,
        trace: &mut T,
        clock: &mut Deadline,
    ) -> Result<(), F4Error>
    where
        F: FieldOps<Coeff = u32> + Sync,
        T: Trace,
    {
        let select = SelectOptions {
            max_batch_pairs: options.max_batch_pairs,
        };
        while !self.basis.is_unit() {
            clock.check()?;
            let Some(selected) = self.pairs.take_lowest_degree(self.basis.table(), &select) else {
                break;
            };
            self.process_batch(selected, field, options, limits, trace, clock)?;
        }
        Ok(())
    }

    fn process_batch<F, T>(
        &mut self,
        selected: Vec<pairs::Pair>,
        field: &F,
        options: &F4Options,
        limits: &ComputeLimits,
        trace: &mut T,
        clock: &mut Deadline,
    ) -> Result<(), F4Error>
    where
        F: FieldOps<Coeff = u32> + Sync,
        T: Trace,
    {
        self.build_batch(selected, options, limits, trace, clock)?;
        self.count_batch();
        if T::RECORDS {
            self.report.send_rows(self.nvars, trace);
        }
        let mut context = kernel::ReductionContext::new(
            &self.basis,
            field,
            trace,
            clock,
            limits.memory,
            options.threads,
        );
        let reduction = kernel::reduce(&mut self.batch, &mut self.ws, &mut context)?;
        self.counters.zero_rows += u64::from(reduction.zero_rows);
        self.counters.new_pivots += u64::from(reduction.new_pivots);
        self.insert_pivots(field, trace, clock)?;
        self.check_memory(limits.memory, trace)
    }

    fn build_batch<T: Trace>(
        &mut self,
        mut selected: Vec<pairs::Pair>,
        options: &F4Options,
        limits: &ComputeLimits,
        trace: &T,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        loop {
            self.build_matrix::<T>(&selected, options, clock)?;
            if within(self.memory_bytes(trace), limits.memory) {
                return Ok(());
            }
            self.retry_batch(&mut selected)?;
        }
    }

    fn build_matrix<T: Trace>(
        &mut self,
        selected: &[pairs::Pair],
        options: &F4Options,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        self.sym_table.clear();
        let mut sym = symbolic::preprocess(
            selected,
            &self.basis,
            &mut self.sym_table,
            &options.strategy,
            clock,
        )?;
        matrix::build(&mut sym, &self.basis, &mut self.batch, &mut self.ws, clock)?;
        if T::RECORDS {
            self.report.collect_rows(&sym, &self.batch);
        }
        Ok(())
    }

    fn retry_batch(&mut self, selected: &mut Vec<pairs::Pair>) -> Result<(), F4Error> {
        if selected.len() < 2 {
            return Err(F4Error::MemoryLimitExceeded);
        }
        let half = selected.len() / 2;
        self.pairs.requeue(selected.split_off(half))?;
        self.sym_table.clear();
        self.sym_table.shrink();
        self.batch.clear();
        self.batch.shrink();
        self.ws.shrink();
        self.counters.batch_retries += 1;
        Ok(())
    }

    fn memory_bytes<T: Trace>(&self, trace: &T) -> usize {
        held_bytes(
            &self.basis,
            &self.pairs,
            &self.update_ws,
            &self.sym_table,
            &self.batch,
            &self.ws,
            trace,
        )
    }

    fn count_batch(&mut self) {
        self.counters.batches += 1;
        self.counters.matrix_rows +=
            self.batch.upper().len() as u64 + self.batch.lower().len() as u64;
        self.counters.matrix_columns += self.batch.ncols() as u64;
        self.counters.matrix_nonzeros += nonzeros(&self.batch);
    }

    fn insert_pivots<F: FieldOps<Coeff = u32>, T: Trace>(
        &mut self,
        field: &F,
        trace: &mut T,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        let fresh = self.adopt_pivots(field)?;
        let slots = self.batch.npiv() + self.batch.new_pivots().len() as u32;
        self.basis.insert_batch(
            fresh,
            &mut self.pairs,
            &mut self.update_ws,
            clock,
            T::RECORDS.then_some(&mut self.report.assigned),
        )?;
        if T::RECORDS {
            self.report.send_insertions(slots, trace);
        }
        Ok(())
    }

    fn adopt_pivots<F: FieldOps<Coeff = u32>>(
        &mut self,
        field: &F,
    ) -> Result<Vec<BasisPoly>, F4Error> {
        let mut fresh = Vec::new();
        basis::reserve(&mut fresh, self.batch.new_pivots().len())?;
        for row in self.batch.new_pivots().iter().rev() {
            self.monos.clear();
            self.coeffs.clear();
            let cols = self.batch.row_columns(row);
            self.monos
                .extend(cols.iter().map(|&col| self.batch.column(col)));
            let (values, _) = self.batch.row_coeffs(row, &self.basis);
            self.coeffs.extend_from_slice(values);
            fresh.push(
                self.basis
                    .adopt_row(&self.sym_table, &self.monos, &self.coeffs, field)?,
            );
        }
        Ok(fresh)
    }

    fn finish<F, T>(
        mut self,
        ring: &PolynomialRing,
        field: &F,
        options: &F4Options,
        limits: &ComputeLimits,
        trace: &mut T,
        clock: &mut Deadline,
    ) -> Result<(Vec<Polynomial>, RunCounters), F4Error>
    where
        F: FieldOps<Coeff = u32>,
        T: Trace,
    {
        self.counters.pairs = self.pairs.counters();
        self.counters.basis_monomials = self.basis.table().len();
        self.check_memory(limits.memory, trace)?;
        if let Some(basis) = self.terminal_basis(ring, trace) {
            return Ok((basis, self.counters));
        }
        let basis = self.interreduce(ring, field, options, limits, trace, clock)?;
        Ok((basis, self.counters))
    }

    fn terminal_basis<T: Trace>(
        &self,
        ring: &PolynomialRing,
        trace: &mut T,
    ) -> Option<Vec<Polynomial>> {
        if self.basis.is_unit() {
            if T::RECORDS {
                trace.returned(&[unit_element(&self.basis, &self.batch, self.seed_unit)]);
            }
            return Some(vec![ring.one()]);
        }
        if self.basis.live().is_empty() {
            if T::RECORDS {
                trace.returned(&[]);
            }
            return Some(Vec::new());
        }
        None
    }

    fn interreduce<F, T>(
        &mut self,
        ring: &PolynomialRing,
        field: &F,
        options: &F4Options,
        limits: &ComputeLimits,
        trace: &mut T,
        clock: &mut Deadline,
    ) -> Result<Vec<Polynomial>, F4Error>
    where
        F: FieldOps<Coeff = u32>,
        T: Trace,
    {
        clock.check()?;
        self.build_interreduction::<T>(options, clock)?;
        if T::RECORDS {
            self.report.send_rows(self.nvars, trace);
        }
        self.counters.batches += 1;
        self.counters.matrix_rows += self.batch.upper().len() as u64;
        self.counters.matrix_columns += self.batch.ncols() as u64;
        self.counters.matrix_nonzeros += nonzeros(&self.batch);
        kernel::interreduce(
            &mut self.batch,
            &self.basis,
            field,
            trace,
            &mut self.ws,
            clock,
        )?;
        self.check_memory(limits.memory, trace)?;
        if T::RECORDS {
            self.report.send_returned(&self.basis, &self.batch, trace);
        }
        Ok(convert_out(ring, &self.basis, &self.batch, &self.sym_table))
    }

    fn build_interreduction<T: Trace>(
        &mut self,
        options: &F4Options,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        self.sym_table.clear();
        let mut sym =
            symbolic::preprocess_basis(&self.basis, &mut self.sym_table, &options.strategy, clock)?;
        matrix::build(&mut sym, &self.basis, &mut self.batch, &mut self.ws, clock)?;
        if T::RECORDS {
            self.report.collect_rows(&sym, &self.batch);
        }
        Ok(())
    }
}

/// The buffers the four reports of contract section 9.4 are built in.
///
/// The run holds one value, so a report costs no allocation after the
/// first batch that needs one.
#[derive(Default)]
struct Report {
    rows: Vec<BatchRow>,
    mults: Vec<u32>,
    exps: Vec<u32>,
    assigned: Vec<Option<u32>>,
    insertions: Vec<Insertion>,
    returned: Vec<Returned>,
}

impl Report {
    /// Collect the rows of the batch that `build` just wrote (report 1).
    ///
    /// The batch keeps no multiplier and no source, and it filters and
    /// sorts the rows, so the entries come from the symbolic rows the
    /// batch names.
    fn collect_rows<L: Lanes>(&mut self, sym: &Symbolic<'_, L>, batch: &Batch) {
        self.rows.clear();
        self.mults.clear();
        let pivots = batch
            .upper()
            .iter()
            .enumerate()
            .map(|(slot, row)| (RowId::Pivot(slot as u32), row));
        let lower = batch
            .lower()
            .iter()
            .enumerate()
            .map(|(slot, row)| (RowId::Lower(slot as u32), row));
        for (place, row) in pivots.chain(lower) {
            let raw = sym.rows[row.raw() as usize];
            self.rows.push(BatchRow {
                place,
                source: raw.source,
            });
            sym.table.unpack(raw.mult, &mut self.exps);
            self.mults.extend_from_slice(&self.exps);
        }
    }

    /// Report the rows of the batch that is about to reduce.
    fn send_rows<T: Trace>(&self, nvars: usize, trace: &mut T) {
        trace.rows(&BatchRows {
            rows: &self.rows,
            mults: &self.mults,
            nvars,
        });
    }

    /// Report the elements the batch added to the basis (report 3).
    ///
    /// `slots` is the number of pivot slots of the batch. Insertion takes
    /// the new pivots in decreasing slot, so entry `k` is slot
    /// `slots - 1 - k`.
    fn send_insertions<T: Trace>(&mut self, slots: u32, trace: &mut T) {
        self.insertions.clear();
        self.insertions.extend(
            self.assigned
                .iter()
                .enumerate()
                .map(|(k, &basis)| Insertion {
                    pivot: slots - 1 - k as u32,
                    basis,
                }),
        );
        trace.inserted(&self.insertions);
    }

    /// Report the basis the run returns (report 4).
    ///
    /// The entries are the pivot slots of the interreduction batch that
    /// hold a live basis element, in the order [`convert_out`] reads them.
    fn send_returned<L: Lanes, T: Trace>(
        &mut self,
        basis: &Basis<L>,
        batch: &Batch,
        trace: &mut T,
    ) {
        let elements = basis.live().len() as u32;
        self.returned.clear();
        self.returned.extend(
            batch
                .upper()
                .iter()
                .enumerate()
                .filter(|(_, row)| row.raw() < elements)
                .map(|(slot, _)| Returned::Pivot(slot as u32)),
        );
        trace.returned(&self.returned);
    }
}

/// The one element of a unit basis (report 4).
///
/// Seeding stops at the first constant generator, and that generator is
/// the live element. A batch instead installs the constant as the pivot of
/// the identity column.
fn unit_element<L: Lanes>(basis: &Basis<L>, batch: &Batch, seed_unit: bool) -> Returned {
    if seed_unit {
        let element = basis.live()[0];
        let input = basis
            .source(element)
            .expect("a seeded element carries its input index");
        return Returned::Input(input);
    }
    let slot = (0..batch.pivot_count())
        .find(|&slot| batch.column(batch.pivot(slot).lead()) == MonomialId::ONE)
        .expect("the batch installed the constant pivot that took the unit basis");
    Returned::Pivot(slot)
}

/// The nonzero entries of every row of one batch.
fn nonzeros(batch: &Batch) -> u64 {
    let rows = batch.upper().iter().chain(batch.lower());
    rows.map(|row| u64::from(row.len())).sum()
}

/// The bytes the run holds, counting capacities (design section 3.9).
///
/// The trace is part of the run, so a recorded run charges its nodes to
/// the same limit the engine holds to.
fn held_bytes<L: Lanes, T: Trace>(
    basis: &Basis<L>,
    pairs: &PairSet,
    update_ws: &MonomialTable<L>,
    sym_table: &MonomialTable<L>,
    batch: &Batch,
    ws: &Workspace,
    trace: &T,
) -> usize {
    basis
        .memory_bytes()
        .saturating_add(pairs.memory_bytes())
        .saturating_add(update_ws.memory_bytes())
        .saturating_add(sym_table.memory_bytes())
        .saturating_add(batch.bytes())
        .saturating_add(ws.bytes())
        .saturating_add(trace.held_bytes())
}

/// Report whether `held` bytes fit the budget.
fn within(held: usize, memory: Option<usize>) -> bool {
    memory.is_none_or(|limit| held <= limit)
}

/// Stop the run when `held` bytes pass the budget.
fn check_memory(held: usize, memory: Option<usize>) -> Result<(), F4Error> {
    if within(held, memory) {
        return Ok(());
    }
    Err(F4Error::MemoryLimitExceeded)
}

/// The largest exponent any generator holds.
fn max_exponent(generators: &[Polynomial]) -> u32 {
    generators
        .iter()
        .flat_map(|poly| poly.terms())
        .flat_map(|(_, exps)| exps.iter().copied())
        .map(u32::from)
        .max()
        .unwrap_or(0)
}

/// Intern the generators into the basis table (design section 2.1).
///
/// One pass, in the caller's order. The index of a generator in the result
/// is its public input index, which [`Basis::seed`] records.
fn intern_generators<L: Lanes>(
    generators: &[Polynomial],
    table: &mut MonomialTable<L>,
) -> Result<Vec<InputPoly>, F4Error> {
    let mut inputs = Vec::new();
    basis::reserve(&mut inputs, generators.len())?;
    let mut exps: Vec<u32> = Vec::new();
    for poly in generators {
        let mut monos = Vec::new();
        let mut coeffs = Vec::new();
        basis::reserve(&mut monos, poly.terms().len())?;
        basis::reserve(&mut coeffs, poly.terms().len())?;
        for (coeff, term) in poly.terms() {
            exps.clear();
            exps.extend(term.iter().map(|&e| u32::from(e)));
            monos.push(table.intern(&exps)?);
            // The ring caps the modulus at 2^31 - 1, so a residue fits a
            // u32.
            coeffs.push(coeff.value() as u32);
        }
        inputs.push(InputPoly { monos, coeffs });
    }
    Ok(inputs)
}

/// Convert the interreduced basis back to [`Polynomial`] (design section
/// 2.1).
///
/// The rows of `batch` are the live basis elements and the reducer
/// multiples symbolic closure added. A row whose `raw` index is below the
/// number of live elements is an element; the rest are dropped. The rows
/// are sorted by leading column ascending, so the result is sorted by
/// leading monomial descending.
fn convert_out<L: Lanes>(
    ring: &PolynomialRing,
    basis: &Basis<L>,
    batch: &Batch,
    table: &MonomialTable<L>,
) -> Vec<Polynomial> {
    let elements = basis.live().len() as u32;
    let mut out = Vec::with_capacity(elements as usize);
    let mut exps: Vec<u32> = Vec::new();
    for row in batch.upper() {
        if row.raw() >= elements {
            continue;
        }
        let cols = batch.row_columns(row);
        let (vals, _) = batch.row_coeffs(row, basis);
        let mut terms: Vec<Term> = Vec::with_capacity(cols.len());
        // A row's columns ascend, so its monomials descend. A polynomial
        // holds its terms ascending.
        for (&col, &val) in cols.iter().zip(vals).rev() {
            table.unpack(batch.column(col), &mut exps);
            let mono: Exps = exps.iter().map(|&e| e as u16).collect();
            terms.push(Term {
                coeff: Felt::from_residue(u64::from(val)),
                mono: Monomial::from_exps(mono),
            });
        }
        out.push(Polynomial::from_sorted_terms(ring.clone(), terms));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring(nvars: usize) -> PolynomialRing {
        let names: Vec<String> = (0..nvars).map(|i| format!("x{i}")).collect();
        PolynomialRing::prime_field(32003, names).expect("32003 is prime")
    }

    fn parse(ring: &PolynomialRing, text: &str) -> Polynomial {
        ring.parse_polynomial(text).expect("the syntax holds")
    }

    #[test]
    fn the_run_counters_report_the_work() {
        let ring = ring(3);
        let generators = vec![
            parse(&ring, "x0^2 + x1*x2"),
            parse(&ring, "x0*x1 - x2^2"),
            parse(&ring, "x1^3 - x0*x2"),
        ];
        let (basis, counters) =
            solve(&ring, &generators, &ComputeLimits::default()).expect("no limit");

        assert!(!basis.is_empty());
        assert!(counters.batches > 1, "{counters:?}");
        assert!(counters.matrix_rows >= counters.batches, "{counters:?}");
        assert!(counters.matrix_columns > 0, "{counters:?}");
        assert!(
            counters.matrix_nonzeros >= counters.matrix_rows,
            "{counters:?}"
        );
        assert!(counters.new_pivots > 0, "{counters:?}");
        assert!(counters.pairs.generated > 0, "{counters:?}");
        assert!(counters.basis_monomials > 1, "{counters:?}");
        assert_eq!(counters.lane_restarts, 0, "{counters:?}");
    }

    #[test]
    fn a_wide_exponent_restarts_the_run_at_the_next_lane_width() {
        // x0^100 and x1^100 multiply to degree 200 in one variable pair,
        // which passes the 127 bound of Lanes8 during the run.
        let ring = ring(2);
        let generators = vec![
            parse(&ring, "x0^100*x1^30 - 1"),
            parse(&ring, "x0^30*x1^100 - 1"),
        ];
        let (basis, counters) =
            solve(&ring, &generators, &ComputeLimits::default()).expect("no limit");
        assert!(!basis.is_empty());
        assert_eq!(counters.lane_restarts, 1, "{counters:?}");
    }

    #[test]
    fn the_starting_width_follows_the_input_exponents() {
        assert_eq!(Width::for_exponent(0), Width::W8);
        assert_eq!(Width::for_exponent(127), Width::W8);
        assert_eq!(Width::for_exponent(128), Width::W16);
        assert_eq!(Width::for_exponent(32767), Width::W16);
        assert_eq!(Width::for_exponent(32768), Width::W32);
        assert_eq!(Width::for_exponent(MAX_EXPONENT), Width::W32);
    }

    #[test]
    fn an_external_cancellation_flag_reaches_every_f4_deadline() {
        let flag = Arc::new(AtomicBool::new(false));
        let limits = ComputeLimits {
            external_cancel: Some(Arc::clone(&flag)),
            ..ComputeLimits::default()
        };
        let mut clock = Deadline::of(&limits);
        for _ in 0..TICK {
            assert_eq!(clock.tick(), Ok(()));
        }
        flag.store(true, Ordering::Relaxed);
        for _ in 0..TICK - 1 {
            assert_eq!(clock.tick(), Ok(()));
        }
        assert_eq!(clock.tick(), Err(F4Error::Cancelled));
        let mut fork = clock.fork();
        assert_eq!(fork.check(), Err(F4Error::Cancelled));
    }
}
