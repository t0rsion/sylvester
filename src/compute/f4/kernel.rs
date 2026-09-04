//! The reduction kernel and new pivot discovery (design section 3.7).
//!
//! One row reduces in three steps. It scatters into a dense `u64`
//! accumulator, one lane per column. It walks the columns in increasing
//! order and, at every column that has a pivot row, adds `p - v` times
//! that pivot's tail. It gathers what is left. The walk covers every
//! column of the row, and it sets every lane it reads to physical zero, so
//! the accumulator is clean for the next row without a clearing pass.
//!
//! The kernel reads column indices and coefficients. It never touches a
//! monomial table.

use std::sync::atomic::{AtomicBool, Ordering};

use super::field::FieldOps;
use super::matrix::{BasisCoeffs, Batch, NO_ROW, PartialRow, Workspace, fill, reserve};
use super::monomial::MonomialId;
use super::trace::{NoTrace, RowId, Trace};
use super::{Deadline, F4Error};

/// What one call to [`reduce`] produced.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Reduction {
    /// The number of pivots installed, all in [`Batch::new_pivots`].
    pub(crate) new_pivots: u32,
    /// The number of rows that reduced to zero.
    pub(crate) zero_rows: u32,
    /// A new pivot is a nonzero constant, so the ideal is the whole ring.
    ///
    /// The run stops at once and returns the unit basis (design section
    /// 3.2). The rows after the constant are not reduced.
    pub(crate) unit: bool,
}

/// The accumulator, the pivot map, and the row the kernel is writing.
struct Scratch<'a> {
    acc: &'a mut [u64],
    pivot_at: &'a mut [u32],
    cols: &'a mut Vec<u32>,
    vals: &'a mut Vec<u32>,
    shoup: &'a mut Vec<u64>,
}

pub(crate) struct ReductionContext<'a, F, B, T> {
    basis: &'a B,
    field: &'a F,
    trace: &'a mut T,
    clock: &'a mut Deadline,
    memory: Option<usize>,
    threads: Option<usize>,
}

impl<'a, F, B, T> ReductionContext<'a, F, B, T> {
    pub(crate) fn new(
        basis: &'a B,
        field: &'a F,
        trace: &'a mut T,
        clock: &'a mut Deadline,
        memory: Option<usize>,
        threads: Option<usize>,
    ) -> Self {
        Self {
            basis,
            field,
            trace,
            clock,
            memory,
            threads,
        }
    }
}

enum LoadedRow {
    Zero,
    Start(usize),
}

/// Reduce the `[C | D]` rows of `batch` and install the new pivots.
///
/// Rows reduce in batch order, and a row that survives becomes the pivot
/// of its leading column at once, so every later row reduces against it
/// (design section 3.7). After the row loop the new pivots are
/// interreduced backward, and [`Batch::new_pivots`] holds them sorted by
/// leading column.
///
/// `threads` is the count the caller asked for, or `None` for the rayon
/// global pool. A count of 1 keeps every row on the calling thread.
///
/// Both field kernels carry `Coeff = u32`. The kernel needs a coefficient
/// multiply for the monic normalization, which it does as one
/// [`FieldOps::reduce_acc`] of the `u64` product, so it names the concrete
/// coefficient type.
pub(crate) fn reduce<F, B, T>(
    batch: &mut Batch,
    ws: &mut Workspace,
    context: &mut ReductionContext<'_, F, B, T>,
) -> Result<Reduction, F4Error>
where
    F: FieldOps<Coeff = u32> + Sync,
    B: BasisCoeffs + Sync,
    T: Trace,
{
    let ncols = batch.ncols();
    prepare(batch, ws, ncols)?;
    batch.debug_check_coeffs(context.basis);
    let npiv = batch.npiv() as usize;
    // The parallel phase writes here, and taking the vector out of the
    // workspace keeps its allocations and frees the borrow the scratch
    // holds.
    let mut partials = core::mem::take(&mut ws.partials);
    let frozen = phase_one::<F, B, T>(batch, ws, &mut partials, context)?;
    let mut sc = scratch(ws, ncols);

    let mut out = reduce_lower_rows(batch, &mut sc, &partials, frozen, context)?;
    finish_reduction(batch, &mut sc, npiv, &mut out, context)?;
    ws.partials = partials;
    Ok(out)
}

fn reduce_lower_rows<F, B, T>(
    batch: &mut Batch,
    sc: &mut Scratch<'_>,
    partials: &[PartialRow],
    frozen: bool,
    context: &mut ReductionContext<'_, F, B, T>,
) -> Result<Reduction, F4Error>
where
    F: FieldOps<Coeff = u32>,
    B: BasisCoeffs,
    T: Trace,
{
    let mut out = Reduction::default();
    for r in 0..batch.lower().len() {
        context.clock.tick()?;
        let row = batch.lower()[r];
        if reduce_lower_row(
            batch,
            sc,
            LowerRow {
                partial: partials.get(r),
                frozen,
                index: r,
                row,
            },
            &mut out,
            context,
        )? {
            break;
        }
    }
    Ok(out)
}

struct LowerRow<'a> {
    partial: Option<&'a PartialRow>,
    frozen: bool,
    index: usize,
    row: super::matrix::RowRef,
}

fn reduce_lower_row<F, B, T>(
    batch: &mut Batch,
    sc: &mut Scratch<'_>,
    row: LowerRow<'_>,
    out: &mut Reduction,
    context: &mut ReductionContext<'_, F, B, T>,
) -> Result<bool, F4Error>
where
    F: FieldOps<Coeff = u32>,
    B: BasisCoeffs,
    T: Trace,
{
    let place = RowId::Lower(row.index as u32);
    let LoadedRow::Start(first) =
        load_row(batch, context.basis, sc, row.partial, row.frozen, row.row)
    else {
        out.zero_rows += 1;
        context.trace.start(place);
        context.trace.end(None);
        return Ok(false);
    };
    context.trace.start(place);
    eliminate(
        batch,
        context.basis,
        context.field,
        context.trace,
        sc,
        first,
        context.clock,
    )?;
    if sc.cols.is_empty() {
        out.zero_rows += 1;
        context.trace.end(None);
        return Ok(false);
    }
    let pivot = install(batch, context.field, context.trace, sc, NO_ROW)?;
    out.new_pivots += 1;
    let row = batch.pivot(pivot);
    let unit = row.len() == 1 && batch.column(row.lead()) == MonomialId::ONE;
    out.unit = unit;
    Ok(unit)
}

fn load_row<B: BasisCoeffs>(
    batch: &Batch,
    basis: &B,
    sc: &mut Scratch<'_>,
    partial: Option<&PartialRow>,
    frozen: bool,
    row: super::matrix::RowRef,
) -> LoadedRow {
    if frozen {
        let partial = partial.expect("the parallel phase wrote every lower row");
        if partial.cols.is_empty() {
            return LoadedRow::Zero;
        }
        scatter(sc.acc, &partial.cols, &partial.vals);
        return LoadedRow::Start(partial.cols[0] as usize);
    }
    let cols = batch.row_columns(&row);
    let (vals, _) = batch.row_coeffs(&row, basis);
    scatter(sc.acc, cols, vals);
    LoadedRow::Start(row.lead() as usize)
}

fn finish_reduction<F, B, T>(
    batch: &mut Batch,
    sc: &mut Scratch<'_>,
    npiv: usize,
    out: &mut Reduction,
    context: &mut ReductionContext<'_, F, B, T>,
) -> Result<(), F4Error>
where
    F: FieldOps<Coeff = u32>,
    B: BasisCoeffs,
    T: Trace,
{
    sort_new_pivots(batch, context.trace, npiv);
    if out.unit {
        return Ok(());
    }
    for i in npiv as u32..batch.pivot_count() {
        sc.pivot_at[batch.pivot(i).lead() as usize] = i;
    }
    backward(
        batch,
        context.basis,
        context.field,
        context.trace,
        sc,
        npiv as u32,
        context.clock,
    )
}

/// Sort the pivots from `npiv` on, and report the permutation.
///
/// Contract section 9.4 report 2: entry `r` of the reported array is the
/// slot the row now at slot `npiv + r` held before the sort. The recorder
/// needs it because every pivot slot the backward pass reports is a
/// post-sort slot, while the row loop bound its nodes to pre-sort slots.
/// Two new pivots never lead in one column, so sorting the slots by
/// leading column gives the same permutation the batch applies.
fn sort_new_pivots<T: Trace>(batch: &mut Batch, trace: &mut T, npiv: usize) {
    if !T::RECORDS {
        batch.sort_pivots_from(npiv);
        return;
    }
    let mut order: Vec<u32> = (npiv as u32..batch.pivot_count()).collect();
    order.sort_by_key(|&slot| batch.pivot(slot).lead());
    batch.sort_pivots_from(npiv);
    trace.pivots_sorted(npiv as u32, &order);
}

/// The columns one walk covers between two reads of the clock.
///
/// Reading the clock costs more than one column, so the walk reads it once
/// per stride of columns. Each column applies at most one pivot, so the
/// stride bounds the overshoot of a deadline: a run past it stops inside
/// the row, after at most this many more pivot applications.
const WALK_STRIDE: usize = 1024;

/// The nonzero count one batch must reach before the parallel phase pays.
///
/// The estimate is the nonzero count of the rows to reduce, which is the
/// row count times the mean row length of design section 5. Below the
/// threshold a batch runs on the calling thread and makes no rayon call.
/// The value comes from the measurement in the M2 report: noon-6 and
/// cyclic-7, whose batches sit on both sides of it.
const PARALLEL_WORK_THRESHOLD: u64 = 2_500;

/// Reduce every row against the frozen pivots, on several threads.
///
/// Returns whether the phase ran. It runs for a trace that records
/// nothing, a `threads` count other than 1, and a batch above
/// [`PARALLEL_WORK_THRESHOLD`]. Each row is a pure function of the row and
/// the frozen pivot set, and every worker owns one accumulator, so the
/// result does not depend on the thread count.
///
/// The frozen pivots are the `[A | B]` reducer rows. A pivot the kernel
/// installs later leads in a column no reducer row leads in, and it holds
/// no entry in a reducer row's column, so applying it after this phase
/// gives what the one-pass walk gives.
fn phase_one<F, B, T>(
    batch: &Batch,
    ws: &mut Workspace,
    partials: &mut Vec<PartialRow>,
    context: &mut ReductionContext<'_, F, B, T>,
) -> Result<bool, F4Error>
where
    F: FieldOps<Coeff = u32> + Sync,
    B: BasisCoeffs + Sync,
    T: Trace,
{
    if !parallel_enabled::<T>(batch, context.threads) {
        return Ok(false);
    }
    context.clock.check()?;
    let ncols = batch.ncols();
    let worker_acc_bytes = rayon::current_num_threads()
        .saturating_mul(ncols)
        .saturating_mul(size_of::<u64>());
    check_worker_memory(worker_acc_bytes, context.memory)?;
    let rows = batch.lower().len();
    prepare_partials(partials, rows)?;
    let forked = context.clock.fork();
    ws.worker_acc_bytes = worker_acc_bytes;
    let pivot_at: &[u32] = &ws.pivot_at[..ncols];
    let stopped = run_parallel_rows(
        batch,
        context.basis,
        context.field,
        partials,
        ncols,
        pivot_at,
        &forked,
    );
    if stopped {
        context.clock.check()?;
        return Err(F4Error::Timeout);
    }
    context.clock.check()?;
    Ok(true)
}

fn run_parallel_rows<F, B>(
    batch: &Batch,
    basis: &B,
    field: &F,
    partials: &mut [PartialRow],
    ncols: usize,
    pivot_at: &[u32],
    template: &Deadline,
) -> bool
where
    F: FieldOps<Coeff = u32> + Sync,
    B: BasisCoeffs + Sync,
{
    use rayon::prelude::*;

    let stopped = AtomicBool::new(false);
    let rows = batch.lower().len();
    batch
        .lower()
        .par_iter()
        .zip(partials[..rows].par_iter_mut())
        .for_each_init(
            || (vec![0u64; ncols], template.fork()),
            |(acc, clock), (row, out)| {
                if stopped.load(Ordering::Relaxed) {
                    return;
                }
                if clock.tick().is_err() {
                    stopped.store(true, Ordering::Relaxed);
                    return;
                }
                let cols = batch.row_columns(row);
                let (vals, _) = batch.row_coeffs(row, basis);
                scatter(acc, cols, vals);
                let walked = walk(
                    batch,
                    basis,
                    field,
                    row.lead() as usize,
                    &mut WalkState {
                        trace: &mut NoTrace,
                        acc,
                        pivot_at,
                        cols_out: &mut out.cols,
                        vals_out: &mut out.vals,
                        clock,
                    },
                );
                if walked.is_err() {
                    stopped.store(true, Ordering::Relaxed);
                }
            },
        );
    stopped.load(Ordering::Relaxed)
}

fn parallel_enabled<T: Trace>(batch: &Batch, threads: Option<usize>) -> bool {
    if T::RECORDS || threads == Some(1) {
        return false;
    }
    let work: u64 = batch.lower().iter().map(|row| u64::from(row.len())).sum();
    work >= PARALLEL_WORK_THRESHOLD && rayon::current_num_threads() >= 2
}

fn check_worker_memory(bytes: usize, memory: Option<usize>) -> Result<(), F4Error> {
    if memory.is_some_and(|limit| bytes > limit) {
        return Err(F4Error::MemoryLimitExceeded);
    }
    Ok(())
}

fn prepare_partials(partials: &mut Vec<PartialRow>, rows: usize) -> Result<(), F4Error> {
    if partials.len() < rows {
        reserve(partials, rows - partials.len())?;
        partials.resize_with(rows, PartialRow::default);
    }
    Ok(())
}

/// Reduce the tail of every pivot row of `batch` (design section 3.8).
///
/// The caller builds a batch whose pivot rows are the whole basis: one
/// raw row per live element with multiplier 1, then symbolic
/// preprocessing, then [`super::matrix::build`]. Symbolic closure adds the
/// reducer multiples, and every element becomes the pivot row of its own
/// leading column.
///
/// Each row keeps its leading term and reduces its tail against the rows
/// with larger leading columns, which are the rows with smaller leading
/// monomials. Those rows are already reduced when the tail meets them,
/// because the pass runs in decreasing leading column order. An element
/// never reduces itself: its own leading monomial divides no monomial of
/// its own tail.
pub(crate) fn interreduce<F, B, T>(
    batch: &mut Batch,
    basis: &B,
    field: &F,
    trace: &mut T,
    ws: &mut Workspace,
    clock: &mut Deadline,
) -> Result<(), F4Error>
where
    F: FieldOps<Coeff = u32>,
    B: BasisCoeffs,
    T: Trace,
{
    let ncols = batch.ncols();
    prepare(batch, ws, ncols)?;
    batch.debug_check_coeffs(basis);
    let mut sc = scratch(ws, ncols);
    backward(batch, basis, field, trace, &mut sc, 0, clock)
}

/// Size the accumulator and map every pivot column to its pivot row.
fn prepare(batch: &Batch, ws: &mut Workspace, ncols: usize) -> Result<(), F4Error> {
    if ws.acc.len() < ncols {
        let extra = ncols - ws.acc.len();
        reserve(&mut ws.acc, extra)?;
        ws.acc.resize(ncols, 0);
    }
    debug_assert!(
        ws.acc[..ncols].iter().all(|&lane| lane == 0),
        "a lane was left dirty"
    );
    fill(&mut ws.pivot_at, ncols, NO_ROW)?;
    for i in 0..batch.pivot_count() {
        ws.pivot_at[batch.pivot(i).lead() as usize] = i;
    }
    Ok(())
}

/// Borrow the accumulator, the pivot map, and the output row.
fn scratch(ws: &mut Workspace, ncols: usize) -> Scratch<'_> {
    Scratch {
        acc: &mut ws.acc[..ncols],
        pivot_at: &mut ws.pivot_at[..],
        cols: &mut ws.cols_out,
        vals: &mut ws.vals_out,
        shoup: &mut ws.shoup_out,
    }
}

/// Write one row into the accumulator.
fn scatter(acc: &mut [u64], cols: &[u32], vals: &[u32]) {
    debug_assert_eq!(cols.len(), vals.len(), "a row has one column per value");
    for (&col, &val) in cols.iter().zip(vals) {
        debug_assert_eq!(acc[col as usize], 0, "a lane was left dirty");
        acc[col as usize] = val as u64;
    }
}

/// Reduce the accumulator from column `first` on, and gather what is left.
///
/// The result is in `sc.cols` and `sc.vals`, ascending by column.
fn eliminate<F, B, T>(
    batch: &Batch,
    basis: &B,
    field: &F,
    trace: &mut T,
    sc: &mut Scratch<'_>,
    first: usize,
    clock: &mut Deadline,
) -> Result<(), F4Error>
where
    F: FieldOps<Coeff = u32>,
    B: BasisCoeffs,
    T: Trace,
{
    walk(
        batch,
        basis,
        field,
        first,
        &mut WalkState {
            trace,
            acc: sc.acc,
            pivot_at: sc.pivot_at,
            cols_out: sc.cols,
            vals_out: sc.vals,
            clock,
        },
    )
}

/// The elimination walk over one accumulator.
///
/// `acc` holds the row, `pivot_at` maps a column to its pivot row, and
/// the result lands in `cols` and `vals`, ascending by column. Every lane
/// the walk reads becomes physical zero, so the accumulator is clean for
/// the next row without a clearing pass.
///
/// One row is unbounded work, so the walk reads the clock once per
/// [`WALK_STRIDE`] columns. A run past its deadline stops inside the row,
/// not after it.
struct WalkState<'a, T> {
    trace: &'a mut T,
    acc: &'a mut [u64],
    pivot_at: &'a [u32],
    cols_out: &'a mut Vec<u32>,
    vals_out: &'a mut Vec<u32>,
    clock: &'a mut Deadline,
}

fn walk<F, B, T>(
    batch: &Batch,
    basis: &B,
    field: &F,
    first: usize,
    state: &mut WalkState<'_, T>,
) -> Result<(), F4Error>
where
    F: FieldOps<Coeff = u32>,
    B: BasisCoeffs,
    T: Trace,
{
    state.cols_out.clear();
    state.vals_out.clear();
    let ncols = state.acc.len();
    let between_sweeps = field.applications_between_sweeps();
    let mut applied = 0u64;
    for col in first..ncols {
        if col.is_multiple_of(WALK_STRIDE) {
            state.clock.check()?;
        }
        let val = field.reduce_acc(state.acc[col]);
        // Every lane the walk reads becomes physical zero, including a
        // lane whose value is a nonzero multiple of p. The walk covers
        // every column of the row, so no clearing pass is needed.
        state.acc[col] = 0;
        if val == 0 {
            continue;
        }
        let index = state.pivot_at[col];
        if index == NO_ROW {
            state.cols_out.push(col as u32);
            state.vals_out.push(val);
            continue;
        }
        let pivot = batch.pivot(index);
        debug_assert_eq!(pivot.lead() as usize, col, "the pivot map is wrong");
        // The lane count is enforced, not assumed: a lane takes at most
        // `between_sweeps` additions before it can pass u64.
        if applied == between_sweeps {
            field.sweep(&mut state.acc[col + 1..]);
            applied = 0;
        }
        // The pivot's leading coefficient is 1 and lane `col` is already
        // zero, so applying the whole pivot row would write the leading
        // coefficient back. Only the tail is applied.
        let factor = field.sub(0, val);
        let cols = tail(batch.row_columns(&pivot));
        let (vals, shoup) = batch.row_coeffs(&pivot, basis);
        field.axpy(state.acc, cols, tail(vals), tail(shoup), factor);
        applied += 1;
        state.trace.step(index, factor);
    }
    Ok(())
}

/// Normalize the gathered row, append it, and install it as a pivot.
fn install<F, T>(
    batch: &mut Batch,
    field: &F,
    trace: &mut T,
    sc: &mut Scratch<'_>,
    raw: u32,
) -> Result<u32, F4Error>
where
    F: FieldOps<Coeff = u32>,
    T: Trace,
{
    let lead = sc.vals[0];
    if lead != 1 {
        let scale = field.inv(lead);
        for val in sc.vals.iter_mut() {
            *val = field.reduce_acc(*val as u64 * scale as u64);
        }
        trace.scale(scale);
    }
    debug_assert_eq!(sc.vals[0], 1, "a new pivot is monic");
    // The precomputation belongs to the coefficients it was built from, so
    // it is built after the normalization and before the row reduces
    // anything.
    field.precompute(sc.vals, sc.shoup);
    let row = batch.push_owned(sc.cols, sc.vals, sc.shoup, raw)?;
    let index = batch.push_pivot(row)?;
    sc.pivot_at[row.lead() as usize] = index;
    trace.end(Some(index));
    Ok(index)
}

/// Reduce the tail of every pivot from `from` on, backward.
///
/// The pivots from `from` on are sorted by leading column. The pass takes
/// them in decreasing leading column order, so a row meets only rows that
/// are already reduced.
fn backward<F, B, T>(
    batch: &mut Batch,
    basis: &B,
    field: &F,
    trace: &mut T,
    sc: &mut Scratch<'_>,
    from: u32,
    clock: &mut Deadline,
) -> Result<(), F4Error>
where
    F: FieldOps<Coeff = u32>,
    B: BasisCoeffs,
    T: Trace,
{
    debug_assert!(
        (from + 1..batch.pivot_count()).all(|i| batch.pivot(i - 1).lead() < batch.pivot(i).lead()),
        "the backward pass needs the pivots sorted by leading column"
    );
    let mut index = batch.pivot_count();
    while index > from {
        clock.tick()?;
        index -= 1;
        let row = batch.pivot(index);
        if row.len() < 2 {
            continue;
        }
        let first = {
            let cols = batch.row_columns(&row);
            let (vals, _) = batch.row_coeffs(&row, basis);
            debug_assert_eq!(vals[0], 1, "a pivot is monic");
            scatter(sc.acc, tail(cols), tail(vals));
            cols[1] as usize
        };
        trace.start(RowId::Pivot(index));
        eliminate(batch, basis, field, trace, sc, first, clock)?;
        let unchanged = {
            let cols = batch.row_columns(&row);
            let (vals, _) = batch.row_coeffs(&row, basis);
            sc.cols[..] == cols[1..] && sc.vals[..] == vals[1..]
        };
        if unchanged {
            trace.end(Some(index));
            continue;
        }
        sc.cols.insert(0, row.lead());
        sc.vals.insert(0, 1);
        field.precompute(sc.vals, sc.shoup);
        let new = batch.push_owned(sc.cols, sc.vals, sc.shoup, row.raw())?;
        batch.set_pivot(index, new);
        trace.end(Some(index));
    }
    Ok(())
}

/// A row without its leading entry.
fn tail<T>(row: &[T]) -> &[T] {
    row.get(1..).unwrap_or(&[])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compute::f4::field::Small31;
    use crate::compute::f4::matrix::build;
    use crate::compute::f4::matrix::fixture::{Case, Row, fixture};
    use crate::compute::f4::trace::NoTrace;
    use crate::compute::f4::trace::{BatchRows, Insertion, Returned};

    const SMALL_P: u32 = 101;
    const BENCH_P: u32 = 1073741827;

    /// A trace that counts what the kernel reports.
    #[derive(Default)]
    struct Counts {
        starts: u32,
        steps: u32,
        scales: u32,
        installed: u32,
    }

    impl Trace for Counts {
        const RECORDS: bool = true;

        fn restart(&mut self) {}

        fn rows(&mut self, _rows: &BatchRows<'_>) {}

        fn pivots_sorted(&mut self, _npiv: u32, _order: &[u32]) {}

        fn inserted(&mut self, _entries: &[Insertion]) {}

        fn returned(&mut self, _entries: &[Returned]) {}

        fn start(&mut self, _row: RowId) {
            self.starts += 1;
        }

        fn step(&mut self, _pivot: u32, factor: u32) {
            assert!(factor != 0, "a step with factor zero is not recorded");
            self.steps += 1;
        }

        fn scale(&mut self, _scalar: u32) {
            self.scales += 1;
        }

        fn end(&mut self, installed: Option<u32>) {
            if installed.is_some() {
                self.installed += 1;
            }
        }
    }

    struct Rng(u64);

    impl Rng {
        fn below(&mut self, n: u64) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 33) % n
        }
    }

    /// A random matrix with `npiv` pivot columns and `nlower` rows to
    /// reduce.
    ///
    /// Every row leads in a pivot column, which is what symbolic
    /// preprocessing produces: a row is a multiple of a basis element, so
    /// its leading monomial has a divisor among the leads.
    fn random_case(rng: &mut Rng, p: u32, ncols: usize, npiv: usize, nlower: usize) -> Case {
        let mut pivot_cols: Vec<u32> = Vec::new();
        while pivot_cols.len() < npiv {
            let col = rng.below(ncols as u64) as u32;
            if !pivot_cols.contains(&col) {
                pivot_cols.push(col);
            }
        }
        pivot_cols.sort_unstable();
        let row = |lead: u32, monic: bool, rng: &mut Rng| {
            let mut cols = vec![lead];
            let mut vals = vec![if monic {
                1
            } else {
                1 + rng.below(p as u64 - 1) as u32
            }];
            for col in lead + 1..ncols as u32 {
                if rng.below(3) == 0 {
                    cols.push(col);
                    vals.push(1 + rng.below(p as u64 - 1) as u32);
                }
            }
            Row { cols, vals }
        };
        let upper = pivot_cols
            .iter()
            .map(|&lead| row(lead, true, rng))
            .collect();
        let lower = (0..nlower)
            .map(|_| {
                let lead = pivot_cols[rng.below(npiv as u64) as usize];
                row(lead, false, rng)
            })
            .collect();
        Case {
            p,
            ncols,
            pivot_cols,
            upper,
            lower,
        }
    }

    /// Reduce a case and report the new pivots and the counters.
    fn run<F: FieldOps<Coeff = u32> + Sync>(
        case: &Case,
        field: &F,
        ws: &mut Workspace,
    ) -> (Vec<Row>, Reduction) {
        let mut held = fixture(case, field);
        let (mut sym, sources) = held.split();
        let mut batch = Batch::default();
        build(&mut sym, sources, &mut batch, ws, &mut Deadline::none()).unwrap();
        let mut clock = Deadline::none();
        let mut trace = NoTrace;
        let mut context =
            ReductionContext::new(sources, field, &mut trace, &mut clock, None, Some(1));
        let out = reduce(&mut batch, ws, &mut context).unwrap();
        let pivots = batch
            .new_pivots()
            .iter()
            .map(|row| Row {
                cols: batch.row_columns(row).to_vec(),
                vals: batch.row_coeffs(row, sources).0.to_vec(),
            })
            .collect();
        (pivots, out)
    }

    fn inverse(a: u64, p: u64) -> u64 {
        let (mut result, mut base, mut exp) = (1u64, a % p, p - 2);
        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base % p;
            }
            base = base * base % p;
            exp >>= 1;
        }
        result
    }

    fn dense(ncols: usize, row: &Row) -> Vec<u64> {
        let mut out = vec![0u64; ncols];
        for (&col, &val) in row.cols.iter().zip(&row.vals) {
            out[col as usize] = val as u64;
        }
        out
    }

    fn sparse(row: &[u64]) -> Row {
        let mut out = Row {
            cols: Vec::new(),
            vals: Vec::new(),
        };
        for (col, &val) in row.iter().enumerate() {
            if val != 0 {
                out.cols.push(col as u32);
                out.vals.push(val as u32);
            }
        }
        out
    }

    /// Dense modular Gaussian elimination, the reference of design 3.7.
    fn reference(case: &Case) -> (Vec<Row>, Reduction) {
        let p = case.p as u64;
        let mut pivots: Vec<Option<Vec<u64>>> = vec![None; case.ncols];
        for row in &case.upper {
            pivots[row.cols[0] as usize] = Some(dense(case.ncols, row));
        }
        let mut out = Reduction::default();
        let mut leads: Vec<usize> = Vec::new();
        for row in &case.lower {
            let Some((lead, acc)) = reduce_reference_row(case, row, &pivots, p) else {
                out.zero_rows += 1;
                continue;
            };
            pivots[lead] = Some(acc);
            leads.push(lead);
            out.new_pivots += 1;
            if lead == case.ncols - 1 {
                out.unit = true;
                break;
            }
        }
        if !out.unit {
            backward_reference(&mut pivots, &mut leads, case.ncols, p);
        }
        leads.sort_unstable();
        let rows = leads
            .iter()
            .map(|&lead| sparse(pivots[lead].as_ref().unwrap()))
            .collect();
        (rows, out)
    }

    fn reduce_reference_row(
        case: &Case,
        row: &Row,
        pivots: &[Option<Vec<u64>>],
        p: u64,
    ) -> Option<(usize, Vec<u64>)> {
        let mut acc = dense(case.ncols, row);
        for col in 0..case.ncols {
            if acc[col] == 0 {
                continue;
            }
            if let Some(pivot) = &pivots[col] {
                subtract_reference(&mut acc, pivot, col, p);
                assert_eq!(acc[col], 0);
            }
        }
        let lead = acc.iter().position(|&value| value != 0)?;
        let scale = inverse(acc[lead], p);
        for value in &mut acc {
            *value = *value * scale % p;
        }
        Some((lead, acc))
    }

    fn backward_reference(
        pivots: &mut [Option<Vec<u64>>],
        leads: &mut [usize],
        ncols: usize,
        p: u64,
    ) {
        leads.sort_unstable();
        for &lead in leads.iter().rev() {
            let mut acc = pivots[lead].take().expect("a lead names a pivot");
            for col in lead + 1..ncols {
                if acc[col] == 0 {
                    continue;
                }
                if let Some(pivot) = &pivots[col] {
                    subtract_reference(&mut acc, pivot, col, p);
                }
            }
            pivots[lead] = Some(acc);
        }
    }

    fn subtract_reference(target: &mut [u64], pivot: &[u64], col: usize, p: u64) {
        let factor = p - target[col];
        for index in col..pivot.len() {
            target[index] = (target[index] + factor * pivot[index]) % p;
        }
    }

    fn check<F: FieldOps<Coeff = u32> + Sync>(case: &Case, field: &F) {
        let (want, want_out) = reference(case);
        let (got, got_out) = run(case, field, &mut Workspace::default());
        assert_eq!(got_out, want_out, "counters differ");
        assert_eq!(got.len(), want.len(), "pivot count differs");
        for (got, want) in got.iter().zip(&want) {
            assert_eq!(got.cols, want.cols, "columns differ");
            assert_eq!(got.vals, want.vals, "coefficients differ");
        }
    }

    #[test]
    fn random_batches_match_the_dense_reference() {
        let mut rng = Rng(0x5eed);
        let field = Small31::new(SMALL_P);
        for _ in 0..200 {
            let ncols = 4 + rng.below(20) as usize;
            let npiv = 1 + rng.below(ncols as u64 - 1) as usize;
            let nlower = 1 + rng.below(8) as usize;
            let case = random_case(&mut rng, SMALL_P, ncols, npiv, nlower);
            check(&case, &field);
        }
    }

    #[test]
    fn random_batches_match_the_reference_at_the_benchmark_modulus() {
        let mut rng = Rng(0xf00d);
        for _ in 0..100 {
            let ncols = 4 + rng.below(20) as usize;
            let npiv = 1 + rng.below(ncols as u64 - 1) as usize;
            let case = random_case(&mut rng, BENCH_P, ncols, npiv, 6);
            check(&case, &Small31::new(BENCH_P));
        }
    }

    #[test]
    fn two_runs_agree() {
        let mut rng = Rng(0xabc);
        let case = random_case(&mut rng, SMALL_P, 24, 9, 7);
        let field = Small31::new(SMALL_P);
        let first = run(&case, &field, &mut Workspace::default());
        let second = run(&case, &field, &mut Workspace::default());
        assert_eq!(first.1, second.1);
        for (a, b) in first.0.iter().zip(&second.0) {
            assert_eq!(a.cols, b.cols);
            assert_eq!(a.vals, b.vals);
        }
    }

    #[test]
    fn a_long_chain_of_applications_stays_inside_the_lane() {
        // One row that meets 24 pivots in a row, each dense to the last
        // column, at the modulus the gate measures.
        let field = Small31::new(BENCH_P);
        let ncols = 30;
        let pivot_cols: Vec<u32> = (0..24).collect();
        let upper = pivot_cols
            .iter()
            .map(|&lead| Row {
                cols: (lead..ncols as u32).collect(),
                vals: core::iter::once(1)
                    .chain((lead + 1..ncols as u32).map(|c| BENCH_P - 1 - c))
                    .collect(),
            })
            .collect();
        let lower = vec![Row {
            cols: (0..ncols as u32).collect(),
            vals: (0..ncols as u32).map(|c| BENCH_P - 2 - c).collect(),
        }];
        let case = Case {
            p: BENCH_P,
            ncols,
            pivot_cols,
            upper,
            lower,
        };
        check(&case, &field);
    }

    #[test]
    fn a_pivot_never_writes_its_leading_coefficient_back() {
        let field = Small31::new(SMALL_P);
        let case = Case {
            p: SMALL_P,
            ncols: 3,
            pivot_cols: vec![0],
            upper: vec![Row {
                cols: vec![0, 2],
                vals: vec![1, 7],
            }],
            lower: vec![Row {
                cols: vec![0, 1],
                vals: vec![4, 5],
            }],
        };
        let (rows, out) = run(&case, &field, &mut Workspace::default());
        assert_eq!(out.new_pivots, 1);
        // 4x^2 + 5x reduced by x^2 + 7 gives 5x - 28, monic x - 28/5.
        assert_eq!(rows[0].cols, [1, 2]);
        assert_eq!(rows[0].vals[0], 1);
        check(&case, &field);
    }

    #[test]
    fn a_lane_holding_a_multiple_of_p_is_cleared() {
        // The first row leaves column 2 at 101, whose residue is zero. The
        // lane must be physically zero for the second row.
        let field = Small31::new(SMALL_P);
        let case = Case {
            p: SMALL_P,
            ncols: 4,
            pivot_cols: vec![0],
            upper: vec![Row {
                cols: vec![0, 2],
                vals: vec![1, 1],
            }],
            lower: vec![
                Row {
                    cols: vec![0, 1, 2],
                    vals: vec![1, 9, 1],
                },
                Row {
                    cols: vec![0, 2],
                    vals: vec![1, 3],
                },
            ],
        };
        let (rows, out) = run(&case, &field, &mut Workspace::default());
        assert_eq!(out.new_pivots, 2);
        assert_eq!(rows[0].cols, [1]);
        assert_eq!(rows[1].cols, [2]);
        assert_eq!(rows[1].vals, [1]);
        check(&case, &field);
    }

    #[test]
    fn a_constant_pivot_stops_the_batch() {
        let field = Small31::new(SMALL_P);
        let case = Case {
            p: SMALL_P,
            ncols: 2,
            pivot_cols: vec![0],
            upper: vec![Row {
                cols: vec![0, 1],
                vals: vec![1, 3],
            }],
            lower: vec![
                Row {
                    cols: vec![0, 1],
                    vals: vec![2, 5],
                },
                Row {
                    cols: vec![0],
                    vals: vec![4],
                },
            ],
        };
        let (rows, out) = run(&case, &field, &mut Workspace::default());
        assert!(out.unit, "a nonzero constant is the unit basis");
        assert_eq!(out.new_pivots, 1, "the rows after the constant are left");
        assert_eq!(rows[0].cols, [1]);
        assert_eq!(rows[0].vals, [1]);
    }

    #[test]
    fn a_reused_workspace_leaves_no_stale_lane() {
        let mut rng = Rng(0x1234);
        let first = random_case(&mut rng, SMALL_P, 18, 6, 5);
        let second = random_case(&mut rng, SMALL_P, 12, 4, 4);
        let field = Small31::new(SMALL_P);

        let mut shared = Workspace::default();
        run(&first, &field, &mut shared);
        let (reused, reused_out) = run(&second, &field, &mut shared);
        let (fresh, fresh_out) = run(&second, &field, &mut Workspace::default());

        assert_eq!(reused_out, fresh_out);
        for (a, b) in reused.iter().zip(&fresh) {
            assert_eq!(a.cols, b.cols);
            assert_eq!(a.vals, b.vals);
        }
    }

    #[test]
    fn the_trace_sees_every_row() {
        let field = Small31::new(SMALL_P);
        let mut rng = Rng(0x77);
        let case = random_case(&mut rng, SMALL_P, 16, 5, 4);
        let mut held = fixture(&case, &field);
        let (mut sym, sources) = held.split();
        let mut batch = Batch::default();
        let mut ws = Workspace::default();
        build(
            &mut sym,
            sources,
            &mut batch,
            &mut ws,
            &mut Deadline::none(),
        )
        .unwrap();
        let mut counts = Counts::default();
        let mut clock = Deadline::none();
        let mut context =
            ReductionContext::new(sources, &field, &mut counts, &mut clock, None, Some(1));
        let out = reduce(&mut batch, &mut ws, &mut context).unwrap();
        assert!(
            counts.starts >= case.lower.len() as u32,
            "every row to reduce starts"
        );
        assert!(counts.steps > 0);
        assert!(
            counts.installed >= out.new_pivots,
            "every new pivot is reported"
        );
        assert!(counts.scales <= case.lower.len() as u32);
    }

    #[test]
    fn interreduction_reduces_every_tail() {
        // x^2 + 5x and x + 2 as pivot rows. The first row's tail holds x,
        // which the second row reduces away.
        let field = Small31::new(SMALL_P);
        let case = Case {
            p: SMALL_P,
            ncols: 3,
            pivot_cols: vec![0, 1],
            upper: vec![
                Row {
                    cols: vec![0, 1],
                    vals: vec![1, 5],
                },
                Row {
                    cols: vec![1, 2],
                    vals: vec![1, 2],
                },
            ],
            lower: Vec::new(),
        };
        let mut held = fixture(&case, &field);
        let (mut sym, sources) = held.split();
        let mut batch = Batch::default();
        let mut ws = Workspace::default();
        build(
            &mut sym,
            sources,
            &mut batch,
            &mut ws,
            &mut Deadline::none(),
        )
        .unwrap();
        interreduce(
            &mut batch,
            sources,
            &field,
            &mut NoTrace,
            &mut ws,
            &mut Deadline::none(),
        )
        .unwrap();

        let rows: Vec<Row> = batch
            .upper()
            .iter()
            .map(|row| Row {
                cols: batch.row_columns(row).to_vec(),
                vals: batch.row_coeffs(row, sources).0.to_vec(),
            })
            .collect();
        assert_eq!(rows[0].cols, [0, 2]);
        assert_eq!(rows[0].vals, [1, SMALL_P - 10]);
        assert_eq!(rows[1].cols, [1, 2]);
        assert_eq!(rows[1].vals, [1, 2]);
    }
}
