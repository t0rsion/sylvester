//! The recording hook of the F4 engine (design sections 3.2 and 6.3).
//!
//! The run is generic over [`Trace`], so the raw path carries no trace
//! branch and no trace allocation. [`NoTrace`] is the raw path's
//! implementation. The recorder that lowers the reports into certificate
//! nodes is `crate::cert::v2`, because the shape of a node belongs to the
//! certificate contract and not to the engine.
//!
//! Four of the methods carry the reports of `docs/certificate-v2.md`
//! section 9.4: [`Trace::rows`], [`Trace::pivots_sorted`],
//! [`Trace::inserted`], and [`Trace::returned`]. The other four report one
//! row as the kernel reduces it. The engine builds a report only when
//! [`Trace::RECORDS`] holds, so the raw path builds none.

/// A row of one batch, for the trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RowId {
    /// Row `index` of [`super::matrix::Batch::lower`].
    Lower(u32),
    /// Row `index` of the batch's pivot list.
    Pivot(u32),
}

/// One row of a batch: where it sits and what it multiplies.
///
/// The multiplier is in [`BatchRows::mults`], at the row's own position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct BatchRow {
    /// The row's place in the batch.
    pub(crate) place: RowId,
    /// The global basis index the row multiplies.
    pub(crate) source: u32,
}

/// The rows of one batch (contract section 9.4, report 1).
///
/// The pivot rows come first, in increasing pivot slot, then the rows to
/// reduce, in increasing lower slot. Every slot of both lists appears
/// exactly once. `mults` holds `nvars` exponents per row, in the order of
/// `rows`.
pub(crate) struct BatchRows<'a> {
    /// One entry per row of the batch.
    pub(crate) rows: &'a [BatchRow],
    /// The multipliers, `nvars` exponents per row.
    pub(crate) mults: &'a [u32],
    /// The number of variables, which is the stride of `mults`.
    pub(crate) nvars: usize,
}

impl BatchRows<'_> {
    /// The multiplier of the row at `index` of [`BatchRows::rows`].
    pub(crate) fn mult(&self, index: usize) -> &[u32] {
        &self.mults[index * self.nvars..(index + 1) * self.nvars]
    }
}

/// One element a batch added to the basis (contract section 9.4, report
/// 3).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Insertion {
    /// The post-sort pivot slot the element came from.
    pub(crate) pivot: u32,
    /// The global basis index the element received, or `None` when the
    /// run stopped before it.
    pub(crate) basis: Option<u32>,
}

/// One element of the basis a run returns (contract section 9.4, report
/// 4).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Returned {
    /// The pivot at this slot of the last batch.
    Pivot(u32),
    /// The generator at this input index, which seeding found to be a
    /// nonzero constant.
    Input(u32),
}

/// What the engine reports while it computes.
///
/// The kernel calls [`Trace::start`], [`Trace::step`], [`Trace::scale`],
/// and [`Trace::end`] in the order it applies the operations, so a
/// recorder builds one node per row. The run loop calls the other four
/// methods at the points contract section 9.4 fixes.
pub(crate) trait Trace {
    /// Whether the recorder keeps what it is told.
    ///
    /// The parallel phase of [`super::kernel::reduce`] reduces rows
    /// against the frozen pivot set on several threads, and it reports no
    /// operation, so a run that records stays sequential. The engine also
    /// skips building a report when this is false.
    const RECORDS: bool;

    /// The run starts again at a wider lane width (design section 2.3).
    ///
    /// The recorder drops every node and every map of the failed attempt
    /// (contract section 9.2).
    fn restart(&mut self);

    /// The rows of the batch that is about to reduce.
    fn rows(&mut self, rows: &BatchRows<'_>);

    /// The pivots from `npiv` on were sorted by leading column.
    ///
    /// `order[r]` is the slot the row now at slot `npiv + r` held before
    /// the sort.
    fn pivots_sorted(&mut self, npiv: u32, order: &[u32]);

    /// The new pivots of the batch went into the basis, in this order.
    fn inserted(&mut self, entries: &[Insertion]);

    /// The run returns this basis.
    fn returned(&mut self, entries: &[Returned]);

    /// A row starts reducing.
    fn start(&mut self, row: RowId);

    /// `factor` times the tail of pivot `pivot` was added to the row.
    fn step(&mut self, pivot: u32, factor: u32);

    /// The row was multiplied by `scalar`.
    fn scale(&mut self, scalar: u32);

    /// The row ended, as pivot `installed` or as zero.
    fn end(&mut self, installed: Option<u32>);

    /// The bytes this trace holds.
    ///
    /// The run's memory estimate counts them, so one memory limit covers
    /// the engine and the recorder together. A trace that records nothing
    /// holds nothing.
    fn held_bytes(&self) -> usize {
        0
    }
}

/// The trace of a run that records nothing.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct NoTrace;

impl Trace for NoTrace {
    const RECORDS: bool = false;

    fn restart(&mut self) {}
    fn rows(&mut self, _rows: &BatchRows<'_>) {}
    fn pivots_sorted(&mut self, _npiv: u32, _order: &[u32]) {}
    fn inserted(&mut self, _entries: &[Insertion]) {}
    fn returned(&mut self, _entries: &[Returned]) {}
    fn start(&mut self, _row: RowId) {}
    fn step(&mut self, _pivot: u32, _factor: u32) {}
    fn scale(&mut self, _scalar: u32) {}
    fn end(&mut self, _installed: Option<u32>) {}
}
