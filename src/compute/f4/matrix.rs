//! Matrix construction for one batch (design section 7).
//!
//! [`build`] turns the symbolic preprocessing output into a [`Batch`]: a
//! column numbering, one flat arena of column indices, and one [`RowRef`]
//! per row. A row that comes from a basis element names that element's
//! coefficient vector instead of copying it. Multiplying a polynomial by
//! a monomial permutes the support and leaves the coefficients alone.
//!
//! Columns descend under the table order. Design section 7 numbers the
//! pivot columns first and the free columns second, both descending, and
//! [`super::symbolic::preprocess`] hands over that numbering. [`build`]
//! merges the two blocks into one descending numbering. The reason is
//! [`super::field::FieldOps::axpy`], which needs strictly increasing
//! column indices: a row's column indices ascend in its source's term
//! order only when the numbering reverses the monomial order over every
//! column. Under the two-block numbering a row that holds a free monomial
//! above a pivot monomial has a descent, and repairing it would permute
//! the coefficients and end the reuse. Both blocks descend, so the merge
//! is one pass.

use core::cmp::Ordering;
use core::ops::Range;

use super::basis::Basis;
use super::monomial::{Lanes, MonomialId, MonomialTable};
use super::symbolic::Symbolic;
use super::{Deadline, F4Error};

/// The absent row or column index.
pub(crate) const NO_ROW: u32 = u32::MAX;

/// Where a matrix row reads its coefficients.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CoeffSrc {
    /// The coefficient vector of basis element `source`.
    Basis(u32),
    /// Slot `slot` of the batch's owned coefficients, written by the
    /// kernel for a row it computed.
    Owned(u32),
}

/// One row of a batch: a column range, a coefficient source, and a lead.
///
/// `cols` names a range of [`Batch::support`]. The columns ascend, so the
/// leading entry is the first one and `lead` is the smallest column index
/// of the row. The coefficients are as long as the range and are in the
/// same order.
#[derive(Clone, Copy, Debug)]
pub(crate) struct RowRef {
    start: u32,
    len: u32,
    coeffs: CoeffSrc,
    lead: u32,
    raw: u32,
}

impl RowRef {
    /// The range of [`Batch::support`] that holds this row's columns.
    pub(crate) fn cols(&self) -> Range<usize> {
        self.start as usize..(self.start + self.len) as usize
    }

    /// The number of nonzero entries.
    pub(crate) fn len(&self) -> u32 {
        self.len
    }

    /// The row's leading column, which is its smallest column index.
    pub(crate) fn lead(&self) -> u32 {
        self.lead
    }

    /// Where the coefficients live.
    #[cfg(test)]
    pub(crate) fn coeffs(&self) -> CoeffSrc {
        self.coeffs
    }

    /// The symbolic row this row came from, or [`NO_ROW`].
    ///
    /// The kernel keeps the value when it rewrites a row, so a reduced
    /// pivot still names the raw row it started as. A row the kernel
    /// created has no raw row.
    pub(crate) fn raw(&self) -> u32 {
        self.raw
    }
}

/// One owned coefficient vector, written by the kernel.
#[derive(Clone, Copy, Debug)]
struct Owned {
    vals: u32,
    shoup: u32,
    shoup_len: u32,
}

/// What the reduction kernel reads from the working basis.
///
/// The kernel reads two slices per basis element and never a monomial.
pub(crate) trait BasisCoeffs {
    /// The coefficients of basis element `source`, descending by
    /// monomial.
    ///
    /// The element is monic. The tail is `coeffs(source)[1..]`. Every
    /// matrix row built from this element has this length.
    fn coeffs(&self, source: u32) -> &[u32];

    /// The multiply precomputation for [`BasisCoeffs::coeffs`].
    ///
    /// [`super::field::FieldOps::precompute`] writes it. It is empty under
    /// a kernel that needs no precomputation.
    fn shoup(&self, source: u32) -> &[u64];
}

/// What [`build`] reads from the working basis.
///
/// A row is a basis element and a multiplier, so the builder needs the
/// element's monomials.
pub(crate) trait BasisRows<L: Lanes>: BasisCoeffs {
    /// The monomials of basis element `source`, strictly descending.
    ///
    /// The order decides the row's column order, so it must be the order
    /// of [`BasisCoeffs::coeffs`].
    fn monomials(&self, source: u32) -> &[MonomialId];
}

impl<L: Lanes> BasisCoeffs for Basis<L> {
    fn coeffs(&self, source: u32) -> &[u32] {
        &self.poly(source).coeffs
    }

    fn shoup(&self, source: u32) -> &[u64] {
        &self.poly(source).shoup
    }
}

impl<L: Lanes> BasisRows<L> for Basis<L> {
    fn monomials(&self, source: u32) -> &[MonomialId] {
        &self.poly(source).monos
    }
}

/// One reduction matrix, with its column numbering and its rows.
///
/// The run holds one value and passes it to [`build`] for every batch, so
/// the arenas keep their capacity. `pivots` holds the `[A | B]` reducer
/// rows first, one per pivot column and sorted by leading column, and the
/// pivots the kernel installs after them. `lower` holds the `[C | D]`
/// rows still to reduce, in symbolic row order.
#[derive(Default)]
pub(crate) struct Batch {
    columns: Vec<MonomialId>,
    npiv: u32,
    nfree: u32,
    support: Vec<u32>,
    pivots: Vec<RowRef>,
    lower: Vec<RowRef>,
    vals: Vec<u32>,
    shoup: Vec<u64>,
    owned: Vec<Owned>,
}

impl Batch {
    /// Drop every row, keeping the allocations.
    pub(crate) fn clear(&mut self) {
        self.columns.clear();
        self.npiv = 0;
        self.nfree = 0;
        self.support.clear();
        self.pivots.clear();
        self.lower.clear();
        self.vals.clear();
        self.shoup.clear();
        self.owned.clear();
    }

    /// Release every byte the batch's arenas hold above their contents.
    ///
    /// The memory retry of design section 11 calls this after
    /// [`Batch::clear`], because the estimate counts capacity.
    pub(crate) fn shrink(&mut self) {
        self.columns.shrink_to_fit();
        self.support.shrink_to_fit();
        self.pivots.shrink_to_fit();
        self.lower.shrink_to_fit();
        self.vals.shrink_to_fit();
        self.shoup.shrink_to_fit();
        self.owned.shrink_to_fit();
    }

    /// The number of columns.
    pub(crate) fn ncols(&self) -> usize {
        self.columns.len()
    }

    /// The number of pivot columns, which is the number of reducer rows.
    pub(crate) fn npiv(&self) -> u32 {
        self.npiv
    }

    /// The number of columns with no reducer row.
    #[cfg(test)]
    pub(crate) fn nfree(&self) -> u32 {
        self.nfree
    }

    /// The column monomials, descending under the table order.
    #[cfg(test)]
    pub(crate) fn columns(&self) -> &[MonomialId] {
        &self.columns
    }

    /// The monomial of column `col`.
    pub(crate) fn column(&self, col: u32) -> MonomialId {
        self.columns[col as usize]
    }

    /// The `[A | B]` reducer rows, sorted by leading column.
    pub(crate) fn upper(&self) -> &[RowRef] {
        &self.pivots[..self.npiv as usize]
    }

    /// The `[C | D]` rows still to reduce, in symbolic row order.
    pub(crate) fn lower(&self) -> &[RowRef] {
        &self.lower
    }

    /// The pivots the kernel installed, sorted by leading column.
    ///
    /// Leading columns ascend, so leading monomials descend. Design
    /// section 3.3 inserts new elements in increasing leading monomial
    /// order, so it reads this slice in reverse.
    pub(crate) fn new_pivots(&self) -> &[RowRef] {
        &self.pivots[self.npiv as usize..]
    }

    /// The columns of `row`, ascending.
    pub(crate) fn row_columns(&self, row: &RowRef) -> &[u32] {
        &self.support[row.cols()]
    }

    /// The coefficients of `row` and their precomputation.
    ///
    /// Both are as long as the row, except the precomputation, which is
    /// empty under a kernel that needs none.
    pub(crate) fn row_coeffs<'a, B: BasisCoeffs>(
        &'a self,
        row: &RowRef,
        basis: &'a B,
    ) -> (&'a [u32], &'a [u64]) {
        match row.coeffs {
            CoeffSrc::Basis(source) => (basis.coeffs(source), basis.shoup(source)),
            CoeffSrc::Owned(slot) => {
                let owned = self.owned[slot as usize];
                let vals = &self.vals[owned.vals as usize..][..row.len as usize];
                let shoup = &self.shoup[owned.shoup as usize..][..owned.shoup_len as usize];
                (vals, shoup)
            }
        }
    }

    /// The bytes the batch's allocations hold, counting capacities.
    pub(crate) fn bytes(&self) -> usize {
        capacity_bytes::<MonomialId>(self.columns.capacity())
            .saturating_add(capacity_bytes::<u32>(self.support.capacity()))
            .saturating_add(capacity_bytes::<RowRef>(self.pivots.capacity()))
            .saturating_add(capacity_bytes::<RowRef>(self.lower.capacity()))
            .saturating_add(capacity_bytes::<u32>(self.vals.capacity()))
            .saturating_add(capacity_bytes::<u64>(self.shoup.capacity()))
            .saturating_add(capacity_bytes::<Owned>(self.owned.capacity()))
    }

    /// Check the structure of the matrix, in debug builds.
    ///
    /// Three of the four conditions of design section 7 are here: one
    /// reducer row per pivot column, that row's lead equal to its pivot
    /// column, and the pivot columns in descending monomial order. The
    /// fourth condition is monic basis elements, which needs the basis
    /// and is checked in [`Batch::debug_check_coeffs`].
    pub(crate) fn debug_check<L: Lanes>(&self, table: &MonomialTable<L>) {
        if !cfg!(debug_assertions) {
            return;
        }
        self.debug_check_columns(table);
        self.debug_check_rows();
    }

    fn debug_check_columns<L: Lanes>(&self, table: &MonomialTable<L>) {
        for pair in self.columns.windows(2) {
            assert_eq!(
                table.cmp(pair[0], pair[1]),
                core::cmp::Ordering::Greater,
                "columns must descend under the table order"
            );
        }
        assert_eq!(
            self.npiv as usize + self.nfree as usize,
            self.columns.len(),
            "every column is a pivot column or a free column"
        );
        assert_eq!(
            self.npiv as usize,
            self.upper().len(),
            "a pivot column takes exactly one reducer row"
        );
        for pair in self.upper().windows(2) {
            assert!(
                pair[0].lead < pair[1].lead,
                "two reducer rows share a pivot column"
            );
        }
    }

    fn debug_check_rows(&self) {
        for row in self.pivots.iter().chain(&self.lower) {
            let cols = self.row_columns(row);
            assert!(!cols.is_empty(), "a matrix row holds at least one entry");
            assert_eq!(cols[0], row.lead, "the lead is the row's first column");
            assert!(
                cols.windows(2).all(|w| w[0] < w[1]),
                "a row's columns must ascend strictly"
            );
            assert!(
                (row.lead as usize) < self.columns.len(),
                "a column index is inside the batch"
            );
        }
    }

    /// Check the coefficients of every row, in debug builds.
    ///
    /// This is the fourth condition of design section 7: every basis
    /// element is monic, so `A` is unit upper triangular. It also checks
    /// that a row is as long as its coefficient vector.
    pub(crate) fn debug_check_coeffs<B: BasisCoeffs>(&self, basis: &B) {
        if !cfg!(debug_assertions) {
            return;
        }
        for row in self.pivots.iter().chain(&self.lower) {
            let (vals, shoup) = self.row_coeffs(row, basis);
            assert_eq!(
                vals.len(),
                row.len as usize,
                "a row is as long as its coefficient vector"
            );
            assert!(
                shoup.is_empty() || shoup.len() == vals.len(),
                "the precomputation is stale"
            );
        }
        for row in self.upper() {
            let (vals, _) = self.row_coeffs(row, basis);
            assert_eq!(vals[0], 1, "a reducer row must be monic");
        }
    }

    /// Append one row the kernel computed, and return it.
    ///
    /// `cols` and `vals` describe the row, ascending by column. `shoup` is
    /// what the field context precomputed for `vals`. The row's leading
    /// column is `cols[0]`.
    pub(super) fn push_owned(
        &mut self,
        cols: &[u32],
        vals: &[u32],
        shoup: &[u64],
        raw: u32,
    ) -> Result<RowRef, F4Error> {
        debug_assert_eq!(cols.len(), vals.len(), "a row has one column per value");
        debug_assert!(!cols.is_empty(), "an empty row is not a pivot");
        reserve(&mut self.support, cols.len())?;
        reserve(&mut self.vals, vals.len())?;
        reserve(&mut self.shoup, shoup.len())?;
        reserve(&mut self.owned, 1)?;
        let row = RowRef {
            start: self.support.len() as u32,
            len: cols.len() as u32,
            coeffs: CoeffSrc::Owned(self.owned.len() as u32),
            lead: cols[0],
            raw,
        };
        self.owned.push(Owned {
            vals: self.vals.len() as u32,
            shoup: self.shoup.len() as u32,
            shoup_len: shoup.len() as u32,
        });
        self.support.extend_from_slice(cols);
        self.vals.extend_from_slice(vals);
        self.shoup.extend_from_slice(shoup);
        Ok(row)
    }

    /// Install `row` as a new pivot and return its index in `pivots`.
    pub(super) fn push_pivot(&mut self, row: RowRef) -> Result<u32, F4Error> {
        reserve(&mut self.pivots, 1)?;
        self.pivots.push(row);
        Ok(self.pivots.len() as u32 - 1)
    }

    /// The pivot row at index `index` of the combined pivot list.
    pub(super) fn pivot(&self, index: u32) -> RowRef {
        self.pivots[index as usize]
    }

    /// Replace the pivot row at index `index`.
    pub(super) fn set_pivot(&mut self, index: u32, row: RowRef) {
        debug_assert_eq!(
            self.pivots[index as usize].lead, row.lead,
            "a pivot keeps its leading column"
        );
        self.pivots[index as usize] = row;
    }

    /// Sort the pivots from `from` on by leading column, ascending.
    pub(super) fn sort_pivots_from(&mut self, from: usize) {
        self.pivots[from..].sort_unstable_by_key(|row| row.lead);
    }

    /// The number of pivot rows: the reducer rows and the installed
    /// pivots.
    pub(super) fn pivot_count(&self) -> u32 {
        self.pivots.len() as u32
    }
}

/// Arenas and accumulators that outlive one batch.
///
/// The run holds one value. Every vector is truncated between batches and
/// never dropped, so the run keeps the allocations it already made.
#[derive(Default)]
pub(crate) struct Workspace {
    /// Maps a symbolic column to the batch column that holds it.
    col_map: Vec<u32>,
    /// One `u64` lane per column. Every lane is zero between rows.
    pub(super) acc: Vec<u64>,
    /// Maps a column to its pivot row, or [`NO_ROW`].
    pub(super) pivot_at: Vec<u32>,
    pub(super) cols_out: Vec<u32>,
    pub(super) vals_out: Vec<u32>,
    /// The precomputation for `vals_out`.
    pub(super) shoup_out: Vec<u64>,
    /// One partly reduced row per row to reduce, written by the parallel
    /// phase of the kernel.
    pub(super) partials: Vec<PartialRow>,
    /// The bytes the accumulators of the last parallel phase held, one
    /// per worker (design section 11).
    pub(super) worker_acc_bytes: usize,
}

/// One row after the parallel phase reduced it against the frozen
/// pivots.
///
/// The row is empty when it reduced to zero. Both vectors keep their
/// allocation between batches.
#[derive(Default)]
pub(super) struct PartialRow {
    pub(super) cols: Vec<u32>,
    pub(super) vals: Vec<u32>,
}

impl Workspace {
    /// Release every byte the workspace holds above its contents.
    ///
    /// The memory retry of design section 11 calls this, because the
    /// estimate counts capacity. The next batch sizes the buffers again.
    pub(crate) fn shrink(&mut self) {
        self.col_map = Vec::new();
        self.acc = Vec::new();
        self.pivot_at = Vec::new();
        self.cols_out = Vec::new();
        self.vals_out = Vec::new();
        self.shoup_out = Vec::new();
        self.partials = Vec::new();
        self.worker_acc_bytes = 0;
    }

    /// The bytes the workspace's allocations hold, counting capacities.
    pub(crate) fn bytes(&self) -> usize {
        capacity_bytes::<u32>(self.col_map.capacity())
            .saturating_add(capacity_bytes::<u64>(self.acc.capacity()))
            .saturating_add(capacity_bytes::<u32>(self.pivot_at.capacity()))
            .saturating_add(capacity_bytes::<u32>(self.cols_out.capacity()))
            .saturating_add(capacity_bytes::<u32>(self.vals_out.capacity()))
            .saturating_add(capacity_bytes::<u64>(self.shoup_out.capacity()))
            .saturating_add(self.worker_acc_bytes)
            .saturating_add(capacity_bytes::<PartialRow>(self.partials.capacity()))
            .saturating_add(
                self.partials
                    .iter()
                    .map(|row| row.cols.capacity() + row.vals.capacity())
                    .fold(0usize, usize::saturating_add)
                    .saturating_mul(size_of::<u32>()),
            )
    }
}

/// Number the columns and fill the row arena of one batch.
///
/// `sym` is the symbolic preprocessing output and `basis` holds the source
/// elements. `batch` is cleared first, so the caller keeps one value for
/// the whole run.
///
/// The pass writes one `u32` column index per nonzero and copies no
/// coefficient. It rewrites the products already interned in
/// [`Symbolic::row_monos`] into batch columns and keeps that vector as
/// the support arena.
pub(crate) fn build<L: Lanes, B: BasisRows<L>>(
    sym: &mut Symbolic<'_, L>,
    basis: &B,
    batch: &mut Batch,
    ws: &mut Workspace,
    clock: &mut Deadline,
) -> Result<(), F4Error> {
    clock.check()?;
    batch.clear();
    let ncols = sym.columns.len();
    let npiv = sym.npiv as usize;
    debug_assert!(npiv <= ncols, "the pivot block is inside the columns");
    merge_columns(sym, batch, ws)?;
    rewrite_support(sym, batch, ws, clock)?;
    build_rows(sym, basis, batch, &ws.col_map, clock)?;
    batch.pivots.sort_unstable_by_key(|row| row.lead);

    debug_assert_eq!(
        batch.pivots.len(),
        npiv,
        "a pivot column takes exactly one reducer row"
    );
    batch.debug_check(sym.table);
    Ok(())
}

fn merge_columns<L: Lanes>(
    sym: &Symbolic<'_, L>,
    batch: &mut Batch,
    ws: &mut Workspace,
) -> Result<(), F4Error> {
    let ncols = sym.columns.len();
    let npiv = sym.npiv as usize;
    reserve(&mut batch.columns, ncols)?;
    fill(&mut ws.col_map, sym.table.len(), NO_ROW)?;
    let (mut pivot, mut free) = (0, npiv);
    while pivot < npiv || free < ncols {
        let take_pivot = if pivot >= npiv {
            false
        } else if free >= ncols {
            true
        } else {
            sym.table.cmp(sym.columns[pivot], sym.columns[free]) == Ordering::Greater
        };
        let from = if take_pivot { &mut pivot } else { &mut free };
        let col = *from;
        *from += 1;
        ws.col_map[sym.columns[col].index()] = batch.columns.len() as u32;
        batch.columns.push(sym.columns[col]);
    }
    batch.npiv = sym.npiv;
    batch.nfree = sym.nfree();
    Ok(())
}

fn rewrite_support<L: Lanes>(
    sym: &mut Symbolic<'_, L>,
    batch: &mut Batch,
    ws: &Workspace,
    clock: &mut Deadline,
) -> Result<(), F4Error> {
    batch.support = core::mem::take(&mut sym.row_monos);
    for slot in batch.support.iter_mut() {
        clock.tick()?;
        *slot = ws.col_map[*slot as usize];
        debug_assert_ne!(
            *slot, NO_ROW,
            "symbolic preprocessing interns every term of every row"
        );
    }
    Ok(())
}

fn build_rows<L: Lanes, B: BasisRows<L>>(
    sym: &Symbolic<'_, L>,
    basis: &B,
    batch: &mut Batch,
    col_map: &[u32],
    clock: &mut Deadline,
) -> Result<(), F4Error> {
    let npiv = sym.npiv as usize;
    reserve(&mut batch.pivots, npiv)?;
    reserve(&mut batch.lower, sym.rows.len().saturating_sub(npiv))?;
    for r in 0..sym.rows.len() {
        clock.tick()?;
        let raw = sym.rows[r];
        let terms = basis.monomials(raw.source).len();
        let start = sym.row_start[r];
        let lead = batch.support[start as usize];
        debug_assert_eq!(
            lead,
            col_map[sym.columns[raw.lead_col as usize].index()],
            "a row leads in the column symbolic preprocessing named"
        );
        let row = RowRef {
            start,
            len: terms as u32,
            coeffs: CoeffSrc::Basis(raw.source),
            lead,
            raw: r as u32,
        };
        if sym.is_pivot_row(r as u32) {
            batch.pivots.push(row);
        } else {
            batch.lower.push(row);
        }
    }
    Ok(())
}

/// Grow `v` by `extra` elements, or report the memory limit.
pub(super) fn reserve<T>(v: &mut Vec<T>, extra: usize) -> Result<(), F4Error> {
    v.try_reserve(extra)
        .map_err(|_| F4Error::MemoryLimitExceeded)
}

/// Replace `v` with `len` copies of `value`, keeping the allocation.
pub(super) fn fill(v: &mut Vec<u32>, len: usize, value: u32) -> Result<(), F4Error> {
    v.clear();
    reserve(v, len)?;
    v.resize(len, value);
    Ok(())
}

/// The bytes `capacity` elements of `T` hold.
fn capacity_bytes<T>(capacity: usize) -> usize {
    capacity.saturating_mul(size_of::<T>())
}

#[cfg(test)]
pub(super) mod fixture {
    //! Hand built batches, for the tests of this module and the kernel.

    use super::*;
    use crate::compute::f4::field::FieldOps;
    use crate::compute::f4::monomial::Lanes8;
    use crate::compute::f4::symbolic::{NO_COLUMN, RawRow};

    /// One row, by batch column, ascending. `cols[0]` is the lead.
    #[derive(Clone, Debug)]
    pub(crate) struct Row {
        pub(crate) cols: Vec<u32>,
        pub(crate) vals: Vec<u32>,
    }

    /// One matrix: the pivot columns, the reducer rows, and the rows to
    /// reduce.
    ///
    /// Batch column `c` holds `x^(ncols - 1 - c)`, so the columns descend
    /// and the last column is the constant monomial. Every row leads in a
    /// pivot column, which is what symbolic preprocessing produces.
    #[derive(Clone, Debug)]
    pub(crate) struct Case {
        pub(crate) p: u32,
        pub(crate) ncols: usize,
        pub(crate) pivot_cols: Vec<u32>,
        pub(crate) upper: Vec<Row>,
        pub(crate) lower: Vec<Row>,
    }

    /// The basis behind a case: one element per row, multiplier 1.
    pub(crate) struct Sources {
        monos: Vec<Vec<MonomialId>>,
        coeffs: Vec<Vec<u32>>,
        shoup: Vec<Vec<u64>>,
        table: MonomialTable<Lanes8>,
    }

    impl BasisCoeffs for Sources {
        fn coeffs(&self, source: u32) -> &[u32] {
            &self.coeffs[source as usize]
        }

        fn shoup(&self, source: u32) -> &[u64] {
            &self.shoup[source as usize]
        }
    }

    impl BasisRows<Lanes8> for Sources {
        fn monomials(&self, source: u32) -> &[MonomialId] {
            &self.monos[source as usize]
        }
    }

    /// A case, ready to build.
    pub(crate) struct Fixture {
        sym_table: MonomialTable<Lanes8>,
        sources: Sources,
        rows: Vec<RawRow>,
        columns: Vec<MonomialId>,
        npiv: u32,
        col_of: Vec<u32>,
        pivot_row_of: Vec<u32>,
        row_monos: Vec<u32>,
        row_start: Vec<u32>,
    }

    impl Fixture {
        /// The symbolic output and the basis, borrowed together.
        pub(crate) fn split(&mut self) -> (Symbolic<'_, Lanes8>, &Sources) {
            let sym = Symbolic {
                table: &mut self.sym_table,
                rows: self.rows.clone(),
                columns: self.columns.clone(),
                npiv: self.npiv,
                col_of: self.col_of.clone(),
                pivot_row_of: self.pivot_row_of.clone(),
                row_monos: self.row_monos.clone(),
                row_start: self.row_start.clone(),
            };
            (sym, &self.sources)
        }
    }

    /// Turn a case into a fixture, with the columns in the two blocks
    /// symbolic preprocessing hands over.
    pub(crate) fn fixture<F: FieldOps<Coeff = u32>>(case: &Case, field: &F) -> Fixture {
        let ncols = case.ncols;
        let mono = |c: u32| [(ncols - 1 - c as usize) as u32];
        let mut basis_table = MonomialTable::<Lanes8>::new(1);
        let mut sym_table = MonomialTable::<Lanes8>::new(1);

        let mut order = case.pivot_cols.clone();
        order.extend((0..ncols as u32).filter(|c| !case.pivot_cols.contains(c)));
        let columns: Vec<MonomialId> = order
            .iter()
            .map(|&c| sym_table.intern(&mono(c)).unwrap())
            .collect();
        let mut col_of = vec![NO_COLUMN; sym_table.len()];
        for (col, id) in columns.iter().enumerate() {
            col_of[id.index()] = col as u32;
        }
        let column_of = |c: u32| order.iter().position(|&o| o == c).unwrap() as u32;

        let mut sources = Sources {
            monos: Vec::new(),
            coeffs: Vec::new(),
            shoup: Vec::new(),
            table: MonomialTable::<Lanes8>::new(1),
        };
        let mut rows = Vec::new();
        let mut row_monos = Vec::new();
        let mut row_start = Vec::new();
        let mut pivot_row_of = vec![u32::MAX; case.pivot_cols.len()];
        for (r, row) in case.upper.iter().chain(&case.lower).enumerate() {
            let lead_col = column_of(row.cols[0]);
            rows.push(RawRow {
                source: r as u32,
                mult: MonomialId::ONE,
                lead_col,
            });
            if r < case.upper.len() {
                pivot_row_of[lead_col as usize] = r as u32;
            }
            sources.monos.push(
                row.cols
                    .iter()
                    .map(|&c| basis_table.intern(&mono(c)).unwrap())
                    .collect(),
            );
            // Every multiplier is 1, so a row's products are its own
            // column monomials in the symbolic table.
            row_start.push(row_monos.len() as u32);
            row_monos.extend(
                row.cols
                    .iter()
                    .map(|&c| sym_table.intern(&mono(c)).unwrap().index() as u32),
            );
            let mut shoup = Vec::new();
            field.precompute(&row.vals, &mut shoup);
            sources.coeffs.push(row.vals.clone());
            sources.shoup.push(shoup);
        }
        sources.table = basis_table;

        Fixture {
            sym_table,
            sources,
            rows,
            columns,
            npiv: case.pivot_cols.len() as u32,
            col_of,
            pivot_row_of,
            row_monos,
            row_start,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::fixture::{Case, Row, fixture};
    use super::*;
    use crate::compute::f4::field::Small31;

    const P: u32 = 101;

    /// A case whose rows are given by batch column.
    fn case(ncols: usize, pivot_cols: &[u32], upper: &[&[u32]], lower: &[&[u32]]) -> Case {
        let row = |cols: &&[u32], monic: bool| Row {
            cols: cols.to_vec(),
            vals: cols
                .iter()
                .enumerate()
                .map(|(k, &c)| if k == 0 && monic { 1 } else { 2 + c % 7 })
                .collect(),
        };
        Case {
            p: P,
            ncols,
            pivot_cols: pivot_cols.to_vec(),
            upper: upper.iter().map(|cols| row(cols, true)).collect(),
            lower: lower.iter().map(|cols| row(cols, false)).collect(),
        }
    }

    #[test]
    fn columns_descend_and_row_supports_ascend() {
        let case = case(5, &[0, 2], &[&[0, 1, 4], &[2, 3]], &[]);
        let mut fixture = fixture(&case, &Small31::new(P));
        let (mut sym, sources) = fixture.split();
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

        assert_eq!(batch.ncols(), 5);
        assert_eq!(batch.npiv(), 2);
        assert_eq!(batch.nfree(), 3);
        for pair in batch.columns().windows(2) {
            assert_eq!(sym.table.cmp(pair[0], pair[1]), Ordering::Greater);
        }
        assert_eq!(batch.row_columns(&batch.upper()[0]), [0, 1, 4]);
        assert_eq!(batch.row_columns(&batch.upper()[1]), [2, 3]);
        assert!(batch.lower().is_empty());
    }

    #[test]
    fn the_reducer_row_of_a_column_is_the_one_symbolic_named() {
        // Three rows lead in column 1. Symbolic preprocessing names the
        // first as the reducer row, and the other two are rows to reduce.
        let case = case(4, &[1], &[&[1, 2]], &[&[1, 3], &[1, 2, 3]]);
        let mut fixture = fixture(&case, &Small31::new(P));
        let (mut sym, sources) = fixture.split();
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

        assert_eq!(batch.npiv(), 1);
        assert_eq!(batch.upper()[0].coeffs(), CoeffSrc::Basis(0));
        assert_eq!(batch.upper()[0].raw(), 0);
        assert_eq!(batch.lower().len(), 2);
        assert_eq!(batch.lower()[0].coeffs(), CoeffSrc::Basis(1));
        assert_eq!(batch.lower()[1].coeffs(), CoeffSrc::Basis(2));
    }

    #[test]
    fn a_second_build_keeps_the_arenas() {
        let case = case(6, &[0, 1], &[&[0, 2, 5], &[1, 3, 4]], &[]);
        let mut fixture = fixture(&case, &Small31::new(P));
        let mut batch = Batch::default();
        let mut ws = Workspace::default();
        {
            let (mut sym, sources) = fixture.split();
            build(
                &mut sym,
                sources,
                &mut batch,
                &mut ws,
                &mut Deadline::none(),
            )
            .unwrap();
        }
        let support = batch.support.capacity();
        let bytes = batch.bytes();
        assert!(support >= 6 && bytes > 0);

        let (mut sym, sources) = fixture.split();
        build(
            &mut sym,
            sources,
            &mut batch,
            &mut ws,
            &mut Deadline::none(),
        )
        .unwrap();
        assert_eq!(batch.support.capacity(), support);
        assert_eq!(batch.bytes(), bytes);
        assert_eq!(batch.support.len(), 6);
        assert!(ws.bytes() > 0);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "a row's columns must ascend strictly")]
    fn a_row_that_does_not_descend_is_caught() {
        let mut case = case(3, &[2], &[&[2, 0]], &[]);
        // Undo the ascending order the case builder keeps.
        case.upper[0].cols = vec![2, 0];
        let mut fixture = fixture(&case, &Small31::new(P));
        let (mut sym, sources) = fixture.split();
        build(
            &mut sym,
            sources,
            &mut Batch::default(),
            &mut Workspace::default(),
            &mut Deadline::none(),
        )
        .unwrap();
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "two reducer rows share a pivot column")]
    fn two_reducer_rows_in_one_column_are_caught() {
        let case = case(4, &[0, 1], &[&[0, 2], &[1, 3]], &[]);
        let mut fixture = fixture(&case, &Small31::new(P));
        let (mut sym, sources) = fixture.split();
        let mut batch = Batch::default();
        build(
            &mut sym,
            sources,
            &mut batch,
            &mut Workspace::default(),
            &mut Deadline::none(),
        )
        .unwrap();
        batch.pivots[1].lead = batch.pivots[0].lead;
        batch.debug_check(sym.table);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "a reducer row must be monic")]
    fn a_reducer_row_that_is_not_monic_is_caught() {
        let mut case = case(3, &[0], &[&[0, 1]], &[]);
        case.upper[0].vals[0] = 2;
        let mut fixture = fixture(&case, &Small31::new(P));
        let (mut sym, sources) = fixture.split();
        let mut batch = Batch::default();
        build(
            &mut sym,
            sources,
            &mut batch,
            &mut Workspace::default(),
            &mut Deadline::none(),
        )
        .unwrap();
        batch.debug_check_coeffs(sources);
    }

    #[test]
    fn an_expired_deadline_stops_build() {
        let case = case(1, &[0], &[&[0]], &[]);
        let mut fixture = fixture(&case, &Small31::new(P));
        let (mut sym, sources) = fixture.split();
        let past = std::time::Instant::now()
            .checked_sub(std::time::Duration::from_secs(1))
            .expect("the clock is at least one second past its epoch");
        let mut clock = super::Deadline::new(Some(past));
        let result = build(
            &mut sym,
            sources,
            &mut Batch::default(),
            &mut Workspace::default(),
            &mut clock,
        );

        assert_eq!(result, Err(F4Error::Timeout));
    }
}
