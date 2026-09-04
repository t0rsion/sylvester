//! Symbolic preprocessing for one F4 batch (design section 3.5).
//!
//! Preprocessing turns the selected pairs into the rows and the columns of
//! one matrix. A row is a basis element and a multiplier, never a copy of
//! the element's coefficients: multiplying by a monomial permutes the
//! support and leaves the coefficients alone (design section 3.6).
//!
//! Every monomial here lives in the batch-local table, which the caller
//! clears once the batch's new basis elements are converted.

use super::basis::{Basis, LeadView, reserve};
use super::monomial::{Lanes, MonomialId, MonomialTable, intern_from, quotient_into};
use super::pairs::Pair;
use super::{Deadline, F4Error};

/// The column index of a monomial that is not a column of the batch.
pub(crate) const NO_COLUMN: u32 = u32::MAX;

/// The row index of a pivot column that has no reducer row.
///
/// Every pivot column of a finished [`Symbolic`] has one, so this value
/// appears only while preprocessing runs.
pub(crate) const NO_ROW: u32 = u32::MAX;

/// One row of the batch, as a basis element and a multiplier.
///
/// The row holds the terms `mult * t` for every term `t` of basis element
/// `source`. The multiplier lives in the batch-local table of
/// [`Symbolic`]. The source's monomials live in the basis table, so the
/// column of a term comes from
/// [`MonomialTable::mul_external`](super::monomial::MonomialTable::mul_external).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RawRow {
    /// The basis element the row multiplies.
    pub(crate) source: u32,
    /// The multiplier, in the batch-local table.
    pub(crate) mult: MonomialId,
    /// The column of `mult * lm(source)`, which is the row's leading
    /// column. It is always a pivot column.
    pub(crate) lead_col: u32,
}

/// The rows and the columns of one batch.
///
/// Preprocessing sorts the columns, so the matrix builder takes the order
/// as given. The order is
///
/// ```text
/// [pivot columns, descending] [non-pivot columns, descending]
/// ```
///
/// with `npiv` the boundary. Each block descends under the table order.
/// The whole vector does not. `columns[c]` is the monomial of column `c`,
/// and `col_of[m.index()]` is the column of monomial `m`, or [`NO_COLUMN`].
///
/// Every row's `lead_col` is a pivot column, and every pivot column has
/// exactly one reducer row, named by `pivot_row_of`. Rows that are not
/// named there are the rows still to reduce.
pub(crate) struct Symbolic<'t, L: Lanes> {
    /// The batch-local monomial table, holding every column and every
    /// multiplier.
    pub(crate) table: &'t mut MonomialTable<L>,
    /// The rows, reducer rows and rows to reduce in one vector.
    pub(crate) rows: Vec<RawRow>,
    /// The column monomials, pivot block then non-pivot block.
    pub(crate) columns: Vec<MonomialId>,
    /// The number of pivot columns, which is the block boundary.
    pub(crate) npiv: u32,
    /// The column of each monomial of the table, or [`NO_COLUMN`].
    pub(crate) col_of: Vec<u32>,
    /// The reducer row of each pivot column, one entry per pivot column.
    pub(crate) pivot_row_of: Vec<u32>,
    /// The monomials of every row, one segment per row, in row order.
    ///
    /// A row's segment holds `mult * t` for every term monomial `t` of
    /// the row's source, in the source's descending term order, as
    /// [`MonomialId::index`] values. The walk computes these products to
    /// close the row set. [`super::matrix::build`] takes the vector,
    /// rewrites each entry into a batch column, and keeps it as the
    /// batch's support arena, so no product is computed or copied twice.
    pub(crate) row_monos: Vec<u32>,
    /// The first index of each row's segment of `row_monos`.
    pub(crate) row_start: Vec<u32>,
}

impl<L: Lanes> Symbolic<'_, L> {
    /// The number of non-pivot columns.
    pub(crate) fn nfree(&self) -> u32 {
        self.columns.len() as u32 - self.npiv
    }

    /// The column of `m`, or `None` when `m` is not a column.
    pub(crate) fn column_of(&self, m: MonomialId) -> Option<u32> {
        match self.col_of.get(m.index()) {
            Some(&NO_COLUMN) | None => None,
            Some(&col) => Some(col),
        }
    }

    /// Report whether column `col` is a pivot column.
    pub(crate) fn is_pivot_column(&self, col: u32) -> bool {
        col < self.npiv
    }

    /// Report whether row `row` is the reducer row of its leading column.
    ///
    /// Two rows can share a leading column, because two basis elements can
    /// share a leading monomial. Only one of them reduces the column.
    pub(crate) fn is_pivot_row(&self, row: u32) -> bool {
        let col = self.rows[row as usize].lead_col;
        self.is_pivot_column(col) && self.pivot_row_of[col as usize] == row
    }

    /// The bytes the result holds, counting capacities and the table.
    pub(crate) fn memory_bytes(&self) -> usize {
        self.table
            .memory_bytes()
            .saturating_add(self.rows.capacity().saturating_mul(size_of::<RawRow>()))
            .saturating_add(
                self.columns
                    .capacity()
                    .saturating_mul(size_of::<MonomialId>()),
            )
            .saturating_add(self.col_of.capacity().saturating_mul(size_of::<u32>()))
            .saturating_add(
                self.pivot_row_of
                    .capacity()
                    .saturating_mul(size_of::<u32>()),
            )
            .saturating_add(self.row_monos.capacity().saturating_mul(size_of::<u32>()))
            .saturating_add(self.row_start.capacity().saturating_mul(size_of::<u32>()))
    }
}

/// Which divisor becomes the reducer of a column.
///
/// The default is [`Reducer::First`]. [`Reducer::Shortest`] takes a
/// reducer of higher degree more often, and its multiple pulls more
/// monomials into the symbolic closure, so the batch carries more rows and
/// more nonzeros. It was slower on every cell M2 chunk F2 measured, whose
/// counts are in that record.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum Reducer {
    /// The first divisor in scan order.
    #[default]
    First,
    /// The divisor with the fewest terms, smallest basis index first.
    Shortest,
}

/// Which row of an lcm group becomes the reducer row of that column.
///
/// The default is [`UpperRow::SmallestSource`]. The rule changes no
/// counter, only which row of an lcm group reduces its column, so it acts
/// through fill alone. Neither rule won the cells M2 chunk F2 measured, so
/// the default stays where M1 put it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum UpperRow {
    /// The row whose source has the smallest basis index.
    #[default]
    SmallestSource,
    /// The row whose source has the fewest terms, smallest basis index
    /// first.
    ShortestSource,
}

/// The strategy choices of one run.
///
/// A third choice is not here: the order the kernel takes the rows of one
/// batch in. Leading column then fewest terms, and leading column then
/// density, both cost time against the batch order and neither changed a
/// counter, on the cells M2 chunk F2 measured. The kernel keeps the batch
/// order and carries no order option.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Strategy {
    /// The reducer choice of step 4.
    pub(crate) reducer: Reducer,
    /// The reducer row choice of step 2.
    pub(crate) upper_row: UpperRow,
}

/// The mark of a monomial the walk has not classified yet.
const UNKNOWN: u8 = 0;
/// The mark of a monomial that has a reducer row.
const PIVOT: u8 = 1;
/// The mark of a monomial that has no reducer.
const FREE: u8 = 2;

/// Build the rows and the columns of one batch (design section 3.5).
///
/// `selected` is one batch of pairs, and `table` is the batch-local table,
/// which must be empty. Every row is a live or retired basis element and a
/// multiplier. Every column monomial that a live basis element's lead
/// divides gets a reducer row, so the pivot block of the result is a
/// complete reducer set for the batch.
pub(crate) fn preprocess<'t, L: Lanes>(
    selected: &[Pair],
    basis: &Basis<L>,
    table: &'t mut MonomialTable<L>,
    strategy: &Strategy,
    clock: &mut Deadline,
) -> Result<Symbolic<'t, L>, F4Error> {
    debug_assert!(!selected.is_empty(), "a batch holds at least one pair");
    debug_assert_eq!(table.len(), 1, "the batch table starts empty");

    let pairs = sorted_pairs(selected, basis)?;
    let view = basis.lead_index().project(basis.table(), table)?;
    let mut walk = Walk::default();
    add_pair_rows(&pairs, basis, table, strategy, clock, &mut walk)?;
    walk.close(table, basis, &view, strategy, clock)?;
    walk.finish(table)
}

fn sorted_pairs<L: Lanes>(selected: &[Pair], basis: &Basis<L>) -> Result<Vec<Pair>, F4Error> {
    let mut pairs: Vec<Pair> = Vec::new();
    reserve(&mut pairs, selected.len())?;
    pairs.extend_from_slice(selected);
    pairs.sort_unstable_by(|a, b| {
        basis
            .table()
            .cmp(a.lcm, b.lcm)
            .then(a.i.cmp(&b.i))
            .then(a.j.cmp(&b.j))
    });
    Ok(pairs)
}

fn add_pair_rows<L: Lanes>(
    pairs: &[Pair],
    basis: &Basis<L>,
    table: &mut MonomialTable<L>,
    strategy: &Strategy,
    clock: &mut Deadline,
    walk: &mut Walk,
) -> Result<(), F4Error> {
    let mut sources: Vec<u32> = Vec::new();
    let mut group = PairGroup {
        basis,
        table,
        strategy,
        walk,
        sources: &mut sources,
    };
    let mut start = 0;
    while start < pairs.len() {
        clock.tick()?;
        let lcm = pairs[start].lcm;
        let end = pair_group_end(pairs, start, lcm);
        add_pair_group(&pairs[start..end], lcm, &mut group)?;
        start = end;
    }
    Ok(())
}

fn pair_group_end(pairs: &[Pair], start: usize, lcm: MonomialId) -> usize {
    let mut end = start;
    while end < pairs.len() && pairs[end].lcm == lcm {
        end += 1;
    }
    end
}

struct PairGroup<'a, L: Lanes> {
    basis: &'a Basis<L>,
    table: &'a mut MonomialTable<L>,
    strategy: &'a Strategy,
    walk: &'a mut Walk,
    sources: &'a mut Vec<u32>,
}

fn add_pair_group<L: Lanes>(
    pairs: &[Pair],
    lcm: MonomialId,
    group: &mut PairGroup<'_, L>,
) -> Result<(), F4Error> {
    pair_sources(pairs, group.sources)?;
    let column = intern_from(group.table, group.basis.table(), lcm)?;
    group.walk.mark(column, PIVOT);
    let first = group.walk.rows.len() as u32;
    for &source in group.sources.iter() {
        let mult = quotient_into(
            group.table,
            group.basis.table(),
            lcm,
            group.basis.table(),
            group.basis.lead(source),
        )?
        .expect("the lead of an endpoint divides the pair's lcm");
        group.walk.push_row(source, mult, column)?;
    }
    let pivot = first
        + upper_row(
            group.basis,
            &group.walk.rows[first as usize..],
            group.strategy.upper_row,
        );
    group.walk.set_pivot_row(column, pivot);
    reserve(&mut group.walk.pivot_cols, 1)?;
    group.walk.pivot_cols.push(column);
    Ok(())
}

fn pair_sources(pairs: &[Pair], sources: &mut Vec<u32>) -> Result<(), F4Error> {
    sources.clear();
    for pair in pairs {
        reserve(sources, 2)?;
        sources.push(pair.i);
        sources.push(pair.j);
    }
    sources.sort_unstable();
    sources.dedup();
    Ok(())
}

/// Build the rows and the columns of the final interreduction.
///
/// Every live basis element becomes one row with multiplier 1, so the
/// element is the reducer row of its own leading column. Symbolic closure
/// then adds one row per reducer multiple, exactly as it does for a pair
/// batch. The caller passes the result to [`super::matrix::build`] and
/// then to [`super::kernel::interreduce`] (design section 3.8).
pub(crate) fn preprocess_basis<'t, L: Lanes>(
    basis: &Basis<L>,
    table: &'t mut MonomialTable<L>,
    strategy: &Strategy,
    clock: &mut Deadline,
) -> Result<Symbolic<'t, L>, F4Error> {
    debug_assert!(!basis.live().is_empty(), "interreduction needs an element");
    debug_assert_eq!(table.len(), 1, "the batch table starts empty");

    let view = basis.lead_index().project(basis.table(), table)?;
    let mut walk = Walk::default();
    for &element in basis.live() {
        let column = intern_from(table, basis.table(), basis.lead(element))?;
        // Two live elements never share a leading monomial, so each
        // element is the only reducer row of its own column.
        debug_assert_eq!(walk.mark_of(column), UNKNOWN, "a lead appears once");
        walk.mark(column, PIVOT);
        let row = walk.rows.len() as u32;
        walk.push_row(element, MonomialId::ONE, column)?;
        walk.set_pivot_row(column, row);
        reserve(&mut walk.pivot_cols, 1)?;
        walk.pivot_cols.push(column);
    }
    walk.close(table, basis, &view, strategy, clock)?;
    walk.finish(table)
}

/// Which row of an lcm group reduces that column.
fn upper_row<L: Lanes>(basis: &Basis<L>, group: &[RawRow], rule: UpperRow) -> u32 {
    match rule {
        // The group's sources are sorted, so the first row has the
        // smallest source index.
        UpperRow::SmallestSource => 0,
        UpperRow::ShortestSource => group
            .iter()
            .enumerate()
            .min_by_key(|(index, row)| (basis.poly(row.source).len(), *index))
            .map(|(index, _)| index as u32)
            .unwrap_or(0),
    }
}

/// The state of one symbolic preprocessing run.
#[derive(Default)]
struct Walk {
    rows: Vec<RawRow>,
    /// The products of every row, one segment per row, in row order.
    row_monos: Vec<u32>,
    /// The first index of each row's segment of `row_monos`.
    row_start: Vec<u32>,
    /// The leading monomial of each row, parallel to `rows`. It becomes
    /// the row's `lead_col` once the columns are numbered.
    row_lead: Vec<MonomialId>,
    /// The mark of each monomial of the table.
    marks: Vec<u8>,
    /// The reducer row of each monomial that is a pivot column.
    pivot_row: Vec<u32>,
    pivot_cols: Vec<MonomialId>,
    free_cols: Vec<MonomialId>,
    stack: Vec<MonomialId>,
}

impl Walk {
    fn mark_of(&self, m: MonomialId) -> u8 {
        self.marks.get(m.index()).copied().unwrap_or(UNKNOWN)
    }

    fn mark(&mut self, m: MonomialId, mark: u8) {
        if self.marks.len() <= m.index() {
            self.marks.resize(m.index() + 1, UNKNOWN);
        }
        self.marks[m.index()] = mark;
    }

    fn set_pivot_row(&mut self, m: MonomialId, row: u32) {
        if self.pivot_row.len() <= m.index() {
            self.pivot_row.resize(m.index() + 1, NO_ROW);
        }
        self.pivot_row[m.index()] = row;
    }

    fn push_row(&mut self, source: u32, mult: MonomialId, lead: MonomialId) -> Result<(), F4Error> {
        reserve(&mut self.rows, 1)?;
        reserve(&mut self.row_lead, 1)?;
        self.rows.push(RawRow {
            source,
            mult,
            lead_col: NO_COLUMN,
        });
        self.row_lead.push(lead);
        Ok(())
    }

    /// Close the row set under reduction (design section 3.5, steps 3 to
    /// 5).
    ///
    /// Every monomial of every row already in the walk goes on the stack.
    /// A monomial with a live divisor becomes a pivot column and adds the
    /// reducer's row, whose monomials go on the stack in turn. A monomial
    /// with no divisor becomes a non-pivot column.
    fn close<L: Lanes>(
        &mut self,
        table: &mut MonomialTable<L>,
        basis: &Basis<L>,
        view: &LeadView,
        strategy: &Strategy,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        self.seed_stack(table, basis, clock)?;

        while let Some(m) = self.stack.pop() {
            clock.tick()?;
            self.classify(table, basis, view, strategy.reducer, m)?;
        }
        Ok(())
    }

    fn seed_stack<L: Lanes>(
        &mut self,
        table: &mut MonomialTable<L>,
        basis: &Basis<L>,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        for index in 0..self.rows.len() {
            clock.tick()?;
            self.push_monomials(table, basis, index)?;
        }
        Ok(())
    }

    fn classify<L: Lanes>(
        &mut self,
        table: &mut MonomialTable<L>,
        basis: &Basis<L>,
        view: &LeadView,
        reducer: Reducer,
        mono: MonomialId,
    ) -> Result<(), F4Error> {
        if self.mark_of(mono) != UNKNOWN {
            return Ok(());
        }
        let Some(divisor) = view.find_divisor(table, mono, reducer)? else {
            return self.add_free(mono);
        };
        self.add_pivot(table, basis, mono, divisor.basis, divisor.mult)
    }

    fn add_pivot<L: Lanes>(
        &mut self,
        table: &mut MonomialTable<L>,
        basis: &Basis<L>,
        mono: MonomialId,
        source: u32,
        mult: MonomialId,
    ) -> Result<(), F4Error> {
        self.mark(mono, PIVOT);
        let row = self.rows.len() as u32;
        self.push_row(source, mult, mono)?;
        self.set_pivot_row(mono, row);
        reserve(&mut self.pivot_cols, 1)?;
        self.pivot_cols.push(mono);
        self.push_monomials(table, basis, row as usize)
    }

    fn add_free(&mut self, mono: MonomialId) -> Result<(), F4Error> {
        self.mark(mono, FREE);
        reserve(&mut self.free_cols, 1)?;
        self.free_cols.push(mono);
        Ok(())
    }

    /// Push every monomial of one row, in the row's descending term order.
    fn push_monomials<L: Lanes>(
        &mut self,
        table: &mut MonomialTable<L>,
        basis: &Basis<L>,
        index: usize,
    ) -> Result<(), F4Error> {
        let row = self.rows[index];
        // The matrix builder reads a row's segment by row index, so the
        // segments must arrive in row order and one per row.
        debug_assert_eq!(self.row_start.len(), index, "a row is closed once");
        let poly = basis.poly(row.source);
        reserve(&mut self.stack, poly.len())?;
        reserve(&mut self.row_monos, poly.len())?;
        reserve(&mut self.row_start, 1)?;
        self.row_start.push(self.row_monos.len() as u32);
        for &mono in &poly.monos {
            let product = table.mul_external(row.mult, basis.table(), mono)?;
            self.row_monos.push(product.index() as u32);
            self.stack.push(product);
        }
        Ok(())
    }

    /// Number the columns and finish the result.
    fn finish<L: Lanes>(
        mut self,
        table: &mut MonomialTable<L>,
    ) -> Result<Symbolic<'_, L>, F4Error> {
        self.pivot_cols.sort_unstable_by(|a, b| table.cmp(*b, *a));
        self.free_cols.sort_unstable_by(|a, b| table.cmp(*b, *a));
        let npiv = self.pivot_cols.len() as u32;
        let mut columns = self.pivot_cols;
        reserve(&mut columns, self.free_cols.len())?;
        columns.append(&mut self.free_cols);

        let mut col_of = Vec::new();
        reserve(&mut col_of, table.len())?;
        col_of.resize(table.len(), NO_COLUMN);
        for (col, &m) in columns.iter().enumerate() {
            col_of[m.index()] = col as u32;
        }

        let mut pivot_row_of = Vec::new();
        reserve(&mut pivot_row_of, npiv as usize)?;
        for &m in &columns[..npiv as usize] {
            let row = self.pivot_row[m.index()];
            debug_assert_ne!(row, NO_ROW, "every pivot column has a reducer row");
            pivot_row_of.push(row);
        }

        for (row, lead) in self.rows.iter_mut().zip(&self.row_lead) {
            row.lead_col = col_of[lead.index()];
            debug_assert!(row.lead_col < npiv, "a row's lead is a pivot column");
        }

        Ok(Symbolic {
            table,
            rows: self.rows,
            columns,
            npiv,
            col_of,
            pivot_row_of,
            row_monos: self.row_monos,
            row_start: self.row_start,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compute::f4::basis::BasisPoly;
    use crate::compute::f4::monomial::Lanes8;
    use crate::compute::f4::pairs::{PairSet, SelectOptions};

    /// A seeded basis and the tables one batch needs.
    struct Fixture {
        basis: Basis<Lanes8>,
        pairs: PairSet,
        batch: MonomialTable<Lanes8>,
    }

    /// Build a basis from monic polynomials, each a list of exponent
    /// vectors in descending order.
    fn build(nvars: usize, polys: &[&[&[u32]]]) -> Fixture {
        let mut table = MonomialTable::<Lanes8>::new(nvars);
        let interned: Vec<Vec<MonomialId>> = polys
            .iter()
            .map(|poly| poly.iter().map(|e| table.intern(e).unwrap()).collect())
            .collect();
        let mut basis = Basis::new(table);
        let mut pairs = PairSet::new();
        let mut ws = MonomialTable::new(nvars);
        for monos in interned {
            let coeffs = vec![1; monos.len()];
            basis
                .insert(
                    BasisPoly {
                        monos,
                        coeffs,
                        shoup: Vec::new(),
                    },
                    &mut pairs,
                    &mut ws,
                    &mut Deadline::none(),
                )
                .unwrap();
        }
        Fixture {
            basis,
            pairs,
            batch: MonomialTable::new(nvars),
        }
    }

    impl Fixture {
        /// Preprocess the lowest degree class.
        fn preprocess(&mut self) -> Symbolic<'_, Lanes8> {
            let batch = self
                .pairs
                .take_lowest_degree(self.basis.table(), &SelectOptions::default())
                .expect("the basis has a pair");
            preprocess(
                &batch,
                &self.basis,
                &mut self.batch,
                &Strategy::default(),
                &mut Deadline::none(),
            )
            .unwrap()
        }
    }

    fn exps<L: Lanes>(table: &MonomialTable<L>, id: MonomialId) -> Vec<u32> {
        let mut out = Vec::new();
        table.unpack(id, &mut out);
        out
    }

    fn columns<L: Lanes>(symbolic: &Symbolic<'_, L>) -> Vec<Vec<u32>> {
        symbolic
            .columns
            .iter()
            .map(|&m| exps(symbolic.table, m))
            .collect()
    }

    /// Check the invariants every result carries.
    fn check<L: Lanes>(symbolic: &Symbolic<'_, L>) {
        let npiv = symbolic.npiv as usize;
        check_column_order(symbolic, npiv);
        check_column_map(symbolic);
        check_pivot_rows(symbolic);
        check_row_leads(symbolic);
        assert_eq!(symbolic.pivot_row_of.len(), npiv);
        assert_eq!(symbolic.nfree() as usize, symbolic.columns.len() - npiv);
    }

    fn check_column_order<L: Lanes>(symbolic: &Symbolic<'_, L>, npiv: usize) {
        for block in [&symbolic.columns[..npiv], &symbolic.columns[npiv..]] {
            for window in block.windows(2) {
                assert_eq!(
                    symbolic.table.cmp(window[0], window[1]),
                    std::cmp::Ordering::Greater,
                    "each block descends"
                );
            }
        }
    }

    fn check_column_map<L: Lanes>(symbolic: &Symbolic<'_, L>) {
        for (col, &m) in symbolic.columns.iter().enumerate() {
            assert_eq!(symbolic.column_of(m), Some(col as u32));
        }
    }

    fn check_pivot_rows<L: Lanes>(symbolic: &Symbolic<'_, L>) {
        for (col, &row) in symbolic.pivot_row_of.iter().enumerate() {
            assert_eq!(symbolic.rows[row as usize].lead_col, col as u32);
            assert!(symbolic.is_pivot_row(row));
        }
    }

    fn check_row_leads<L: Lanes>(symbolic: &Symbolic<'_, L>) {
        for row in 0..symbolic.rows.len() {
            assert!(symbolic.is_pivot_column(symbolic.rows[row].lead_col));
        }
    }

    #[test]
    fn the_m1_baseline_is_first_and_smallest_source() {
        let strategy = Strategy::default();
        assert_eq!(strategy.reducer, Reducer::First);
        assert_eq!(strategy.upper_row, UpperRow::SmallestSource);
    }

    #[test]
    fn one_pair_gives_two_rows_and_one_pivot_column() {
        // f = x^2 + y and g = x*y + 1. The pair has lcm x^2*y, so the rows
        // are y*f and x*g, and neither tail monomial has a reducer.
        let mut fixture = build(2, &[&[&[2, 0], &[0, 1]], &[&[1, 1], &[0, 0]]]);
        let symbolic = fixture.preprocess();
        check(&symbolic);

        assert_eq!(symbolic.rows.len(), 2);
        assert_eq!(symbolic.npiv, 1);
        assert_eq!(
            columns(&symbolic),
            vec![vec![2, 1], vec![0, 2], vec![1, 0]],
            "x^2*y, then y^2 and x"
        );
        assert_eq!(symbolic.rows[0].source, 0);
        assert_eq!(exps(symbolic.table, symbolic.rows[0].mult), vec![0, 1]);
        assert_eq!(symbolic.rows[1].source, 1);
        assert_eq!(exps(symbolic.table, symbolic.rows[1].mult), vec![1, 0]);
        assert_eq!(symbolic.pivot_row_of, vec![0], "the smallest source");
    }

    #[test]
    fn a_tail_monomial_with_a_divisor_gets_a_reducer_row() {
        // f = x^2*y + x*y^2 and g = x*y^2 + 1. The pair has lcm x^2*y^2,
        // so the rows are y*f and x*g. The tail monomial x*y^3 of y*f is a
        // multiple of lm(g), which gives it the reducer row y*g.
        let mut fixture = build(2, &[&[&[2, 1], &[1, 2]], &[&[1, 2], &[0, 0]]]);
        let symbolic = fixture.preprocess();
        check(&symbolic);

        assert_eq!(symbolic.npiv, 2);
        assert_eq!(symbolic.rows.len(), 3);
        assert_eq!(symbolic.nfree(), 2);
        assert_eq!(symbolic.rows[2].source, 1);
        assert_eq!(exps(symbolic.table, symbolic.rows[2].mult), vec![0, 1]);
        assert_eq!(symbolic.pivot_row_of[1], 2);
    }

    #[test]
    fn two_sources_with_one_leading_monomial_both_emit_a_row() {
        // Inserting x*y + y retires it under x*y + x, and the pair of the
        // two equal leads stays in the queue. Deduplication is on the
        // source index, so both rows appear with the same multiplier.
        let mut fixture = build(2, &[&[&[1, 1], &[1, 0]], &[&[1, 1], &[0, 1]]]);
        let symbolic = fixture.preprocess();
        check(&symbolic);

        assert_eq!(symbolic.rows.len(), 2);
        assert_eq!(symbolic.npiv, 1);
        assert_eq!(symbolic.rows[0].source, 0);
        assert_eq!(symbolic.rows[1].source, 1);
        assert_eq!(symbolic.rows[0].mult, MonomialId::ONE);
        assert_eq!(symbolic.rows[1].mult, MonomialId::ONE);
    }

    #[test]
    fn two_runs_give_the_same_rows_and_columns() {
        let system: &[&[&[u32]]] = &[
            &[&[2, 0, 0], &[0, 1, 0]],
            &[&[1, 1, 0], &[0, 0, 1]],
            &[&[0, 2, 0], &[1, 0, 0]],
        ];
        let mut first = build(3, system);
        let mut second = build(3, system);
        let one = first.preprocess();
        let two = second.preprocess();

        assert_eq!(columns(&one), columns(&two));
        assert_eq!(one.npiv, two.npiv);
        assert_eq!(one.pivot_row_of, two.pivot_row_of);
        let rows_of = |symbolic: &Symbolic<'_, Lanes8>| -> Vec<(u32, Vec<u32>, u32)> {
            symbolic
                .rows
                .iter()
                .map(|row| (row.source, exps(symbolic.table, row.mult), row.lead_col))
                .collect()
        };
        assert_eq!(rows_of(&one), rows_of(&two));
    }

    #[test]
    fn the_shortest_upper_row_rule_picks_the_shorter_source() {
        // Both elements have the leading monomial x*y. The first has three
        // terms and the second has two, so the two rules disagree.
        let mut fixture = build(2, &[&[&[1, 1], &[1, 0], &[0, 0]], &[&[1, 1], &[0, 1]]]);
        let batch = fixture
            .pairs
            .take_lowest_degree(fixture.basis.table(), &SelectOptions::default())
            .unwrap();
        let strategy = Strategy {
            reducer: Reducer::First,
            upper_row: UpperRow::ShortestSource,
        };
        let symbolic = preprocess(
            &batch,
            &fixture.basis,
            &mut fixture.batch,
            &strategy,
            &mut Deadline::none(),
        )
        .unwrap();
        check(&symbolic);

        assert_eq!(symbolic.pivot_row_of, vec![1]);
        assert_eq!(symbolic.rows[1].source, 1);
    }
}
