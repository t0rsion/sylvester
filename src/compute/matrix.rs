//! Macaulay-matrix F5: critical pairs batched by degree, reduced by
//! F4-style sparse elimination under the F5 syzygy criterion.

use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::mem::size_of;
use std::time::Instant;

use super::interreduce::reduced_groebner_basis_checked;
use super::signature::{LabeledPoly, Signature, add_syzygy_rule, is_sig_redundant, is_syzygy};
use super::{ComputeError, REDUCE_DEADLINE_STRIDE, poll_deadline};
use crate::poly::{Monomial, Polynomial, Term, key_divides};
use crate::ring::PolynomialRing;
use crate::ring::field::{Felt, Modulus};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// The row count between resource checks inside one degree batch's
/// elimination.
///
/// `Instant::now()` and a byte count both cost something, so checking
/// every row would tax elimination itself. This stride keeps that cost
/// off the hot path while still bounding how far one batch runs an
/// exhausted budget past its limit.
const ROW_CHECK_STRIDE: usize = 64;

/// The smallest group of equal-signature rows that makes parallel
/// reduction worth its coordination cost.
///
/// Rows that share one signature reduce against the same pivots, so the
/// `parallel` feature reduces such a group with a rayon loop once it
/// reaches this many rows. Below it, the single-threaded loop wins. Three
/// is a measured value: on katsura-7, katsura-8, cyclic-6, and noon-5 at
/// p = 1073741827 it beats 2, 4, 6, 8, and 16.
#[cfg(feature = "parallel")]
const PARALLEL_ROW_THRESHOLD: usize = 3;

#[derive(Clone, Debug)]
struct CriticalPair {
    i: usize,
    j: usize,
    sig: Signature,
    origin: usize,
}

#[derive(Default)]
struct PairManager {
    by_degree: BTreeMap<u32, Vec<CriticalPair>>,
}

impl PairManager {
    fn push_pair(&mut self, basis: &[LabeledPoly], i: usize, j: usize) -> Result<(), ComputeError> {
        let Some(lm_i) = basis[i].poly.lm() else {
            return Ok(());
        };
        let Some(lm_j) = basis[j].poly.lm() else {
            return Ok(());
        };
        let lcm = lm_i.lcm(lm_j);
        if lcm.deg > super::DEGREE_LIMIT {
            return Err(ComputeError::DegreeLimit {
                limit: super::DEGREE_LIMIT,
            });
        }
        // lcm is a multiple of lm_i and lm_j by definition.
        let m_i = lcm.quotient(lm_i).expect("lm_i divides lcm");
        // lcm is a multiple of lm_i and lm_j by definition.
        let m_j = lcm.quotient(lm_j).expect("lm_j divides lcm");
        // A signature's degree is unrelated to the lcm's, so the lcm gate
        // above does not bound this product and it reports the limit itself.
        let sig_i = Signature {
            index: basis[i].sig.index,
            term: basis[i].sig.term.checked_mul(&m_i)?,
        };
        let sig_j = Signature {
            index: basis[j].sig.index,
            term: basis[j].sig.term.checked_mul(&m_j)?,
        };
        let (sig, origin) = match sig_i.cmp(&sig_j) {
            Ordering::Greater => (sig_i, i),
            Ordering::Less => (sig_j, j),
            // Non-regular pair: the component signatures agree as module
            // monomials (coefficients are not tracked, so this covers both
            // the singular case, where the module leading terms cancel and
            // the S-polynomial's true signature is strictly smaller, and
            // the super-regular case). Either way the common value is not a
            // trusted signature for the S-polynomial: keeping the pair
            // could mislabel it and, on a zero reduction, record an
            // inflated syzygy signature that later discards necessary
            // pairs. The signature-based Buchberger criterion only requires
            // regular S-pairs, so rejecting is sound.
            Ordering::Equal => return Ok(()),
        };
        let degree = lcm.deg;
        self.by_degree
            .entry(degree)
            .or_default()
            .push(CriticalPair { i, j, sig, origin });
        Ok(())
    }

    fn pop_smallest_degree(&mut self) -> Option<(u32, Vec<CriticalPair>)> {
        let degree = *self.by_degree.keys().next()?;
        let pairs = self.by_degree.remove(&degree).unwrap_or_default();
        Some((degree, pairs))
    }
}

fn check_limits(
    budget: RunBudget,
    basis: &[LabeledPoly],
    syzygy_rules: &[Signature],
    pair_manager: &PairManager,
    rewriters: &[Vec<RewriteRule>],
    extra_pairs: usize,
    nvars: usize,
) -> Result<(), ComputeError> {
    let RunBudget {
        deadline,
        max_memory_bytes,
    } = budget;
    poll_deadline(deadline)?;

    if let Some(limit) = max_memory_bytes {
        let per_term = super::per_term_bytes(nvars);
        let mut bytes = 0usize;
        bytes = bytes.saturating_add(std::mem::size_of_val(basis));
        bytes = bytes.saturating_add(
            basis
                .iter()
                .map(|lp| lp.poly.terms.len() * per_term)
                .sum::<usize>(),
        );
        bytes = bytes.saturating_add(std::mem::size_of_val(syzygy_rules));
        let pair_count = extra_pairs
            + pair_manager
                .by_degree
                .values()
                .map(|pairs| pairs.len())
                .sum::<usize>();
        bytes = bytes.saturating_add(pair_count * size_of::<CriticalPair>());
        let rewrite_count: usize = rewriters.iter().map(|rules| rules.len()).sum();
        bytes = bytes
            .saturating_add(std::mem::size_of_val(rewriters))
            .saturating_add(rewrite_count * size_of::<RewriteRule>());
        if bytes > limit {
            return Err(ComputeError::MemoryLimitExceeded);
        }
    }

    Ok(())
}

#[derive(Clone, Debug)]
struct MatrixRow {
    sig: Signature,
    cols: Vec<usize>,
    coeffs: Vec<Felt>,
}

impl MatrixRow {
    fn leading_col(&self) -> Option<usize> {
        self.cols.first().copied()
    }

    fn is_zero(&self) -> bool {
        self.cols.is_empty()
    }
}

#[derive(Clone, Debug)]
struct RewriteRule {
    term: Monomial,
    basis_index: usize,
}

/// The canonical rewriter of a signature: the most recently added basis
/// element whose signature divides it, if that element was added after the
/// pair's own origin. Returns the basis index of that element.
///
/// F5's rewritten criterion must not simply delete such a pair: with pairs
/// batched by degree the signature-order induction that deletion relies on
/// does not hold, and deleting loses necessary S-polynomials. Instead the
/// pair's row is replaced by the rewriter's multiple, which carries exactly
/// the same signature, so the signature is still handled.
fn canonical_rewriter(
    sig: &Signature,
    origin: usize,
    rewriters: &[Vec<RewriteRule>],
) -> Option<usize> {
    rewriters.get(sig.index).and_then(|rules| {
        rules
            .iter()
            .filter(|rule| rule.basis_index > origin && rule.term.divides(&sig.term))
            .map(|rule| rule.basis_index)
            .max()
    })
}

fn add_rewriter(rewriters: &mut [Vec<RewriteRule>], sig: &Signature, basis_index: usize) {
    let rules = &mut rewriters[sig.index];
    if rules
        .iter()
        .any(|rule| rule.basis_index > basis_index && rule.term.divides(&sig.term))
    {
        return;
    }
    rules.retain(|rule| !(sig.term.divides(&rule.term) && rule.basis_index < basis_index));
    rules.push(RewriteRule {
        term: sig.term.clone(),
        basis_index,
    });
}

/// One row of the matrix held as a dense accumulator over the columns.
///
/// Every entry is a representative in `[0, p)`, and every entry below
/// `lead` is zero. Reduction only ever writes at or above the column it
/// cancels, so `lead` never moves back. A ring modulus is at most
/// 2^31 - 1, so an entry fits a `u32` and a batch's accumulator takes half
/// the cache it would take in `u64`.
struct Accumulator {
    values: Vec<u32>,
    lead: usize,
    modulus: Modulus,
}

impl Accumulator {
    fn new(num_columns: usize, modulus: Modulus) -> Self {
        Accumulator {
            values: vec![0; num_columns],
            lead: 0,
            modulus,
        }
    }

    fn scatter(&mut self, cols: &[usize], coeffs: &[Felt]) {
        for (&col, &coeff) in cols.iter().zip(coeffs) {
            // a coefficient is below the modulus, which is below 2^31.
            self.values[col] = coeff.value() as u32;
        }
        self.lead = cols.first().copied().unwrap_or(self.values.len());
    }

    /// Subtract `factor * pivot` from the row.
    ///
    /// The pivot's own column holds `factor * lead(pivot)` before the
    /// step, so the step clears it. Every entry stays the representative
    /// in `[0, p)`, which the accumulator relies on: a column counts as
    /// occupied exactly when its entry is not zero.
    #[inline]
    fn sub_scaled(&mut self, pivot: &Pivot, factor: u64) {
        let modulus = self.modulus;
        // a ring modulus is below 2^31, so the sum below stays in a u32.
        let p = modulus.value() as u32;
        for (&col, &coeff) in pivot.cols.iter().zip(&pivot.coeffs) {
            let t = modulus.mul(factor, coeff.value()) as u32;
            let value = self.values[col] + (p - t);
            self.values[col] = if value >= p { value - p } else { value };
        }
    }

    /// Advance `lead` to the first occupied column, or past the last
    /// column when the row is zero.
    #[inline]
    fn seek(&mut self) -> Option<usize> {
        while self.lead < self.values.len() && self.values[self.lead] == 0 {
            self.lead += 1;
        }
        (self.lead < self.values.len()).then_some(self.lead)
    }

    /// Move the row out as a sparse column list and leave the accumulator
    /// zero.
    fn gather(&mut self) -> (Vec<usize>, Vec<Felt>) {
        let mut cols: Vec<usize> = Vec::new();
        let mut coeffs: Vec<Felt> = Vec::new();
        for col in self.lead..self.values.len() {
            let value = self.values[col];
            if value != 0 {
                cols.push(col);
                coeffs.push(Felt::from_residue(value as u64));
                self.values[col] = 0;
            }
        }
        (cols, coeffs)
    }
}

/// A row that owns its leading column, with the inverse of its leading
/// coefficient.
struct Pivot {
    sig: Signature,
    cols: Vec<usize>,
    coeffs: Vec<Felt>,
    inv: Felt,
}

/// The columns of one matrix, largest monomial first, with the column of
/// each.
///
/// One table serves both: the pass that collects the monomials is the pass
/// that ends up mapping them, so no monomial is hashed twice. The table
/// grows one row at a time and is charged to `charge` at the
/// [`ROW_CHECK_STRIDE`] cadence, so a batch whose columns alone pass the
/// limit stops there instead of after the whole set is built.
fn collect_columns(
    rows: &[(&Signature, &Polynomial)],
    charge: &SymbolicCharge,
) -> Result<(Vec<Monomial>, FxHashMap<Monomial, usize>), ComputeError> {
    let mut col_map: FxHashMap<Monomial, usize> = FxHashMap::default();
    let mut since_check = 0usize;
    for (_, poly) in rows {
        for term in &poly.terms {
            if !col_map.contains_key(&term.mono) {
                col_map.insert(term.mono.clone(), 0);
            }
        }
        since_check += 1;
        if since_check >= ROW_CHECK_STRIDE {
            since_check = 0;
            charge.check(charge.bytes(0, 0, col_map.len()))?;
        }
    }
    charge.check(charge.bytes(0, 0, col_map.len()))?;
    let mut columns: Vec<Monomial> = col_map.keys().cloned().collect();
    columns.sort_by(|a, b| b.cmp(a));
    for (index, mono) in columns.iter().enumerate() {
        // every column came out of the table.
        *col_map.get_mut(mono).expect("the monomial is a column") = index;
    }
    Ok((columns, col_map))
}

fn poly_to_row(poly: &Polynomial, col_map: &FxHashMap<Monomial, usize>) -> (Vec<usize>, Vec<Felt>) {
    // col_map is built from all monomials appearing in the matrix.
    let mut entries: Vec<(usize, Felt)> = poly
        .terms
        .iter()
        .map(|t| {
            (
                *col_map.get(&t.mono).expect("monomial in column map"),
                t.coeff,
            )
        })
        .collect();
    entries.sort_by_key(|(col, _)| *col);
    let mut cols = Vec::with_capacity(entries.len());
    let mut coeffs = Vec::with_capacity(entries.len());
    for (col, coeff) in entries {
        cols.push(col);
        coeffs.push(coeff);
    }
    (cols, coeffs)
}

fn row_to_poly(
    ring: &PolynomialRing,
    cols: &[usize],
    coeffs: &[Felt],
    columns: &[Monomial],
) -> Polynomial {
    let terms = cols
        .iter()
        .zip(coeffs)
        .map(|(col, coeff)| Term {
            coeff: *coeff,
            mono: columns[*col].clone(),
        })
        .collect();
    Polynomial::from_terms(ring.clone(), terms)
}

fn make_pivot_rows_monic(pivot_rows: &mut [MatrixRow], p: u64) {
    if pivot_rows.is_empty() {
        return;
    }

    let lead_coeffs: Vec<Felt> = pivot_rows
        .iter()
        .map(|row| row.coeffs.first().copied().unwrap_or_else(Felt::zero))
        .collect();
    let inverses = Felt::batch_inv(&lead_coeffs, p);

    for (row, inv) in pivot_rows.iter_mut().zip(inverses) {
        if row.is_zero() || inv.is_zero() {
            continue;
        }
        for coeff in &mut row.coeffs {
            *coeff = coeff.mul(inv, p);
        }
    }
}

/// Reduce one row against the pivots and gather what survives.
///
/// The result is the sparse form of the reduced row, or `None` when the
/// row reduces to zero.
///
/// Two rules fix which reductions happen, and both are needed for a
/// signature-safe result. A pivot may reduce a row only when its
/// signature is strictly smaller. Those pivots are the prefix `eligible`
/// of the pivot list, because the rows arrive in ascending signature
/// order. The pivots then run in creation order, not in column order: a
/// pivot can write back into the leading column of a pivot created
/// earlier, so another visiting order would cancel entries this one
/// keeps. After that sweep only the leading column is reduced, and only
/// while the pivot that owns it is eligible.
fn reduce_row(
    row: &MatrixRow,
    pivots: &[Pivot],
    pivot_of_col: &[Option<usize>],
    eligible: usize,
    acc: &mut Accumulator,
) -> Option<(Vec<usize>, Vec<Felt>)> {
    acc.scatter(&row.cols, &row.coeffs);

    for pivot in &pivots[..eligible] {
        // A pivot owns its leading column, so cols is not empty.
        let value = acc.values[pivot.cols[0]];
        if value != 0 {
            let factor = acc.modulus.mul(value as u64, pivot.inv.value());
            acc.sub_scaled(pivot, factor);
        }
    }

    loop {
        let lead = acc.seek()?;
        match pivot_of_col[lead] {
            Some(index) if index < eligible => {
                let pivot = &pivots[index];
                let factor = acc.modulus.mul(acc.values[lead] as u64, pivot.inv.value());
                acc.sub_scaled(pivot, factor);
            }
            _ => return Some(acc.gather()),
        }
    }
}

/// The number of rows at the front of `rows` that share one signature.
fn signature_group(rows: &[MatrixRow]) -> usize {
    let Some(first) = rows.first() else { return 0 };
    1 + rows[1..]
        .iter()
        .take_while(|row| row.sig == first.sig)
        .count()
}

/// The live bytes one degree batch's matrix holds before elimination.
///
/// `rows` and `terms` count the polynomial rows the symbolic phase built
/// and their terms; the matrix holds one sparse entry per term. The rest is
/// the dense accumulator, the sorted column list, and the column-to-index
/// table, one entry each per column. The count is taken from the counts,
/// not from the built matrix, so the caller can check it before it
/// allocates.
fn batch_bytes(rows: usize, terms: usize, num_columns: usize, nvars: usize) -> usize {
    let heap_exps = crate::poly::heap_exps_bytes(nvars);
    let mut bytes = rows.saturating_mul(size_of::<MatrixRow>());
    bytes = bytes.saturating_add(terms.saturating_mul(size_of::<usize>() + size_of::<Felt>()));
    // the dense accumulator holds one u32 per column.
    bytes = bytes.saturating_add(num_columns.saturating_mul(size_of::<u32>()));
    bytes.saturating_add(
        num_columns.saturating_mul(size_of::<Monomial>() + heap_exps + size_of::<usize>()),
    )
}

/// The memory limit as the symbolic phase of one degree batch sees it.
///
/// The phase builds every row and every column monomial before elimination
/// starts, so a limit checked only on the whole batch is a limit the batch
/// has already passed. This charges what is built so far, at the
/// [`ROW_CHECK_STRIDE`] cadence the deadline uses, which bounds the
/// overshoot by one stride instead of one batch.
#[derive(Clone, Copy)]
struct SymbolicCharge {
    limit: Option<usize>,
    per_term: usize,
    per_mono: usize,
    held: usize,
}

impl SymbolicCharge {
    fn new(limit: Option<usize>, nvars: usize) -> Self {
        SymbolicCharge {
            limit,
            per_term: super::per_term_bytes(nvars),
            per_mono: size_of::<Monomial>() + crate::poly::heap_exps_bytes(nvars),
            held: 0,
        }
    }

    /// The bytes `rows` rows of `terms` terms together with `monomials`
    /// interned monomials cost.
    fn bytes(&self, rows: usize, terms: usize, monomials: usize) -> usize {
        rows.saturating_mul(size_of::<Signature>() + size_of::<Polynomial>())
            .saturating_add(terms.saturating_mul(self.per_term))
            .saturating_add(monomials.saturating_mul(self.per_mono))
    }

    /// Return an error when `bytes` on top of what is held passes the limit.
    fn check(&self, bytes: usize) -> Result<(), ComputeError> {
        match self.limit {
            Some(limit) if self.held.saturating_add(bytes) > limit => {
                Err(ComputeError::MemoryLimitExceeded)
            }
            _ => Ok(()),
        }
    }

    /// Charge `bytes` and return an error when the limit is passed.
    fn hold(&mut self, bytes: usize) -> Result<(), ComputeError> {
        self.held = self.held.saturating_add(bytes);
        self.check(0)
    }
}

/// The live bytes the pivots built so far hold, plus the dense
/// accumulator they reduce against.
fn pivot_bytes(pivots: &[Pivot], num_columns: usize) -> usize {
    let mut bytes = std::mem::size_of_val(pivots);
    bytes = bytes.saturating_add(
        pivots
            .iter()
            .map(|pivot| pivot.cols.len() * (size_of::<usize>() + size_of::<Felt>()))
            .sum::<usize>(),
    );
    bytes.saturating_add(num_columns * size_of::<u32>())
}

fn eliminate_rows(
    mut rows: Vec<MatrixRow>,
    p: u64,
    num_columns: usize,
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
) -> Result<(Vec<MatrixRow>, Vec<Signature>), ComputeError> {
    rows.sort_by(|a, b| {
        a.sig
            .cmp(&b.sig)
            .then_with(|| a.leading_col().cmp(&b.leading_col()))
    });

    let modulus = Modulus::new(p);
    let mut acc = Accumulator::new(num_columns, modulus);
    let mut pivots: Vec<Pivot> = Vec::new();
    let mut pivot_of_col: Vec<Option<usize>> = vec![None; num_columns];
    let mut zero_sigs: Vec<Signature> = Vec::new();
    // Pivot signatures rise with the row order, so the pivots one row may
    // use are a prefix of the list, and the bound never moves back.
    let mut eligible = 0usize;
    let mut since_check = 0usize;

    let mut start = 0usize;
    while start < rows.len() {
        let group = &rows[start..start + signature_group(&rows[start..])];
        while eligible < pivots.len() && pivots[eligible].sig < group[0].sig {
            eligible += 1;
        }

        // Every row of the group reduces against the same pivots: a pivot
        // this group creates has the group's signature, which is not
        // strictly smaller, so no row of the group may use it.
        #[cfg(feature = "parallel")]
        let reduced: Vec<Option<(Vec<usize>, Vec<Felt>)>> = if group.len() >= PARALLEL_ROW_THRESHOLD
        {
            group
                .par_iter()
                .map_init(
                    || Accumulator::new(num_columns, modulus),
                    |acc, row| reduce_row(row, &pivots, &pivot_of_col, eligible, acc),
                )
                .collect()
        } else {
            group
                .iter()
                .map(|row| reduce_row(row, &pivots, &pivot_of_col, eligible, &mut acc))
                .collect()
        };
        #[cfg(not(feature = "parallel"))]
        let reduced: Vec<Option<(Vec<usize>, Vec<Felt>)>> = group
            .iter()
            .map(|row| reduce_row(row, &pivots, &pivot_of_col, eligible, &mut acc))
            .collect();

        for (row, result) in group.iter().zip(reduced) {
            let Some((cols, coeffs)) = result else {
                zero_sigs.push(row.sig.clone());
                continue;
            };
            // The gathered row is not empty, so it has a leading column.
            let lead = cols[0];
            if pivot_of_col[lead].is_some() {
                continue;
            }
            let inv = coeffs[0].inv(p);
            pivot_of_col[lead] = Some(pivots.len());
            pivots.push(Pivot {
                sig: row.sig.clone(),
                cols,
                coeffs,
                inv,
            });
        }

        start += group.len();
        since_check += group.len();
        if since_check >= ROW_CHECK_STRIDE {
            since_check = 0;
            poll_deadline(deadline)?;
            if let Some(limit) = max_memory_bytes {
                if pivot_bytes(&pivots, num_columns) > limit {
                    return Err(ComputeError::MemoryLimitExceeded);
                }
            }
        }
    }

    let mut pivot_rows: Vec<MatrixRow> = pivots
        .into_iter()
        .map(|pivot| MatrixRow {
            sig: pivot.sig,
            cols: pivot.cols,
            coeffs: pivot.coeffs,
        })
        .collect();
    make_pivot_rows_monic(&mut pivot_rows, p);

    Ok((pivot_rows, zero_sigs))
}

/// The degree and the divisor key of one basis lead, with the place of
/// that element in the basis.
struct LeadFilter {
    deg: u32,
    key: u64,
    basis_pos: u32,
}

/// The basis, with the lead filter of every element beside it.
///
/// Symbolic preprocessing and [`f5_reduce`] both test the whole basis
/// against every monomial they meet, so the test that rejects a candidate
/// must not reach into the basis. The filters sit in one array here, and
/// only a candidate that passes reads its monomial. Growing the basis
/// through [`Basis::push`] is what keeps the two in step.
#[derive(Default)]
struct Basis {
    elements: Vec<LabeledPoly>,
    filters: Vec<LeadFilter>,
}

impl Basis {
    fn push(&mut self, element: LabeledPoly) {
        if let Some(lm) = element.poly.lm() {
            self.filters.push(LeadFilter {
                deg: lm.deg,
                key: lm.divisor_key(),
                // a basis holds far fewer than 2^32 elements.
                basis_pos: self.elements.len() as u32,
            });
        }
        self.elements.push(element);
    }

    fn filters(&self) -> &[LeadFilter] {
        &self.filters
    }

    fn into_polynomials(self) -> Vec<Polynomial> {
        self.elements.into_iter().map(|lp| lp.poly).collect()
    }
}

impl std::ops::Deref for Basis {
    type Target = [LabeledPoly];

    fn deref(&self) -> &[LabeledPoly] {
        &self.elements
    }
}

/// Build the reducer rows of one matrix.
///
/// Every monomial of the matrix that a basis lead divides gets one row per
/// such lead, which is the multiple that cancels it. The rows come out in
/// the order the scan finds them, and `eliminate_rows` sorts them, so the
/// scan must keep basis order.
fn build_reducer_rows(
    basis: &Basis,
    mut monomials: FxHashSet<Monomial>,
    deadline: Option<Instant>,
    charge: &SymbolicCharge,
) -> Result<(Vec<(Signature, Polynomial)>, usize), ComputeError> {
    let mut rows: Vec<(Signature, Polynomial)> = Vec::new();
    let mut row_terms = 0usize;
    let mut worklist: Vec<Monomial> = monomials.iter().cloned().collect();
    let filters = basis.filters();
    let mut since_check = 0usize;

    // A monomial reaches the worklist once, so one (basis element,
    // monomial) pair is handled once and no row can repeat.
    while let Some(mono) = worklist.pop() {
        since_check += 1;
        if since_check >= ROW_CHECK_STRIDE {
            since_check = 0;
            poll_deadline(deadline)?;
            charge.check(charge.bytes(rows.len(), row_terms, monomials.len()))?;
        }
        let mono_key = mono.divisor_key();
        for filter in filters {
            if filter.deg > mono.deg || !key_divides(filter.key, mono_key) {
                continue;
            }
            let g = &basis[filter.basis_pos as usize];
            // a filter is only recorded for an element that has a lead.
            let lm_g = g.poly.lm().expect("a filtered element has a lead");
            if !lm_g.divides(&mono) {
                continue;
            }
            // divides() implies a quotient exists.
            let m = mono.quotient(lm_g).expect("divides() implies quotient()");

            // A signature's degree is unrelated to the monomial's, so this
            // product reports the limit itself.
            let sig = Signature {
                index: g.sig.index,
                term: g.sig.term.checked_mul(&m)?,
            };
            let poly = g.poly.shift_monomial(&m)?;
            for term in &poly.terms {
                if !monomials.contains(&term.mono) {
                    monomials.insert(term.mono.clone());
                    worklist.push(term.mono.clone());
                }
            }
            row_terms += poly.terms.len();
            rows.push((sig, poly));
        }
    }

    charge.check(charge.bytes(rows.len(), row_terms, monomials.len()))?;
    Ok((rows, row_terms))
}

/// Post-process one surviving pivot row: finish any regular (signature-safe)
/// reduction the matrix could not perform, record a syzygy on a zero result,
/// drop sig-redundant results, and otherwise insert the row into the basis.
/// Rows whose lead is divisible by a basis lead but not covered at the
/// signature level carry genuinely new information and must be kept;
/// dropping them loses the Gröbner property of the output.
///
/// Returns the basis index of an inserted row, or `None` when the row was a
/// syzygy or sig-redundant.
fn insert_candidate(
    sig: Signature,
    poly: Polynomial,
    basis: &mut Basis,
    syzygy_rules: &mut Vec<Signature>,
    rewriters: &mut [Vec<RewriteRule>],
    modulus: u64,
    deadline: Option<Instant>,
) -> Result<Option<usize>, ComputeError> {
    let labeled = LabeledPoly {
        sig,
        poly,
        index: basis.len(),
    };
    let mut reduced = f5_reduce(labeled, basis, modulus, deadline)?;
    if reduced.poly.is_zero() {
        add_syzygy_rule(syzygy_rules, reduced.sig);
        return Ok(None);
    }
    if is_sig_redundant(&reduced, basis) {
        return Ok(None);
    }
    reduced.poly = reduced.poly.make_monic(modulus);
    add_rewriter(rewriters, &reduced.sig, reduced.index);
    let index = reduced.index;
    basis.push(reduced);
    Ok(Some(index))
}

/// The deadline and the memory cap one degree batch runs under.
#[derive(Clone, Copy)]
struct RunBudget {
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
}

fn process_degree(
    ring: &PolynomialRing,
    pairs: Vec<CriticalPair>,
    basis: &mut Basis,
    syzygy_rules: &mut Vec<Signature>,
    rewriters: &mut [Vec<RewriteRule>],
    pair_manager: &mut PairManager,
    budget: RunBudget,
) -> Result<(), ComputeError> {
    let RunBudget {
        deadline,
        max_memory_bytes,
    } = budget;
    let modulus = ring.modulus();
    let nvars = ring.nvars();
    let mut charge = SymbolicCharge::new(max_memory_bytes, nvars);
    let mut s_rows: Vec<(Signature, Polynomial)> = Vec::new();
    let mut s_terms = 0usize;
    let mut since_check = 0usize;
    for pair in pairs {
        since_check += 1;
        if since_check >= ROW_CHECK_STRIDE {
            since_check = 0;
            charge.check(charge.bytes(s_rows.len(), s_terms, 0))?;
        }
        if is_syzygy(&pair.sig, syzygy_rules) {
            continue;
        }
        if let Some(rewriter) = canonical_rewriter(&pair.sig, pair.origin, rewriters) {
            // Replace the S-polynomial row with the canonical rewriter's
            // multiple R = m * g. Both the original S-row P and R have
            // signature monomial pair.sig, so after scaling to cancel the
            // common module leading term, P - cR has strictly smaller
            // signature and it suffices to handle R: elimination only
            // subtracts rows of strictly smaller signature, so R keeps its
            // exact signature; a zero result records an exact syzygy; a
            // surviving row goes through insert_candidate, which inserts it
            // unless it is sig-redundant (in which case a basis element
            // already covers its signature region and lead). Bare deletion
            // of the pair, by contrast, is unsound under degree batching.
            let g = &basis[rewriter];
            // the rewrite rule's term is g's signature term, and
            // divisibility was checked by canonical_rewriter.
            let m = pair
                .sig
                .term
                .quotient(&g.sig.term)
                .expect("rewriter signature divides pair signature");
            let row = g.poly.shift_monomial(&m)?;
            if !row.is_zero() {
                s_terms += row.terms.len();
                s_rows.push((pair.sig, row));
            }
            continue;
        }
        let s = basis[pair.i]
            .poly
            .s_polynomial(&basis[pair.j].poly, modulus)?;
        if s.is_zero() {
            add_syzygy_rule(syzygy_rules, pair.sig);
            continue;
        }
        s_terms += s.terms.len();
        s_rows.push((pair.sig, s));
    }

    if s_rows.is_empty() {
        return Ok(());
    }
    charge.hold(charge.bytes(s_rows.len(), s_terms, 0))?;

    let mut spoly_monomials: FxHashSet<Monomial> = FxHashSet::default();
    let mut since_check = 0usize;
    for (_, poly) in &s_rows {
        for term in &poly.terms {
            spoly_monomials.insert(term.mono.clone());
        }
        since_check += 1;
        if since_check >= ROW_CHECK_STRIDE {
            since_check = 0;
            charge.check(charge.bytes(0, 0, spoly_monomials.len()))?;
        }
    }
    charge.check(charge.bytes(0, 0, spoly_monomials.len()))?;
    let (reducer_rows, reducer_terms) =
        build_reducer_rows(basis, spoly_monomials, deadline, &charge)?;
    charge.hold(charge.bytes(reducer_rows.len(), reducer_terms, 0))?;

    let mut all_rows: Vec<(Signature, Polynomial)> =
        Vec::with_capacity(s_rows.len() + reducer_rows.len());
    all_rows.extend(s_rows);
    all_rows.extend(reducer_rows);
    let all_terms = s_terms + reducer_terms;

    let row_refs: Vec<(&Signature, &Polynomial)> =
        all_rows.iter().map(|(sig, poly)| (sig, poly)).collect();
    let mut new_indices: Vec<usize> = Vec::new();
    {
        let (columns, col_map) = collect_columns(&row_refs, &charge)?;

        // The sparse matrix is built from counts already known, so the
        // limit is checked before it allocates.
        if let Some(limit) = max_memory_bytes {
            if batch_bytes(all_rows.len(), all_terms, columns.len(), nvars) > limit {
                return Err(ComputeError::MemoryLimitExceeded);
            }
        }

        let matrix_rows: Vec<MatrixRow> = all_rows
            .into_iter()
            .map(|(sig, poly)| {
                let (cols, coeffs) = poly_to_row(&poly, &col_map);
                MatrixRow { sig, cols, coeffs }
            })
            .collect();

        let (pivot_rows, zero_sigs) = eliminate_rows(
            matrix_rows,
            modulus,
            columns.len(),
            deadline,
            max_memory_bytes,
        )?;
        for sig in zero_sigs {
            add_syzygy_rule(syzygy_rules, sig);
        }

        for row in pivot_rows {
            if row.is_zero() {
                continue;
            }
            let poly = row_to_poly(ring, &row.cols, &row.coeffs, &columns);
            if let Some(index) = insert_candidate(
                row.sig,
                poly,
                basis,
                syzygy_rules,
                rewriters,
                modulus,
                deadline,
            )? {
                new_indices.push(index);
            }
        }
    }

    for &idx in &new_indices {
        for i in 0..idx {
            pair_manager.push_pair(basis, i, idx)?;
        }
    }
    Ok(())
}

fn f5_reduce(
    mut p: LabeledPoly,
    basis: &Basis,
    modulus: u64,
    deadline: Option<Instant>,
) -> Result<LabeledPoly, ComputeError> {
    // Every step lowers the leading term, so the terms arrive in descending
    // order. One reverse at the end restores the ascending order a
    // polynomial holds.
    let mut remainder: Vec<Term> = Vec::new();
    let mut since_check = 0usize;

    while let Some(lt_p) = p.poly.lt() {
        since_check += p.poly.terms.len();
        if since_check >= REDUCE_DEADLINE_STRIDE {
            since_check = 0;
            poll_deadline(deadline)?;
        }
        let lead_mono = &lt_p.mono;
        let lead_coeff = lt_p.coeff;
        let lead_key = lead_mono.divisor_key();

        let mut reducer: Option<usize> = None;
        for filter in basis.filters() {
            if filter.deg > lead_mono.deg || !key_divides(filter.key, lead_key) {
                continue;
            }
            let index = filter.basis_pos as usize;
            let g = &basis[index];
            // a filter is only recorded for an element that has a lead.
            let lt_g = g.poly.lt().expect("a filtered element has a lead");
            if !lt_g.mono.divides(lead_mono) {
                continue;
            }
            if g.sig.shifted_is_below(lead_mono, &lt_g.mono, &p.sig) {
                reducer = Some(index);
                break;
            }
        }

        match reducer {
            Some(index) => {
                let g = &basis[index];
                // the scan only picks an element that has a lead.
                let lt_g = g.poly.lt().expect("the reducer has a lead");
                let m = lead_mono
                    .quotient(&lt_g.mono)
                    // divides() implies a quotient exists.
                    .expect("divides() implies quotient()");
                let scale = lead_coeff.div(lt_g.coeff, modulus);
                p.poly = p.poly.sub_scaled(&g.poly, scale, &m, modulus)?;
            }
            None => {
                let term = p
                    .poly
                    .pop_lt()
                    // lt was Some, so the polynomial is not empty here.
                    .expect("polynomial should not be empty");
                remainder.push(term);
            }
        }
    }

    remainder.reverse();
    p.poly = Polynomial::from_sorted_terms(p.poly.ring().clone(), remainder);
    Ok(p)
}

fn solve_raw_checked(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
) -> Result<Vec<Polynomial>, ComputeError> {
    super::check_input_degrees(generators)?;
    let budget = RunBudget {
        deadline,
        max_memory_bytes,
    };
    let modulus = ring.modulus();
    let nvars = ring.nvars();
    let mut basis = Basis::default();
    let mut syzygy_rules: Vec<Signature> = Vec::new();
    let mut pair_manager: PairManager = PairManager::default();
    let mut rewriters: Vec<Vec<RewriteRule>> = vec![Vec::new(); generators.len()];

    for (k, generator) in generators.iter().enumerate() {
        // The rewriter table, the basis, and the syzygy rules exist before
        // the first pair, so a run whose pair loop never starts must still
        // pay for them. Without this check a limit of zero bytes passes.
        check_limits(
            budget,
            &basis,
            &syzygy_rules,
            &pair_manager,
            &rewriters,
            0,
            nvars,
        )?;
        let sig = Signature {
            index: k,
            term: Monomial::one(nvars),
        };
        let labeled = LabeledPoly {
            sig,
            poly: generator.make_monic(modulus),
            index: basis.len(),
        };

        let reduced = f5_reduce(labeled, &basis, modulus, deadline)?;
        if reduced.poly.is_zero() {
            add_syzygy_rule(&mut syzygy_rules, reduced.sig);
            continue;
        }
        let mut labeled = reduced;

        let idx = basis.len();
        labeled.index = idx;
        add_rewriter(&mut rewriters, &labeled.sig, labeled.index);
        basis.push(labeled);

        for i in 0..idx {
            pair_manager.push_pair(&basis, i, idx)?;
        }

        while let Some((_degree, pairs)) = pair_manager.pop_smallest_degree() {
            check_limits(
                budget,
                &basis,
                &syzygy_rules,
                &pair_manager,
                &rewriters,
                pairs.len(),
                nvars,
            )?;
            process_degree(
                ring,
                pairs,
                &mut basis,
                &mut syzygy_rules,
                &mut rewriters,
                &mut pair_manager,
                budget,
            )?;
        }
    }

    check_limits(
        budget,
        &basis,
        &syzygy_rules,
        &pair_manager,
        &rewriters,
        0,
        nvars,
    )?;

    Ok(basis.into_polynomials())
}

/// Run the matrix backend under a deadline and a memory cap.
pub(super) fn solve_checked(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
) -> Result<Vec<Polynomial>, ComputeError> {
    let raw = solve_raw_checked(ring, generators, deadline, max_memory_bytes)?;
    reduced_groebner_basis_checked(ring, raw, deadline, max_memory_bytes)
}

/// Run the matrix backend with no budget.
#[cfg(test)]
pub(crate) fn solve(ring: &PolynomialRing, generators: &[Polynomial]) -> Vec<Polynomial> {
    // unlimited runs should only fail on internal invariants.
    solve_checked(ring, generators, None, None).expect("a run without a budget cannot stop early")
}

#[cfg(test)]
mod tests {
    use super::{solve, solve_checked};
    use crate::compute::ComputeError;
    use crate::compute::interreduce::is_reduced_basis;
    use crate::poly::Polynomial;
    use crate::ring::PolynomialRing;

    #[test]
    fn a_pair_past_the_degree_limit_is_a_typed_error() {
        let ring = PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime");
        let f = ring
            .polynomial([(1, [65535, 0]), (1, [0, 65535])])
            .expect("fits");
        let g = ring.polynomial([(1, [0, 65535])]).expect("fits");
        assert_eq!(
            solve_checked(&ring, &[f, g], None, None),
            Err(ComputeError::DegreeLimit { limit: 65535 })
        );
    }

    fn system(ring: &PolynomialRing, texts: &[&str]) -> Vec<Polynomial> {
        texts
            .iter()
            .map(|text| ring.parse_polynomial(text).expect("the text parses"))
            .collect()
    }

    #[test]
    fn the_matrix_backend_reduces_a_monomial_ideal() {
        let ring = PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime");
        let gb = solve(&ring, &system(&ring, &["x^2", "x*y"]));
        assert!(is_reduced_basis(&gb, 7));
    }

    #[test]
    fn the_matrix_backend_reduces_cyclic_3() {
        let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"]).expect("32003 is prime");
        let gb = solve(
            &ring,
            &system(&ring, &["x + y + z", "x*y + y*z + z*x", "x*y*z - 1"]),
        );
        assert!(is_reduced_basis(&gb, 32003));
    }

    #[test]
    fn the_matrix_backend_reduces_cyclic_4() {
        let ring =
            PolynomialRing::prime_field(32003, ["x", "y", "z", "w"]).expect("32003 is prime");
        let gb = solve(
            &ring,
            &system(
                &ring,
                &[
                    "x + y + z + w",
                    "x*y + y*z + z*w + w*x",
                    "x*y*z + y*z*w + z*w*x + w*x*y",
                    "x*y*z*w - 1",
                ],
            ),
        );
        assert!(is_reduced_basis(&gb, 32003));
    }

    #[test]
    fn the_matrix_backend_reduces_katsura_4() {
        let ring =
            PolynomialRing::prime_field(32003, ["a", "b", "c", "d", "e"]).expect("32003 is prime");
        let gb = solve(
            &ring,
            &system(
                &ring,
                &[
                    "a + 2*b + 2*c + 2*d + 2*e - 1",
                    "2*a*b - b",
                    "2*a*c + b^2 - c",
                    "2*a*d + 2*b*c - d",
                    "2*a*e + 2*b*d + c^2 - e",
                ],
            ),
        );
        assert!(is_reduced_basis(&gb, 32003));
    }

    /// The katsura-`n` system over `n + 1` variables `v_0..v_n`.
    fn katsura(ring: &PolynomialRing, n: usize) -> Vec<Polynomial> {
        let mut generators = Vec::with_capacity(n + 1);

        let mut linear: Vec<(i64, Vec<u16>)> = Vec::with_capacity(n + 1);
        let mut lead = vec![0u16; n + 1];
        lead[0] = 1;
        linear.push((1, lead));
        for i in 1..=n {
            let mut exps = vec![0u16; n + 1];
            exps[i] = 1;
            linear.push((2, exps));
        }
        linear.push((-1, vec![0u16; n + 1]));
        generators.push(ring.polynomial(linear).expect("the exponents fit the ring"));

        for k in 1..=n {
            let mut terms: Vec<(i64, Vec<u16>)> = Vec::new();
            for i in 0..=k / 2 {
                let j = k - i;
                let coeff = if i == j { 1 } else { 2 };
                let mut exps = vec![0u16; n + 1];
                exps[i] += 1;
                exps[j] += 1;
                terms.push((coeff, exps));
            }
            let mut lead = vec![0u16; n + 1];
            lead[k] = 1;
            terms.push((-1, lead));
            generators.push(ring.polynomial(terms).expect("the exponents fit the ring"));
        }

        generators
    }

    fn katsura_6_ring() -> PolynomialRing {
        PolynomialRing::prime_field(32003, ["a", "b", "c", "d", "e", "f", "g"])
            .expect("32003 is prime")
    }

    #[test]
    fn a_tight_memory_limit_stops_matrix_elimination_inside_a_batch() {
        // The limit passes every outer check between degree batches, where
        // the basis and the pair queue are still small: the whole run stays
        // under 46,000 bytes of that accounting. But katsura-6 builds a
        // matrix inside one batch that holds close to 80,000 bytes at its
        // largest, so a limit between the two only stops the run if the
        // batch itself is accounted for.
        const LIMIT: usize = 50_000;
        let ring = katsura_6_ring();
        let gens = katsura(&ring, 6);
        solve(&ring, &gens);
        assert_eq!(
            solve_checked(&ring, &gens, None, Some(LIMIT)),
            Err(ComputeError::MemoryLimitExceeded)
        );
    }

    #[test]
    fn a_tight_timeout_stops_matrix_elimination_inside_a_batch() {
        let ring = katsura_6_ring();
        let gens = katsura(&ring, 6);

        let start = std::time::Instant::now();
        solve(&ring, &gens);
        let full = start.elapsed();

        let start = std::time::Instant::now();
        let error = solve_checked(&ring, &gens, Some(start + full / 20), None)
            .expect_err("the deadline stops the run");
        let elapsed = start.elapsed();

        assert_eq!(error, ComputeError::Timeout);
        assert!(
            elapsed < full,
            "the run must stop near the deadline: {elapsed:?} of a full {full:?}"
        );
    }
}
