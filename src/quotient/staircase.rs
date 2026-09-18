//! Standard monomials of a finite monomial quotient.

use std::collections::{HashMap, HashSet, VecDeque};
use std::mem::size_of;

use crate::compute::{ComputeError, ComputeLimits, DEGREE_LIMIT, RunError};
use crate::ideal::GroebnerBasis;
use crate::poly::{Monomial, Polynomial, heap_exps_bytes};
use crate::ring::Domain;

use super::QuotientError;

/// The standard monomials and their coordinate lookup.
#[derive(Clone, Debug)]
pub(crate) struct Staircase {
    /// Monomials in ascending grevlex order.
    pub(crate) monomials: Vec<Monomial>,
    /// Exponents in the same order as `monomials`.
    pub(crate) basis: Vec<Vec<u16>>,
    /// Index by exponent vector.
    pub(crate) index: HashMap<Vec<u16>, usize>,
}

/// Estimate the heap bytes retained by a completed staircase.
///
/// The estimate includes element slots and vector capacities. It excludes
/// allocator metadata and hash-table control storage.
pub(crate) fn retained_staircase_bytes(staircase: &Staircase, nvars: usize) -> usize {
    let monomials = staircase
        .monomials
        .capacity()
        .saturating_mul(size_of::<Monomial>())
        .saturating_add(
            staircase
                .monomials
                .len()
                .saturating_mul(heap_exps_bytes(nvars)),
        );
    let basis = basis_storage_bytes(&staircase.basis, staircase.basis.capacity());
    let index_keys = index_storage_bytes(&staircase.index);
    monomials.saturating_add(basis).saturating_add(index_keys)
}

/// Build the finite staircase of `source` under `limits`.
///
/// The walk starts at 1 and extends a standard monomial by one variable at
/// a time. Upward closure of a monomial ideal makes this frontier complete.
pub(crate) fn build<D: Domain>(
    source: &GroebnerBasis<D>,
    limits: &ComputeLimits,
) -> Result<Staircase, QuotientError> {
    let nvars = source.ring().nvars();
    let leading = leading_monomials(source, limits)?;
    if !zero_dimensional_from_leading(&leading, nvars) {
        return Err(QuotientError::NotFinite);
    }
    StaircaseWalk::new(source, leading, limits)?.complete()
}

struct StaircaseWalk<'a> {
    limits: &'a ComputeLimits,
    nvars: usize,
    retained: usize,
    leading: Vec<Monomial>,
    leading_bytes: usize,
    seen: HashSet<Monomial>,
    frontier: VecDeque<Monomial>,
    standard: Vec<Monomial>,
    seen_bytes: usize,
    frontier_bytes: usize,
    standard_bytes: usize,
}

impl<'a> StaircaseWalk<'a> {
    fn new<D: Domain>(
        source: &GroebnerBasis<D>,
        leading: Vec<Monomial>,
        limits: &'a ComputeLimits,
    ) -> Result<Self, QuotientError> {
        let nvars = source.ring().nvars();
        let retained = retained_basis_bytes(source);
        let leading_bytes = leading_bytes(&leading, leading.capacity(), nvars);
        let one = Monomial::one(nvars);
        let one_bytes = monomial_bytes(&one, nvars);
        check(
            limits,
            retained
                .saturating_add(leading_bytes)
                .saturating_add(one_bytes.saturating_mul(2))
                .saturating_add(size_of::<Staircase>()),
        )?;
        let mut seen = HashSet::new();
        seen.insert(one.clone());
        let mut frontier = VecDeque::new();
        frontier.push_back(one);
        Ok(StaircaseWalk {
            limits,
            nvars,
            retained,
            leading,
            leading_bytes,
            seen,
            frontier,
            standard: Vec::new(),
            seen_bytes: one_bytes,
            frontier_bytes: one_bytes,
            standard_bytes: 0,
        })
    }

    fn complete(mut self) -> Result<Staircase, QuotientError> {
        while let Some(monomial) = self.frontier.pop_front() {
            self.frontier_bytes = self
                .frontier_bytes
                .saturating_sub(monomial_bytes(&monomial, self.nvars));
            self.process(monomial)?;
        }
        self.finish()
    }

    fn process(&mut self, monomial: Monomial) -> Result<(), QuotientError> {
        self.check(0)?;
        if self.leading.iter().any(|lead| lead.divides(&monomial)) {
            return Ok(());
        }
        let standard_extra = monomial_bytes(&monomial, self.nvars);
        self.check(standard_extra.saturating_add(size_of::<Monomial>()))?;
        self.standard.push(monomial.clone());
        self.standard_bytes = self.standard_bytes.saturating_add(standard_extra);
        self.enqueue_children(&monomial)
    }

    fn enqueue_children(&mut self, monomial: &Monomial) -> Result<(), QuotientError> {
        for variable in 0..self.nvars {
            self.check(0)?;
            self.enqueue_child(monomial, variable)?;
        }
        Ok(())
    }

    fn enqueue_child(&mut self, monomial: &Monomial, variable: usize) -> Result<(), QuotientError> {
        let mut child = monomial.clone();
        let exponent = child.exps[variable]
            .checked_add(1)
            .ok_or(QuotientError::ExponentLimit {
                limit: DEGREE_LIMIT,
            })?;
        child.exps[variable] = exponent;
        child.deg = child
            .deg
            .checked_add(1)
            .ok_or(QuotientError::ExponentLimit {
                limit: DEGREE_LIMIT,
            })?;
        if self.seen.contains(&child) {
            return Ok(());
        }

        let child_bytes = monomial_bytes(&child, self.nvars);
        self.check(
            child_bytes
                .saturating_mul(4)
                .saturating_add(size_of::<Monomial>()),
        )?;
        self.seen.insert(child.clone());
        self.frontier.push_back(child);
        self.seen_bytes = self.seen_bytes.saturating_add(child_bytes);
        self.frontier_bytes = self.frontier_bytes.saturating_add(child_bytes);
        Ok(())
    }

    fn finish(mut self) -> Result<Staircase, QuotientError> {
        self.standard.sort_unstable();
        let basis_payload = exponent_payload(&self.standard, size_of::<u16>());
        let basis_slots = self.standard.len().saturating_mul(size_of::<Vec<u16>>());
        self.check(basis_payload.saturating_add(basis_slots))?;
        let mut basis = Vec::with_capacity(self.standard.len());
        basis.extend(self.standard.iter().map(|monomial| monomial.exps.to_vec()));
        let index_key_payload = basis
            .iter()
            .map(|exponents| exponents.capacity().saturating_mul(size_of::<u16>()))
            .fold(0usize, usize::saturating_add);
        let index_slots = basis.len().saturating_mul(size_of::<(Vec<u16>, usize)>());
        self.check(
            basis_storage_bytes(&basis, basis.capacity())
                .saturating_add(index_key_payload)
                .saturating_add(index_slots),
        )?;
        let mut index = HashMap::with_capacity(basis.len());
        for (position, exps) in basis.iter().enumerate() {
            index.insert(exps.clone(), position);
        }
        self.check(
            basis_storage_bytes(&basis, basis.capacity())
                .saturating_add(index_storage_bytes(&index)),
        )?;
        Ok(Staircase {
            monomials: self.standard,
            basis,
            index,
        })
    }

    fn check(&self, extra: usize) -> Result<(), QuotientError> {
        check(
            self.limits,
            self.retained
                .saturating_add(self.leading_bytes)
                .saturating_add(self.seen_bytes)
                .saturating_add(self.seen_spare_bytes())
                .saturating_add(self.frontier_bytes)
                .saturating_add(self.frontier_spare_bytes())
                .saturating_add(self.standard_bytes)
                .saturating_add(self.standard_spare_bytes())
                .saturating_add(extra),
        )
    }

    fn seen_spare_bytes(&self) -> usize {
        self.seen
            .capacity()
            .saturating_sub(self.seen.len())
            .saturating_mul(size_of::<Monomial>())
    }

    fn frontier_spare_bytes(&self) -> usize {
        self.frontier
            .capacity()
            .saturating_sub(self.frontier.len())
            .saturating_mul(size_of::<Monomial>())
    }

    fn standard_spare_bytes(&self) -> usize {
        self.standard
            .capacity()
            .saturating_sub(self.standard.len())
            .saturating_mul(size_of::<Monomial>())
    }
}

/// Report whether the leading monomial ideal is zero-dimensional.
pub(crate) fn zero_dimensional<D: Domain>(source: &GroebnerBasis<D>) -> bool {
    let nvars = source.ring().nvars();
    let leading: Vec<Monomial> = source
        .iter()
        .filter_map(|poly| polynomial_monomial(poly))
        .collect();
    zero_dimensional_from_leading(&leading, nvars)
}

fn zero_dimensional_from_leading(leading: &[Monomial], nvars: usize) -> bool {
    if leading.iter().any(|m| m.exps.iter().all(|&exp| exp == 0)) {
        return true;
    }
    if nvars == 0 {
        return true;
    }
    (0..nvars).all(|variable| {
        leading.iter().any(|m| {
            m.exps[variable] > 0
                && m.exps
                    .iter()
                    .enumerate()
                    .all(|(index, &exp)| index == variable || exp == 0)
        })
    })
}

fn leading_monomials<D: Domain>(
    source: &GroebnerBasis<D>,
    limits: &ComputeLimits,
) -> Result<Vec<Monomial>, QuotientError> {
    let nvars = source.ring().nvars();
    let count = source
        .iter()
        .filter(|poly| poly.leading_term().is_some())
        .count();
    check(
        limits,
        retained_basis_bytes(source)
            .saturating_add(count.saturating_mul(size_of::<Monomial>()))
            .saturating_add(count.saturating_mul(heap_exps_bytes(nvars))),
    )?;
    Ok(source
        .iter()
        .filter_map(|poly| polynomial_monomial(poly))
        .collect())
}

fn polynomial_monomial<D: Domain>(poly: &Polynomial<D>) -> Option<Monomial> {
    poly.leading_term()
        .map(|(_, exps)| Monomial::from_exps(exps.iter().copied().collect()))
}

/// Estimate bytes retained by a basis clone and its polynomial data.
///
/// The basis vector is charged at its element count because `GroebnerBasis`
/// exposes its vector only as a slice. Polynomial capacities and coefficients
/// are charged exactly by their retained estimates.
pub(crate) fn retained_basis_bytes<D: Domain>(source: &GroebnerBasis<D>) -> usize {
    size_of::<GroebnerBasis<D>>()
        .saturating_add(source.len().saturating_mul(size_of::<Polynomial<D>>()))
        .saturating_add(
            source
                .iter()
                .map(Polynomial::retained_bytes)
                .fold(0usize, usize::saturating_add),
        )
}

fn leading_bytes(leading: &[Monomial], capacity: usize, nvars: usize) -> usize {
    leading
        .len()
        .saturating_mul(heap_exps_bytes(nvars))
        .saturating_add(capacity.saturating_mul(size_of::<Monomial>()))
}

fn exponent_payload<T>(values: &[T], item_bytes: usize) -> usize
where
    T: ExponentLength,
{
    values
        .iter()
        .map(|value| value.exponent_len().saturating_mul(item_bytes))
        .fold(0usize, usize::saturating_add)
}

trait ExponentLength {
    fn exponent_len(&self) -> usize;
}

impl ExponentLength for Monomial {
    fn exponent_len(&self) -> usize {
        self.exps.len()
    }
}

impl ExponentLength for Vec<u16> {
    fn exponent_len(&self) -> usize {
        self.capacity()
    }
}

fn basis_storage_bytes(basis: &[Vec<u16>], capacity: usize) -> usize {
    capacity
        .saturating_mul(size_of::<Vec<u16>>())
        .saturating_add(
            basis
                .iter()
                .map(|exponents| exponents.capacity().saturating_mul(size_of::<u16>()))
                .fold(0usize, usize::saturating_add),
        )
}

fn index_storage_bytes(index: &HashMap<Vec<u16>, usize>) -> usize {
    index
        .keys()
        .map(|exponents| exponents.capacity().saturating_mul(size_of::<u16>()))
        .fold(0usize, usize::saturating_add)
        .saturating_add(
            index
                .capacity()
                .saturating_mul(size_of::<(Vec<u16>, usize)>()),
        )
}

fn monomial_bytes(monomial: &Monomial, nvars: usize) -> usize {
    size_of::<Monomial>().saturating_add(if monomial.exps.len() > 11 {
        heap_exps_bytes(nvars)
    } else {
        0
    })
}

fn check(limits: &ComputeLimits, bytes: usize) -> Result<(), QuotientError> {
    if let Some(stop) = limits.stop() {
        return Err(from_run_error(stop));
    }
    if limits.memory.is_some_and(|limit| bytes > limit) {
        return Err(QuotientError::MemoryLimitExceeded);
    }
    Ok(())
}

fn from_run_error(error: RunError) -> QuotientError {
    match error {
        RunError::Compute(ComputeError::MemoryLimitExceeded) => QuotientError::MemoryLimitExceeded,
        RunError::Compute(ComputeError::Timeout) | RunError::Cancelled => QuotientError::Timeout,
        RunError::Compute(ComputeError::ExponentLimit { limit }) => {
            QuotientError::ExponentLimit { limit }
        }
        RunError::Compute(
            ComputeError::DegreeLimit { .. }
            | ComputeError::TableFull
            | ComputeError::PrimesExhausted,
        ) => unreachable!("staircase budget checks do not report this stop"),
    }
}
