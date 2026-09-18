//! Budgeted arithmetic for polynomials.

use std::any::TypeId;
use std::cmp::Ordering;
use std::mem::size_of;

use crate::compute::{Budget, ComputeError, ComputeLimits, DEGREE_LIMIT, RunError};
use crate::poly::{ExponentOverflow, Monomial, Polynomial, Term, heap_exps_bytes};
use crate::ring::rational::limb_bytes;
use crate::ring::{Coefficient, Domain, DomainOps, Rationals, RingError};

const RATIONAL_WORKSPACE_FACTOR: usize = 6;
const RATIONAL_LIMBS_PER_OPERATION: usize = 8;

/// Why a polynomial arithmetic operation stopped.
///
/// The operation checks its ring and its budget before it allocates a result.
/// A coefficient conversion error is reported only by [`Polynomial::try_scale`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ArithmeticError {
    /// The two operands belong to different rings.
    RingMismatch,
    /// A supplied coefficient cannot be read by the polynomial's domain.
    CoefficientConversion(RingError),
    /// A product needs an exponent above `limit`.
    ExponentLimit {
        /// The largest value one exponent can hold.
        limit: u32,
    },
    /// The deadline passed.
    Timeout,
    /// The operation's live data passed the memory limit.
    MemoryLimitExceeded,
}

impl std::fmt::Display for ArithmeticError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ArithmeticError::RingMismatch => {
                f.write_str("the polynomial operands belong to different rings")
            }
            ArithmeticError::CoefficientConversion(error) => {
                write!(f, "coefficient conversion failed: {error}")
            }
            ArithmeticError::ExponentLimit { limit } => {
                write!(f, "a product reached an exponent above {limit}")
            }
            ArithmeticError::Timeout => {
                f.write_str("the polynomial arithmetic passed its deadline")
            }
            ArithmeticError::MemoryLimitExceeded => {
                f.write_str("the polynomial arithmetic passed its memory limit")
            }
        }
    }
}

impl std::error::Error for ArithmeticError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ArithmeticError::CoefficientConversion(error) => Some(error),
            _ => None,
        }
    }
}

impl From<ExponentOverflow> for ArithmeticError {
    fn from(_: ExponentOverflow) -> Self {
        ArithmeticError::ExponentLimit {
            limit: DEGREE_LIMIT,
        }
    }
}

impl<D: Domain> Polynomial<D> {
    /// Add two polynomials under `budget`.
    ///
    /// The operands must belong to equal rings. The result preserves the
    /// ring's ascending term order. The memory limit charges both operands
    /// and the result capacity.
    pub fn try_add(&self, other: &Self, budget: Budget) -> Result<Self, ArithmeticError> {
        let limits = ComputeLimits::of_budget(&budget);
        self.try_add_with_limits(other, &limits)
    }

    /// Subtract `other` from this polynomial under `budget`.
    ///
    /// The operands must belong to equal rings. The memory limit charges both
    /// operands and the result capacity.
    pub fn try_sub(&self, other: &Self, budget: Budget) -> Result<Self, ArithmeticError> {
        let limits = ComputeLimits::of_budget(&budget);
        self.try_sub_with_limits(other, &limits)
    }

    /// Multiply two polynomials under `budget`.
    ///
    /// The operands must belong to equal rings. A product that needs an
    /// exponent above `u16::MAX` returns [`ArithmeticError::ExponentLimit`].
    pub fn try_mul(&self, other: &Self, budget: Budget) -> Result<Self, ArithmeticError> {
        let limits = ComputeLimits::of_budget(&budget);
        self.try_mul_with_limits(other, &limits)
    }

    /// Negate a polynomial under `budget`.
    ///
    /// The memory limit charges the input and the result capacity.
    pub fn try_neg(&self, budget: Budget) -> Result<Self, ArithmeticError> {
        let limits = ComputeLimits::of_budget(&budget);
        self.try_neg_with_limits(&limits)
    }

    /// Raise a polynomial to a nonnegative integer power under `budget`.
    ///
    /// The exponent zero returns the constant polynomial one. Every product
    /// uses the same deadline and memory limit as this call.
    pub fn try_pow(&self, exponent: u32, budget: Budget) -> Result<Self, ArithmeticError> {
        let limits = ComputeLimits::of_budget(&budget);
        self.try_pow_with_limits(exponent, &limits)
    }

    /// Scale a polynomial by a coefficient under `budget`.
    ///
    /// The coefficient is converted through the polynomial's domain. A
    /// failed conversion returns [`ArithmeticError::CoefficientConversion`].
    pub fn try_scale(
        &self,
        value: impl Into<Coefficient>,
        budget: Budget,
    ) -> Result<Self, ArithmeticError> {
        let limits = ComputeLimits::of_budget(&budget);
        self.try_scale_with_limits(value, &limits)
    }

    /// Add two polynomials under absolute crate limits.
    pub(crate) fn try_add_with_limits(
        &self,
        other: &Self,
        limits: &ComputeLimits,
    ) -> Result<Self, ArithmeticError> {
        check_same_ring(self, other)?;
        check_stop(limits)?;
        let capacity = checked_sum(self.terms.len(), other.terms.len())?;
        let output = term_storage::<D>(self.ring().nvars(), capacity)?;
        let coefficients = sum_output_coefficients(self, other)?;
        let left_held = retained_bytes(self)?;
        let right_held = retained_bytes(other)?;
        check_memory(limits, &[left_held, right_held, output, coefficients])?;
        let terms = combine_terms(self, other, false, limits)?;
        finish(self, terms, limits)
    }

    /// Subtract two polynomials under absolute crate limits.
    pub(crate) fn try_sub_with_limits(
        &self,
        other: &Self,
        limits: &ComputeLimits,
    ) -> Result<Self, ArithmeticError> {
        check_same_ring(self, other)?;
        check_stop(limits)?;
        let capacity = checked_sum(self.terms.len(), other.terms.len())?;
        let output = term_storage::<D>(self.ring().nvars(), capacity)?;
        let coefficients = sum_output_coefficients(self, other)?;
        let left_held = retained_bytes(self)?;
        let right_held = retained_bytes(other)?;
        check_memory(limits, &[left_held, right_held, output, coefficients])?;
        let terms = combine_terms(self, other, true, limits)?;
        finish(self, terms, limits)
    }

    /// Multiply two polynomials under absolute crate limits.
    pub(crate) fn try_mul_with_limits(
        &self,
        other: &Self,
        limits: &ComputeLimits,
    ) -> Result<Self, ArithmeticError> {
        check_same_ring(self, other)?;
        mul_with_extra(self, other, limits, 0)
    }

    /// Negate a polynomial under absolute crate limits.
    pub(crate) fn try_neg_with_limits(
        &self,
        limits: &ComputeLimits,
    ) -> Result<Self, ArithmeticError> {
        check_stop(limits)?;
        let capacity = self.terms.len();
        let output = term_storage::<D>(self.ring().nvars(), capacity)?;
        let held = retained_bytes(self)?;
        check_memory(limits, &[held, output, self.coefficient_bytes()])?;
        let mut terms = Vec::with_capacity(capacity);
        for (index, term) in self.terms.iter().enumerate() {
            check_work(limits, index)?;
            terms.push(Term {
                coeff: self.ring().ops().neg(&term.coeff),
                mono: term.mono.clone(),
            });
        }
        finish(self, terms, limits)
    }

    /// Raise a polynomial under absolute crate limits.
    pub(crate) fn try_pow_with_limits(
        &self,
        exponent: u32,
        limits: &ComputeLimits,
    ) -> Result<Self, ArithmeticError> {
        check_stop(limits)?;
        if exponent == 0 {
            return one_with_limits(self, limits);
        }
        if self.is_zero() {
            check_memory(limits, &[retained_bytes(self)?])?;
            return Ok(self.zero_like());
        }
        if exponent == 1 {
            let held = retained_bytes(self)?;
            check_memory(limits, &[held, held])?;
            return Ok(self.clone());
        }
        check_power_exponents(self, exponent, limits)?;
        pow_nontrivial(self, exponent, limits)
    }

    /// Scale a polynomial under absolute crate limits.
    pub(crate) fn try_scale_with_limits(
        &self,
        value: impl Into<Coefficient>,
        limits: &ComputeLimits,
    ) -> Result<Self, ArithmeticError> {
        check_stop(limits)?;
        check_memory(limits, &[retained_bytes(self)?])?;
        let value = value.into();
        let Some(scalar) = converted_scalar(self, &value, limits)? else {
            return Ok(self.zero_like());
        };
        scale_nonzero(self, &scalar, limits)
    }
}

fn check_same_ring<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
) -> Result<(), ArithmeticError> {
    if left.ring() == right.ring() {
        Ok(())
    } else {
        Err(ArithmeticError::RingMismatch)
    }
}

fn checked_sum(left: usize, right: usize) -> Result<usize, ArithmeticError> {
    left.checked_add(right)
        .ok_or(ArithmeticError::MemoryLimitExceeded)
}

fn term_storage<D: Domain>(nvars: usize, terms: usize) -> Result<usize, ArithmeticError> {
    let per_term = size_of::<Term<D>>()
        .checked_add(heap_exps_bytes(nvars))
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    per_term
        .checked_mul(terms)
        .ok_or(ArithmeticError::MemoryLimitExceeded)
}

fn retained_bytes<D: Domain>(poly: &Polynomial<D>) -> Result<usize, ArithmeticError> {
    let vector = size_of::<Term<D>>()
        .checked_mul(poly.terms.capacity())
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    let exponents = heap_exps_bytes(poly.ring().nvars())
        .checked_mul(poly.terms.len())
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    vector
        .checked_add(exponents)
        .and_then(|bytes| bytes.checked_add(poly.coefficient_bytes()))
        .ok_or(ArithmeticError::MemoryLimitExceeded)
}

fn rational_workspace<D: Domain>(
    input_bytes: usize,
    operations: usize,
) -> Result<usize, ArithmeticError> {
    if TypeId::of::<D>() != TypeId::of::<Rationals>() || input_bytes == 0 || operations == 0 {
        return Ok(0);
    }
    // Rational operations retain cross products, an lcm, and reduction data
    // while the result is built. The factor leaves room for those values.
    let scaled = input_bytes
        .checked_mul(RATIONAL_WORKSPACE_FACTOR)
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    let limb_overhead = operations
        .checked_mul(RATIONAL_LIMBS_PER_OPERATION * u64::BITS as usize)
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    scaled
        .checked_add(limb_overhead)
        .ok_or(ArithmeticError::MemoryLimitExceeded)
}

fn sum_output_coefficients<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
) -> Result<usize, ArithmeticError> {
    let overlap = left.terms.len().min(right.terms.len());
    let bytes = left
        .coefficient_bytes()
        .checked_add(right.coefficient_bytes())
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    bytes
        .checked_add(rational_workspace::<D>(bytes, overlap)?)
        .ok_or(ArithmeticError::MemoryLimitExceeded)
}

fn scaled_output_coefficients<D: Domain>(
    poly: &Polynomial<D>,
    scalar: usize,
) -> Result<usize, ArithmeticError> {
    let terms = poly.terms.len();
    let bytes = poly
        .coefficient_bytes()
        .checked_add(
            scalar
                .checked_mul(terms)
                .ok_or(ArithmeticError::MemoryLimitExceeded)?,
        )
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    bytes
        .checked_add(rational_workspace::<D>(bytes, terms)?)
        .ok_or(ArithmeticError::MemoryLimitExceeded)
}

fn scalar_heap_bound(value_bytes: usize, domain_floor: usize) -> usize {
    if domain_floor == 0 {
        0
    } else {
        value_bytes.max(domain_floor)
    }
}

fn converted_scalar<D: Domain>(
    poly: &Polynomial<D>,
    value: &Coefficient,
    limits: &ComputeLimits,
) -> Result<Option<D::Coeff>, ArithmeticError> {
    check_stop(limits)?;
    let held = retained_bytes(poly)?;
    check_memory(limits, &[held])?;
    let value_bytes = coefficient_input_bytes(value);
    let scalar_bound = scalar_heap_bound(value_bytes, domain_heap_floor::<D>());
    let output = term_storage::<D>(poly.ring().nvars(), poly.terms.len())?;
    let coefficients = scaled_output_coefficients(poly, scalar_bound)?;
    check_memory(
        limits,
        &[held, value_bytes, scalar_bound, output, coefficients],
    )?;
    let scalar = poly
        .ring()
        .ops()
        .convert(value)
        .map_err(ArithmeticError::CoefficientConversion)?;
    check_stop(limits)?;
    Ok(scalar)
}

fn scale_nonzero<D: Domain>(
    poly: &Polynomial<D>,
    scalar: &D::Coeff,
    limits: &ComputeLimits,
) -> Result<Polynomial<D>, ArithmeticError> {
    if poly.is_zero() {
        return Ok(poly.zero_like());
    }
    let ops = poly.ring().ops();
    let mut terms = Vec::with_capacity(poly.terms.len());
    for (index, term) in poly.terms.iter().enumerate() {
        check_work(limits, index)?;
        terms.push(Term {
            coeff: ops.mul(&term.coeff, scalar),
            mono: term.mono.clone(),
        });
    }
    finish(poly, terms, limits)
}

fn domain_heap_floor<D: Domain>() -> usize {
    if TypeId::of::<D>() == TypeId::of::<Rationals>() {
        2 * size_of::<u64>()
    } else {
        0
    }
}

fn product_output_coefficients<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
    pairs: usize,
) -> Result<usize, ArithmeticError> {
    let left_bytes = left
        .coefficient_bytes()
        .checked_mul(right.terms.len())
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    let right_bytes = right
        .coefficient_bytes()
        .checked_mul(left.terms.len())
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    let bytes = left_bytes
        .checked_add(right_bytes)
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    let carry = rational_workspace::<D>(bytes, pairs)?;
    bytes
        .checked_add(carry)
        .ok_or(ArithmeticError::MemoryLimitExceeded)
}

fn check_stop(limits: &ComputeLimits) -> Result<(), ArithmeticError> {
    if let Some(stop) = limits.stop() {
        return Err(stop_error(stop));
    }
    Ok(())
}

fn check_work(limits: &ComputeLimits, index: usize) -> Result<(), ArithmeticError> {
    if let Some(stop) = limits.stop_every(index) {
        return Err(stop_error(stop));
    }
    Ok(())
}

fn stop_error(stop: RunError) -> ArithmeticError {
    match stop.reported() {
        ComputeError::MemoryLimitExceeded => ArithmeticError::MemoryLimitExceeded,
        _ => ArithmeticError::Timeout,
    }
}

fn check_memory(limits: &ComputeLimits, parts: &[usize]) -> Result<(), ArithmeticError> {
    check_stop(limits)?;
    let mut bytes = 0usize;
    for &part in parts {
        bytes = bytes
            .checked_add(part)
            .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    }
    if limits.memory.is_some_and(|limit| bytes > limit) {
        return Err(ArithmeticError::MemoryLimitExceeded);
    }
    Ok(())
}

fn finish<D: Domain>(
    source: &Polynomial<D>,
    terms: Vec<Term<D>>,
    limits: &ComputeLimits,
) -> Result<Polynomial<D>, ArithmeticError> {
    let result = Polynomial::from_sorted_terms(source.ring().clone(), terms);
    check_memory(limits, &[retained_bytes(source)?, retained_bytes(&result)?])?;
    Ok(result)
}

fn combine_terms<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
    subtract: bool,
    limits: &ComputeLimits,
) -> Result<Vec<Term<D>>, ArithmeticError> {
    let capacity = checked_sum(left.terms.len(), right.terms.len())?;
    let mut out = Vec::with_capacity(capacity);
    let mut left_index = 0usize;
    let mut right_index = 0usize;
    let mut work = 0usize;
    let ops = left.ring().ops();
    while left_index < left.terms.len() && right_index < right.terms.len() {
        check_work(limits, work)?;
        let left_term = &left.terms[left_index];
        let right_term = &right.terms[right_index];
        match left_term.mono.cmp(&right_term.mono) {
            Ordering::Less => {
                out.push(left_term.clone());
                left_index += 1;
            }
            Ordering::Greater => {
                out.push(signed_term(right_term, subtract, ops));
                right_index += 1;
            }
            Ordering::Equal => {
                if let Some(coeff) = combined_coefficient(left_term, right_term, subtract, ops) {
                    out.push(Term {
                        coeff,
                        mono: left_term.mono.clone(),
                    });
                }
                left_index += 1;
                right_index += 1;
            }
        }
        work = work.saturating_add(1);
    }
    copy_tail(
        &mut out,
        &left.terms[left_index..],
        false,
        ops,
        limits,
        work,
    )?;
    let right_work = work.saturating_add(left.terms.len().saturating_sub(left_index));
    copy_tail(
        &mut out,
        &right.terms[right_index..],
        subtract,
        ops,
        limits,
        right_work,
    )?;
    Ok(out)
}

fn signed_term<D: Domain>(term: &Term<D>, negate: bool, ops: &D::Ops) -> Term<D> {
    let coeff = if negate {
        ops.neg(&term.coeff)
    } else {
        term.coeff.clone()
    };
    Term {
        coeff,
        mono: term.mono.clone(),
    }
}

fn combined_coefficient<D: Domain>(
    left: &Term<D>,
    right: &Term<D>,
    subtract: bool,
    ops: &D::Ops,
) -> Option<D::Coeff> {
    if subtract {
        ops.sub(&left.coeff, &right.coeff)
    } else {
        ops.add(&left.coeff, &right.coeff)
    }
}

fn copy_tail<D: Domain>(
    out: &mut Vec<Term<D>>,
    terms: &[Term<D>],
    negate: bool,
    ops: &D::Ops,
    limits: &ComputeLimits,
    mut work: usize,
) -> Result<(), ArithmeticError> {
    for term in terms {
        check_work(limits, work)?;
        out.push(signed_term(term, negate, ops));
        work = work.saturating_add(1);
    }
    Ok(())
}

fn mul_with_extra<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
    limits: &ComputeLimits,
    extra: usize,
) -> Result<Polynomial<D>, ArithmeticError> {
    let pairs = product_preflight(left, right, limits, extra)?;
    if pairs == 0 {
        return Ok(left.zero_like());
    }
    let terms = product_terms(left, right, limits)?;
    let terms = merge_product_terms(terms, left.ring().ops(), limits)?;
    finish(left, terms, limits)
}

fn product_preflight<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
    limits: &ComputeLimits,
    extra: usize,
) -> Result<usize, ArithmeticError> {
    check_stop(limits)?;
    let pairs = pair_count(left, right)?;
    if pairs == 0 {
        empty_product_preflight(left, right, limits, extra)?;
        return Ok(0);
    }
    let term_bytes = term_storage::<D>(left.ring().nvars(), pairs)?;
    let coefficient_bytes = product_output_coefficients(left, right, pairs)?;
    product_memory_preflight(left, right, limits, extra, term_bytes, coefficient_bytes)?;
    check_products(left, right, limits)?;
    Ok(pairs)
}

fn pair_count<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
) -> Result<usize, ArithmeticError> {
    left.terms
        .len()
        .checked_mul(right.terms.len())
        .ok_or(ArithmeticError::MemoryLimitExceeded)
}

fn empty_product_preflight<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
    limits: &ComputeLimits,
    extra: usize,
) -> Result<(), ArithmeticError> {
    check_memory(
        limits,
        &[retained_bytes(left)?, retained_bytes(right)?, extra],
    )
}

fn product_memory_preflight<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
    limits: &ComputeLimits,
    extra: usize,
    term_bytes: usize,
    coefficient_bytes: usize,
) -> Result<(), ArithmeticError> {
    let left_held = retained_bytes(left)?;
    let right_held = retained_bytes(right)?;
    check_memory(
        limits,
        &[left_held, right_held, extra, term_bytes, coefficient_bytes],
    )
}

fn check_products<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<(), ArithmeticError> {
    let mut work = 0usize;
    for left_term in &left.terms {
        for right_term in &right.terms {
            check_work(limits, work)?;
            check_monomial_product(&left_term.mono, &right_term.mono)?;
            work = work.saturating_add(1);
        }
    }
    Ok(())
}

fn check_monomial_product(left: &Monomial, right: &Monomial) -> Result<(), ArithmeticError> {
    if left
        .exps
        .iter()
        .zip(&right.exps)
        .any(|(&a, &b)| a.checked_add(b).is_none())
    {
        Err(ArithmeticError::ExponentLimit {
            limit: DEGREE_LIMIT,
        })
    } else {
        Ok(())
    }
}

fn check_power_exponents<D: Domain>(
    poly: &Polynomial<D>,
    exponent: u32,
    limits: &ComputeLimits,
) -> Result<(), ArithmeticError> {
    let mut work = 0usize;
    for variable in 0..poly.ring().nvars() {
        let mut maximum = 0u16;
        for term in &poly.terms {
            check_work(limits, work)?;
            maximum = maximum.max(term.mono.exps[variable]);
            work = work.saturating_add(1);
        }
        if u64::from(maximum) * u64::from(exponent) > u64::from(DEGREE_LIMIT) {
            return Err(ArithmeticError::ExponentLimit {
                limit: DEGREE_LIMIT,
            });
        }
    }
    Ok(())
}

fn product_terms<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<Vec<Term<D>>, ArithmeticError> {
    let pairs = left
        .terms
        .len()
        .checked_mul(right.terms.len())
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    let mut terms = Vec::with_capacity(pairs);
    let ops = left.ring().ops();
    let mut work = 0usize;
    for left_term in &left.terms {
        for right_term in &right.terms {
            check_work(limits, work)?;
            terms.push(Term {
                coeff: ops.mul(&left_term.coeff, &right_term.coeff),
                mono: left_term
                    .mono
                    .checked_mul(&right_term.mono)
                    .map_err(ArithmeticError::from)?,
            });
            work = work.saturating_add(1);
        }
    }
    Ok(terms)
}

fn merge_product_terms<D: Domain>(
    mut terms: Vec<Term<D>>,
    ops: &D::Ops,
    limits: &ComputeLimits,
) -> Result<Vec<Term<D>>, ArithmeticError> {
    check_stop(limits)?;
    // The standard in-place sort cannot poll during its rearrangement. The
    // checks before and after it bound that uninterruptible section.
    terms.sort_unstable_by(|left, right| left.mono.cmp(&right.mono));
    check_stop(limits)?;
    let mut write = 0usize;
    for read in 0..terms.len() {
        check_work(limits, read)?;
        if write > 0 && terms[write - 1].mono == terms[read].mono {
            if !merge_duplicate(&mut terms, write, read, ops) {
                write -= 1;
            }
        } else {
            if write != read {
                terms.swap(write, read);
            }
            write += 1;
        }
    }
    terms.truncate(write);
    Ok(terms)
}

fn merge_duplicate<D: Domain>(
    terms: &mut [Term<D>],
    write: usize,
    read: usize,
    ops: &D::Ops,
) -> bool {
    if let Some(coeff) = ops.add(&terms[write - 1].coeff, &terms[read].coeff) {
        terms[write - 1].coeff = coeff;
        true
    } else {
        false
    }
}

fn one_with_limits<D: Domain>(
    poly: &Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<Polynomial<D>, ArithmeticError> {
    let one_bytes = term_storage::<D>(poly.ring().nvars(), 1)?
        .checked_add(domain_heap_floor::<D>())
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    check_memory(limits, &[retained_bytes(poly)?, one_bytes])?;
    Ok(poly.ring().one())
}

fn pow_nontrivial<D: Domain>(
    poly: &Polynomial<D>,
    exponent: u32,
    limits: &ComputeLimits,
) -> Result<Polynomial<D>, ArithmeticError> {
    let (mut result, mut base) = power_initial(poly, limits)?;
    let mut power = exponent;
    let mut step = 0usize;
    while power != 0 {
        check_work(limits, step)?;
        pow_step(poly, &mut result, &mut base, &mut power, limits)?;
        step = step.saturating_add(1);
    }
    drop(base);
    check_memory(limits, &[retained_bytes(poly)?, retained_bytes(&result)?])?;
    Ok(result)
}

fn power_initial<D: Domain>(
    poly: &Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<(Polynomial<D>, Polynomial<D>), ArithmeticError> {
    let one_bytes = term_storage::<D>(poly.ring().nvars(), 1)?
        .checked_add(domain_heap_floor::<D>())
        .ok_or(ArithmeticError::MemoryLimitExceeded)?;
    let poly_bytes = retained_bytes(poly)?;
    check_memory(limits, &[poly_bytes, one_bytes, poly_bytes])?;
    let result = poly.ring().one();
    let base = poly.clone();
    check_memory(
        limits,
        &[poly_bytes, retained_bytes(&result)?, retained_bytes(&base)?],
    )?;
    Ok((result, base))
}

fn pow_step<D: Domain>(
    poly: &Polynomial<D>,
    result: &mut Polynomial<D>,
    base: &mut Polynomial<D>,
    power: &mut u32,
    limits: &ComputeLimits,
) -> Result<(), ArithmeticError> {
    if *power & 1 == 1 {
        *result = mul_with_extra(result, base, limits, retained_bytes(poly)?)?;
    }
    *power >>= 1;
    if *power != 0 {
        let extra = retained_bytes(poly)?
            .checked_add(retained_bytes(result)?)
            .ok_or(ArithmeticError::MemoryLimitExceeded)?;
        *base = mul_with_extra(base, base, limits, extra)?;
    }
    Ok(())
}

fn coefficient_input_bytes(value: &Coefficient) -> usize {
    match value {
        Coefficient::Small(_) => 0,
        Coefficient::Integer(integer) => limb_bytes(integer),
        Coefficient::Fraction {
            numerator,
            denominator,
        } => limb_bytes(numerator).saturating_add(limb_bytes(denominator)),
        Coefficient::Rational(value) => {
            limb_bytes(value.numer()).saturating_add(limb_bytes(value.denom()))
        }
    }
}
