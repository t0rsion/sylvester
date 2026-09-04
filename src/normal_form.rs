//! Division by a Gröbner basis, and the check that a list is one.
//!
//! One algorithm serves every caller: [`GroebnerBasis::normal_form`],
//! [`GroebnerBasis::contains`], and the checked constructor
//! `GroebnerBasis::from_polynomials`. It is generic over the domain, so
//! `F_p` and `Q` run the same code. The verifiers keep their own division
//! and call nothing here, which `tests/isolation.rs` holds.
//!
//! [`GroebnerBasis::normal_form`]: crate::GroebnerBasis::normal_form
//! [`GroebnerBasis::contains`]: crate::GroebnerBasis::contains

use std::fmt;
use std::mem::size_of;

use crate::compute::{ComputeError, ComputeLimits, DEGREE_LIMIT, RunError};
use crate::poly::{ExponentOverflow, Polynomial, Term, divide_coefficients, heap_exps_bytes};
use crate::ring::{Domain, DomainOps, PolynomialRing};

/// Why a normal form is not computed.
///
/// [`NormalFormError::Timeout`] and
/// [`NormalFormError::MemoryLimitExceeded`] report an exhausted budget.
/// [`NormalFormError::ExponentLimit`] reports a monomial the division
/// needs and one exponent cannot hold. None of them says anything about
/// the ideal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NormalFormError {
    /// The polynomial and the basis belong to different rings.
    RingMismatch,
    /// A monomial the division needs holds an exponent above `limit`.
    ///
    /// Dividing `x^65535*y` by `x^65535 + y^65535` reaches it: the
    /// quotient monomial is `y`, and the tail multiple is `y^65536`.
    ExponentLimit {
        /// The largest value one exponent holds.
        limit: u32,
    },
    /// The deadline passed.
    Timeout,
    /// The live data of the division passed the memory limit.
    MemoryLimitExceeded,
}

impl fmt::Display for NormalFormError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NormalFormError::RingMismatch => {
                f.write_str("the polynomial and the basis belong to different rings")
            }
            NormalFormError::ExponentLimit { limit } => write!(
                f,
                "the division reached an exponent above {limit}, which is the largest one exponent holds"
            ),
            NormalFormError::Timeout => f.write_str("the division passed its deadline"),
            NormalFormError::MemoryLimitExceeded => {
                f.write_str("the division passed its memory limit")
            }
        }
    }
}

impl std::error::Error for NormalFormError {}

/// Why a list of polynomials is not accepted as a reduced Gröbner basis.
///
/// `GroebnerBasis::from_polynomials` reports it. The first six variants
/// name a property the list does not have, with the index or the index
/// pair that shows it. The last three report a limit of the check itself
/// and say nothing about the list.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BasisError {
    /// A polynomial belongs to a different ring than the basis.
    RingMismatch,
    /// The polynomial at `index` is zero.
    ZeroPolynomial {
        /// The position in the list the caller gave.
        index: usize,
    },
    /// The polynomial at `index` has a leading coefficient other than 1.
    NotMonic {
        /// The position in the list the caller gave.
        index: usize,
    },
    /// The leading monomial at `index` is not smaller than the one before
    /// it.
    ///
    /// A reduced basis runs strictly descending, so two equal leading
    /// monomials are reported here too.
    NotSorted {
        /// The position in the list the caller gave.
        index: usize,
    },
    /// A monomial of the polynomial at `index` is divisible by the
    /// leading monomial of another element.
    NotInterreduced {
        /// The position in the list the caller gave.
        index: usize,
    },
    /// The S-polynomial of the pair does not reduce to zero, so the list
    /// is not a Gröbner basis.
    ///
    /// The pair is the first one in ascending order over the caller's own
    /// indices, so the report is a function of the list and not of the
    /// order the check walks.
    NotGroebner {
        /// The smaller position of the pair.
        left: usize,
        /// The larger position of the pair.
        right: usize,
    },
    /// A monomial the check needs holds an exponent above `limit`.
    ///
    /// Building an S-polynomial multiplies the tail of each side by a
    /// quotient monomial, and that product can pass the width even when
    /// both leading monomials and their least common multiple fit.
    ExponentLimit {
        /// The largest value one exponent holds.
        limit: u32,
    },
    /// The deadline passed.
    Timeout,
    /// The live data of the check passed the memory limit.
    MemoryLimitExceeded,
}

impl fmt::Display for BasisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BasisError::RingMismatch => {
                f.write_str("a polynomial belongs to a different ring than the basis")
            }
            BasisError::ZeroPolynomial { index } => {
                write!(f, "the polynomial at index {index} is zero")
            }
            BasisError::NotMonic { index } => write!(
                f,
                "the polynomial at index {index} has a leading coefficient other than 1"
            ),
            BasisError::NotSorted { index } => write!(
                f,
                "the leading monomial at index {index} is not smaller than the one before it"
            ),
            BasisError::NotInterreduced { index } => write!(
                f,
                "a monomial of the polynomial at index {index} is divisible by the leading monomial of another element"
            ),
            BasisError::NotGroebner { left, right } => write!(
                f,
                "the S-polynomial of the elements at indices {left} and {right} does not reduce to zero"
            ),
            BasisError::ExponentLimit { limit } => write!(
                f,
                "the check reached an exponent above {limit}, which is the largest one exponent holds"
            ),
            BasisError::Timeout => f.write_str("the check passed its deadline"),
            BasisError::MemoryLimitExceeded => f.write_str("the check passed its memory limit"),
        }
    }
}

impl std::error::Error for BasisError {}

/// Why a division stops before it has a remainder.
///
/// The public errors name the same cases. The division returns this one
/// and each caller maps it: [`NormalFormError`] for a normal form,
/// [`BasisError`] for a basis check.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum DivisionStop {
    /// The budget ran out, or the cancellation flag was set.
    Run(RunError),
    /// A monomial the division needs holds an exponent past the width.
    ExponentLimit,
}

impl From<ExponentOverflow> for DivisionStop {
    fn from(_: ExponentOverflow) -> Self {
        DivisionStop::ExponentLimit
    }
}

impl From<DivisionStop> for NormalFormError {
    fn from(stop: DivisionStop) -> Self {
        match stop {
            // A cancelled division reports a timeout, as
            // `RunError::reported` does: it stopped before it had a
            // remainder.
            DivisionStop::Run(RunError::Compute(ComputeError::MemoryLimitExceeded)) => {
                NormalFormError::MemoryLimitExceeded
            }
            DivisionStop::Run(_) => NormalFormError::Timeout,
            DivisionStop::ExponentLimit => NormalFormError::ExponentLimit {
                limit: DEGREE_LIMIT,
            },
        }
    }
}

impl From<DivisionStop> for BasisError {
    fn from(stop: DivisionStop) -> Self {
        match NormalFormError::from(stop) {
            NormalFormError::MemoryLimitExceeded => BasisError::MemoryLimitExceeded,
            NormalFormError::ExponentLimit { limit } => BasisError::ExponentLimit { limit },
            _ => BasisError::Timeout,
        }
    }
}

impl From<ExponentOverflow> for BasisError {
    fn from(_: ExponentOverflow) -> Self {
        BasisError::ExponentLimit {
            limit: DEGREE_LIMIT,
        }
    }
}

/// Reduce `f` modulo `basis` and return the remainder.
///
/// The algorithm is full division. While the largest monomial of the
/// working polynomial is divisible by the leading monomial of some
/// element, the matching multiple of that element is subtracted. A
/// monomial no element divides moves to the remainder. `basis` must be a
/// Gröbner basis, and then the remainder is the unique normal form
/// whatever order the divisors are picked in. This takes the first
/// divisor in the order of `basis`, so the steps are deterministic as
/// well.
///
/// The limits stop the division between two reduction steps. One step
/// over `Q` is one exact subtraction, so that step is the granularity of
/// the deadline. The memory limit is read twice per step: once on the
/// working value and the remainder, and once more on the replacement the
/// step is about to allocate next to them.
pub(crate) fn normal_form<D: Domain>(
    basis: &[Polynomial<D>],
    f: &Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<Polynomial<D>, DivisionStop> {
    let ring = f.ring();
    let ops = ring.ops();
    let mut working = f.clone();
    // The remainder grows by the largest term left, so it is built
    // descending and turned around once.
    let mut remainder: Vec<Term<D>> = Vec::new();
    let mut remainder_bytes = 0usize;

    while let Some(lead) = working.lt().cloned() {
        let held = working.heap_bytes().saturating_add(remainder_bytes);
        check(limits, held)?;

        let mut reducer = None;
        for g in basis {
            let Some(lead_g) = g.lt() else { continue };
            if lead_g.mono.divides(&lead.mono) {
                reducer = Some((g, lead_g));
                break;
            }
        }

        match reducer {
            Some((g, lead_g)) => {
                let multiple = lead
                    .mono
                    .quotient(&lead_g.mono)
                    // divides() implies a quotient exists.
                    .expect("divides() implies quotient()");
                let scale = divide_coefficients::<D>(ops, &lead.coeff, &lead_g.coeff);
                // The step builds a whole new polynomial while the working
                // one is still live, so the replacement is charged before
                // it is allocated.
                check(
                    limits,
                    held.saturating_add(working.sub_scaled_bytes(g, &scale, ops)),
                )?;
                working = working.sub_scaled_checked(g, &scale, &multiple, ops)?;
            }
            None => {
                let term = working
                    .pop_lt()
                    // The loop condition read the leading term.
                    .expect("the working polynomial holds a leading term");
                remainder_bytes =
                    remainder_bytes.saturating_add(term_bytes::<D>(&term, ops, ring.nvars()));
                remainder.push(term);
            }
        }
    }

    remainder.reverse();
    Ok(Polynomial::from_sorted_terms(ring.clone(), remainder))
}

/// Check that `polynomials` is the reduced Gröbner basis of the ideal it
/// generates.
///
/// The checks run in one order, so the report is a function of the list:
/// the ring of every element, then no zero element, then monic, then
/// strictly descending leading monomials, then interreduced, then every
/// S-polynomial reduces to zero. The last check is the Buchberger
/// criterion over the whole list, which is what makes the list a Gröbner
/// basis of its own ideal.
///
/// A list this accepts is a check and not a certificate. It says nothing
/// about the ideal a caller meant, only about the list.
pub(crate) fn check_basis<D: Domain>(
    ring: &PolynomialRing<D>,
    polynomials: &[Polynomial<D>],
    limits: &ComputeLimits,
) -> Result<(), BasisError> {
    let ops = ring.ops();
    check_elements(ring, polynomials, limits, ops)?;
    check_order(polynomials)?;
    check_interreduced(polynomials, limits)?;
    check_pairs(polynomials, limits, ops)
}

fn check_elements<D: Domain>(
    ring: &PolynomialRing<D>,
    polynomials: &[Polynomial<D>],
    limits: &ComputeLimits,
    ops: &D::Ops,
) -> Result<(), BasisError> {
    for (index, poly) in polynomials.iter().enumerate() {
        stop(limits)?;
        if poly.ring() != ring {
            return Err(BasisError::RingMismatch);
        }
        if poly.is_zero() {
            return Err(BasisError::ZeroPolynomial { index });
        }
        let lead = poly
            .lt()
            .expect("a nonzero polynomial holds a leading term");
        if !ops.is_one(&lead.coeff) {
            return Err(BasisError::NotMonic { index });
        }
    }
    Ok(())
}

fn check_order<D: Domain>(polynomials: &[Polynomial<D>]) -> Result<(), BasisError> {
    for index in 1..polynomials.len() {
        let previous = polynomials[index - 1]
            .lm()
            .expect("a nonzero polynomial holds a leading monomial");
        let current = polynomials[index]
            .lm()
            .expect("a nonzero polynomial holds a leading monomial");
        if current >= previous {
            return Err(BasisError::NotSorted { index });
        }
    }
    Ok(())
}

fn check_interreduced<D: Domain>(
    polynomials: &[Polynomial<D>],
    limits: &ComputeLimits,
) -> Result<(), BasisError> {
    for (index, poly) in polynomials.iter().enumerate() {
        stop(limits)?;
        for (other, divisor) in polynomials.iter().enumerate() {
            if other == index {
                continue;
            }
            let lead = divisor
                .lm()
                .expect("a nonzero polynomial holds a leading monomial");
            if poly.terms.iter().any(|term| lead.divides(&term.mono)) {
                return Err(BasisError::NotInterreduced { index });
            }
        }
    }
    Ok(())
}

fn check_pairs<D: Domain>(
    polynomials: &[Polynomial<D>],
    limits: &ComputeLimits,
    ops: &D::Ops,
) -> Result<(), BasisError> {
    for left in 0..polynomials.len() {
        for right in (left + 1)..polynomials.len() {
            stop(limits)?;
            let spoly = polynomials[left].s_polynomial_checked(&polynomials[right], ops)?;
            if !normal_form(polynomials, &spoly, limits)?.is_zero() {
                return Err(BasisError::NotGroebner { left, right });
            }
        }
    }
    Ok(())
}

/// Report why the division stops now, or `None` to carry on.
fn check(limits: &ComputeLimits, bytes: usize) -> Result<(), DivisionStop> {
    if let Some(stop) = limits.stop() {
        return Err(DivisionStop::Run(stop));
    }
    if let Some(limit) = limits.memory
        && bytes > limit
    {
        return Err(DivisionStop::Run(RunError::Compute(
            ComputeError::MemoryLimitExceeded,
        )));
    }
    Ok(())
}

/// Report why the check stops now, or `None` to carry on.
///
/// The checks before the S-polynomials hold no data of their own, so the
/// memory limit is read inside the division alone.
fn stop(limits: &ComputeLimits) -> Result<(), BasisError> {
    match limits.stop() {
        Some(stop) => Err(BasisError::from(DivisionStop::Run(stop))),
        None => Ok(()),
    }
}

/// The bytes one term of the remainder holds.
fn term_bytes<D: Domain>(term: &Term<D>, ops: &D::Ops, nvars: usize) -> usize {
    size_of::<Term<D>>()
        .saturating_add(heap_exps_bytes(nvars))
        .saturating_add(ops.heap_bytes(&term.coeff))
}
