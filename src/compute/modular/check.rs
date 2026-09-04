//! The exact tests T1 and T2 of a lifted basis over `Q`.

use crate::compute::{ComputeError, ComputeLimits, DEGREE_LIMIT, RunError};
use crate::normal_form::{DivisionStop, normal_form};
use crate::poly::Polynomial;
use crate::ring::Rationals;

/// Report whether `candidate` passes T1 and T2 over `Q`.
///
/// T1 reduces every generator modulo `candidate` and T2 reduces every
/// S-polynomial of `candidate` modulo itself. Both remainders must be
/// zero. Together they say that `candidate` is the reduced Gröbner basis
/// of an ideal that contains the ideal of `generators`. They do not say
/// that the two ideals are equal: the basis `{1}` passes both.
///
/// A nonzero remainder is `Ok(false)`, which is evidence about the
/// candidate. An exhausted budget or a monomial past the exponent width
/// stops the whole run and is reported as itself.
pub(crate) fn contains_input(
    generators: &[Polynomial<Rationals>],
    candidate: &[Polynomial<Rationals>],
    limits: &ComputeLimits,
) -> Result<bool, RunError> {
    let ops = match candidate.first() {
        Some(element) => *element.ring().ops(),
        None => return Ok(generators.iter().all(Polynomial::is_zero)),
    };
    if !generators_reduce_to_zero(generators, candidate, limits)? {
        return Ok(false);
    }
    pairs_reduce_to_zero(candidate, &ops, limits)
}

fn generators_reduce_to_zero(
    generators: &[Polynomial<Rationals>],
    candidate: &[Polynomial<Rationals>],
    limits: &ComputeLimits,
) -> Result<bool, RunError> {
    for generator in generators {
        if !normal_form(candidate, generator, limits)
            .map_err(stopped)?
            .is_zero()
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn pairs_reduce_to_zero(
    candidate: &[Polynomial<Rationals>],
    ops: &crate::ring::rational::RationalOps,
    limits: &ComputeLimits,
) -> Result<bool, RunError> {
    for left in 0..candidate.len() {
        for right in (left + 1)..candidate.len() {
            if let Some(stop) = limits.stop() {
                return Err(stop);
            }
            let spoly = candidate[left]
                .s_polynomial_checked(&candidate[right], ops)
                .map_err(|_| exponent_limit())?;
            if !normal_form(candidate, &spoly, limits)
                .map_err(stopped)?
                .is_zero()
            {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Map a division stop onto the report the driver passes on.
fn stopped(stop: DivisionStop) -> RunError {
    match stop {
        DivisionStop::Run(error) => error,
        DivisionStop::ExponentLimit => exponent_limit(),
    }
}

fn exponent_limit() -> RunError {
    RunError::Compute(ComputeError::ExponentLimit {
        limit: DEGREE_LIMIT,
    })
}
