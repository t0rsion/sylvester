//! Membership and S-pair representations over the final basis.

use crate::certificate::{CertifyError, EmitterFault};
use crate::poly::Polynomial;
use crate::ring::PrimeOps;

use super::WriterBudget;
use super::divide::divide;

/// One S-pair entry: the pair and one cofactor per basis element.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SpairEntry {
    pub(crate) i: usize,
    pub(crate) j: usize,
    pub(crate) cofactors: Vec<Polynomial>,
}

/// Divide every input polynomial by the basis and keep the quotients.
///
/// Entry i holds one cofactor per basis element, so that
/// f_i = sum_j q_ij * g_j. A nonzero remainder is a defect of the
/// candidate: it reports [`EmitterFault::InputHasRemainder`] and writes
/// nothing.
///
/// The work holds to `budget`, which also carries the quotients of the
/// entries already built.
pub(crate) fn membership_representations(
    input: &[Polynomial],
    basis: &[Polynomial],
    ops: &PrimeOps,
    budget: &mut WriterBudget,
) -> Result<Vec<Vec<Polynomial>>, CertifyError> {
    let mut out = Vec::with_capacity(input.len());
    for (index, f) in input.iter().enumerate() {
        budget.check_stop()?;
        let (quotients, remainder) = divide(f, basis, ops, budget)?;
        if !remainder.is_zero() {
            return Err(CertifyError::Emitter(EmitterFault::InputHasRemainder {
                input: index,
            }));
        }
        budget.hold_polys(&quotients)?;
        out.push(quotients);
    }
    Ok(out)
}

/// Build one entry per pair whose leading monomials share a variable.
///
/// A pair with coprime leading monomials carries no entry. The product
/// criterion justifies the omission and the verifier re-enumerates the
/// pairs itself. A nonzero remainder proves the basis is not a Gröbner
/// basis: it reports [`EmitterFault::NotAGroebnerBasis`] with the pair.
///
/// The work is quadratic in the size of the basis, and each entry holds one
/// cofactor per basis element. It holds to `budget` at every pair.
pub(crate) fn spair_representations(
    basis: &[Polynomial],
    ops: &PrimeOps,
    budget: &mut WriterBudget,
) -> Result<Vec<SpairEntry>, CertifyError> {
    let mut out = Vec::new();
    for i in 0..basis.len() {
        budget.check_stop()?;
        for j in (i + 1)..basis.len() {
            let (Some(lm_i), Some(lm_j)) = (basis[i].lm(), basis[j].lm()) else {
                continue;
            };
            if lm_i.nvars() != lm_j.nvars() || lm_i.is_coprime(lm_j) {
                continue;
            }
            let spoly = basis[i].s_polynomial(&basis[j], ops);
            let (cofactors, remainder) = divide(&spoly, basis, ops, budget)?;
            if !remainder.is_zero() {
                return Err(CertifyError::Emitter(EmitterFault::NotAGroebnerBasis {
                    i,
                    j,
                }));
            }
            debug_assert!(
                bounded_by(&cofactors, basis, &spoly),
                "every summand of an S-pair identity stays at or below the S-polynomial"
            );
            budget.hold_polys(&cofactors)?;
            out.push(SpairEntry { i, j, cofactors });
        }
    }
    Ok(out)
}

/// Report whether every nonzero product q_l * g_l stays at or below the
/// leading monomial of `target`.
fn bounded_by(cofactors: &[Polynomial], basis: &[Polynomial], target: &Polynomial) -> bool {
    cofactors
        .iter()
        .zip(basis)
        .all(|(cofactor, element)| match (cofactor.lm(), element.lm()) {
            (Some(lm_q), Some(lm_g)) => target.lm().is_some_and(|bound| lm_q.mul(lm_g) <= *bound),
            _ => true,
        })
}
