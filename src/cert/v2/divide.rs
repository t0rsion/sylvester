//! Division traces over the final basis (contract section 4.5).

use crate::certificate::CertifyError;
use crate::poly::{Monomial, Polynomial};
use crate::ring::PrimeOps;

use super::WriterBudget;

/// Steps between two reads of the clock.
///
/// One step is one reduction. The clock is the expensive part of the
/// check, so the division reads it once per stride.
const DEADLINE_STRIDE: usize = 64;

/// One step of a division trace.
///
/// The step subtracts `lc(residual) * mono * basis[element]`. The
/// certificate carries no coefficient, because the leading coefficient of
/// the residual is the only one that cancels the lead.
pub(super) struct Step {
    /// The multiplier.
    pub(super) mono: Monomial,
    /// The index of the element in the basis.
    pub(super) element: u32,
}

/// Divide `f` by `basis` and record the steps.
///
/// The reducer at each step is the element with the smallest index whose
/// leading monomial divides the leading monomial of the residual (contract
/// section 9.5). `None` reports a residual that no element reduces, which
/// is a defect in the candidate the writer read.
///
/// Every element of `basis` must be monic. The caller checks that once,
/// because the step subtracts the leading coefficient of the residual and
/// nothing else.
pub(super) fn divide(
    f: &Polynomial,
    basis: &[Polynomial],
    ops: &PrimeOps,
    budget: &WriterBudget,
) -> Result<Option<Vec<Step>>, CertifyError> {
    let mut residual = f.clone();
    let mut steps: Vec<Step> = Vec::new();
    let mut count = 0usize;
    while let Some(lead) = residual.lt().cloned() {
        if count.is_multiple_of(DEADLINE_STRIDE) {
            budget.check_stop()?;
        }
        count += 1;
        budget.check(1, residual.terms.len() + steps.len())?;

        let mut reducer = None;
        for (index, g) in basis.iter().enumerate() {
            let Some(lm) = g.lm() else { continue };
            if lm.nvars() != lead.mono.nvars() {
                continue;
            }
            if let Some(mono) = lead.mono.quotient(lm) {
                reducer = Some((index, mono));
                break;
            }
        }
        let Some((index, mono)) = reducer else {
            return Ok(None);
        };
        residual = residual.sub_scaled(&basis[index], &lead.coeff, &mono, ops);
        debug_assert!(
            residual.lm().is_none_or(|lm| *lm < lead.mono),
            "a division step must lower the leading monomial"
        );
        steps.push(Step {
            mono,
            element: index as u32,
        });
    }
    Ok(Some(steps))
}
