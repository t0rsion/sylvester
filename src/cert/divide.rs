//! Division by a list of polynomials, with the quotients recorded.

use crate::certificate::CertifyError;
use crate::poly::{Polynomial, Term};
use crate::ring::PrimeOps;

use super::WriterBudget;

/// Steps between two reads of the clock.
///
/// One step is a reduction or a move to the remainder. The clock is the
/// expensive part of the check, so the division reads it once per stride.
const DEADLINE_STRIDE: usize = 64;

/// Divide `f` by `divisors` and record one quotient per divisor.
///
/// Returns the quotients q_1..q_k and the remainder r of the identity
/// f = sum_j q_j * g_j + r.
///
/// This is the full division algorithm: every step reduces the leading
/// term of the working polynomial by a divisor, or moves that term to the
/// remainder. Both steps lower the leading monomial of the working
/// polynomial, so the division terminates and every recorded quotient term
/// m of divisor j satisfies m * lm(g_j) <= lm(f). A product of a quotient
/// with its divisor cannot cancel its own leading term, so
/// lm(q_j * g_j) <= lm(f) as well. That bound is what the S-pair entries of
/// the certificate need.
///
/// The division holds to `budget`. It charges the quotients, the working
/// polynomial, and the remainder against the data the emitter already
/// holds, and it stops with [`CertifyError::WriterExhausted`] when the
/// deadline or the memory limit is passed.
///
/// A divisor whose exponent vectors do not match the width of `f` reduces
/// nothing. The function never panics.
pub(crate) fn divide(
    f: &Polynomial,
    divisors: &[Polynomial],
    ops: &PrimeOps,
    budget: &WriterBudget,
) -> Result<(Vec<Polynomial>, Polynomial), CertifyError> {
    budget.check(divisors.len() + 2, f.terms.len())?;
    let mut quotients = vec![f.zero_like(); divisors.len()];
    let mut remainder = f.zero_like();
    let mut work = f.clone();
    let mut quotient_terms = 0usize;
    let mut steps = 0usize;

    while let Some(lead) = work.lt().cloned() {
        if steps.is_multiple_of(DEADLINE_STRIDE) {
            budget.check_stop()?;
        }
        steps += 1;
        let live = quotient_terms + work.terms.len() + remainder.terms.len();
        budget.check(divisors.len() + 2, live)?;

        let mut step = None;
        for (index, g) in divisors.iter().enumerate() {
            let Some(lead_g) = g.lt() else { continue };
            if lead_g.mono.nvars() != lead.mono.nvars() {
                continue;
            }
            let Some(mono) = lead.mono.quotient(&lead_g.mono) else {
                continue;
            };
            let coeff = lead.coeff.div(lead_g.coeff, ops.modulus());
            step = Some((index, Term { coeff, mono }));
            break;
        }

        match step {
            Some((index, term)) => {
                work = work.sub_scaled(&divisors[index], &term.coeff, &term.mono, ops);
                debug_assert!(
                    work.lm().is_none_or(|lm| *lm < lead.mono),
                    "a reduction step must lower the leading monomial"
                );
                let before = quotients[index].terms.len();
                quotients[index].push_term(term, ops);
                quotient_terms = (quotient_terms + quotients[index].terms.len()) - before;
            }
            None => {
                let Some(term) = work.pop_lt() else { break };
                remainder.push_term(term, ops);
            }
        }
    }

    Ok((quotients, remainder))
}
