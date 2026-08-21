//! Division by a list of polynomials, with the quotients recorded.

use std::cmp::Ordering;

use smallvec::SmallVec;

use crate::compute::ComputeError;
use crate::poly::{ExponentOverflow, Monomial, Polynomial, Term};
use crate::ring::field::Felt;

use super::Budget;

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
/// holds, and it stops with [`ComputeError`] when the deadline or the
/// memory limit is passed.
///
/// A divisor whose exponent vectors do not match the width of `f` reduces
/// nothing. The function never panics.
pub(crate) fn divide(
    f: &Polynomial,
    divisors: &[Polynomial],
    modulus: u64,
    budget: &Budget,
) -> Result<(Vec<Polynomial>, Polynomial), ComputeError> {
    budget.check(divisors.len() + 2, f.terms.len())?;
    let mut quotients = vec![f.zero_like(); divisors.len()];
    let mut remainder = f.zero_like();
    let mut work = f.terms.clone();
    let mut scaled: Vec<Term> = Vec::new();
    let mut spare: Vec<Term> = Vec::new();
    let mut quotient_terms = 0usize;
    let mut steps = 0usize;

    while let Some(lead) = work.last().cloned() {
        if steps % DEADLINE_STRIDE == 0 {
            budget.check_deadline()?;
        }
        steps += 1;
        let live = quotient_terms + work.len() + remainder.terms.len();
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
            let coeff = lead.coeff.div(lead_g.coeff, modulus);
            step = Some((index, Term { coeff, mono }));
            break;
        }

        match step {
            Some((index, term)) => {
                negated_multiple(
                    &mut scaled,
                    &divisors[index].terms,
                    term.coeff,
                    &term.mono,
                    modulus,
                )?;
                merge_into(&mut work, &scaled, &mut spare, modulus);
                debug_assert!(
                    work.last().is_none_or(|top| top.mono < lead.mono),
                    "a reduction step must lower the leading monomial"
                );
                let before = quotients[index].terms.len();
                quotients[index].push_term(term, modulus);
                quotient_terms += quotients[index].terms.len() - before;
            }
            None => {
                let Some(term) = work.pop() else { break };
                remainder.push_term(term, modulus);
            }
        }
    }

    Ok((quotients, remainder))
}

/// Write `-coeff * mono * divisor` into `out`, ascending under grevlex.
///
/// One term multiplies every monomial by the same factor, so the order of
/// `divisor` carries over and `out` needs no sort. `coeff` is nonzero and
/// so is every divisor coefficient, so no term drops out. An exponent past
/// the width a `u16` holds returns [`ExponentOverflow`], as `sub_scaled`
/// does.
fn negated_multiple(
    out: &mut Vec<Term>,
    divisor: &[Term],
    coeff: Felt,
    mono: &Monomial,
    p: u64,
) -> Result<(), ExponentOverflow> {
    out.clear();
    out.reserve(divisor.len());
    for term in divisor {
        out.push(Term {
            coeff: term.coeff.mul(coeff, p).neg(p),
            mono: term.mono.checked_mul(mono)?,
        });
    }
    Ok(())
}

/// Apply `work += other`, both term lists ascending under grevlex.
///
/// The merge writes into `spare` and swaps, so the two buffers carry over
/// from step to step and the merge allocates nothing once they are large
/// enough.
fn merge_into(work: &mut Vec<Term>, other: &[Term], spare: &mut Vec<Term>, p: u64) {
    spare.clear();
    spare.reserve(work.len() + other.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < work.len() && j < other.len() {
        match work[i].mono.cmp(&other[j].mono) {
            Ordering::Less => {
                spare.push(copy_term(&work[i]));
                i += 1;
            }
            Ordering::Greater => {
                spare.push(copy_term(&other[j]));
                j += 1;
            }
            Ordering::Equal => {
                let coeff = work[i].coeff.add(other[j].coeff, p);
                if !coeff.is_zero() {
                    spare.push(Term {
                        coeff,
                        mono: copy_mono(&work[i].mono),
                    });
                }
                i += 1;
                j += 1;
            }
        }
    }
    spare.extend(work[i..].iter().map(copy_term));
    spare.extend(other[j..].iter().map(copy_term));
    std::mem::swap(work, spare);
}

/// Copy one monomial by block.
///
/// The merge runs this on its hot path. `from_slice` copies the exponent
/// block in one move; `Clone` would walk it element by element.
fn copy_mono(mono: &Monomial) -> Monomial {
    Monomial {
        exps: SmallVec::from_slice(&mono.exps),
        deg: mono.deg,
    }
}

fn copy_term(term: &Term) -> Term {
    Term {
        coeff: term.coeff,
        mono: copy_mono(&term.mono),
    }
}
