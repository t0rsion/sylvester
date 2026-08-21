//! The verifier obligations, in the order the contract lists them.

use std::cmp::Ordering;

use super::algebra::{self, Mono, Poly};
use super::error::{BasisFault, Location, PolyFault, VerifyError};
use super::json::{MAX_NVARS, RawCert};
use super::limits::{Budget, Limits};

pub(crate) const SCHEMA: &str = "sylv-gb-cert-v1";
pub(crate) const ORDER: &str = "grevlex-v1";

/// The largest modulus the contract allows.
const MAX_MODULUS: u64 = (1 << 31) - 1;

/// The decoded certificate, after every obligation holds.
pub(crate) struct Accepted {
    pub modulus: u64,
    pub nvars: usize,
    pub input: Vec<Poly>,
    pub basis: Vec<Poly>,
}

/// Run every obligation on a decoded certificate.
pub(crate) fn check(raw: RawCert, limits: &Limits) -> Result<Accepted, VerifyError> {
    if raw.schema != SCHEMA {
        return Err(VerifyError::Schema { found: raw.schema });
    }
    if raw.order != ORDER {
        return Err(VerifyError::Order { found: raw.order });
    }
    let modulus = check_modulus(raw.modulus)?;
    let nvars = check_nvars(raw.nvars)?;
    // One budget runs every obligation: it charges the arithmetic buffers
    // in bytes at this certificate's width, and it polls the deadline at a
    // fixed stride of work inside every loop.
    let mut budget = limits.budget(nvars);
    check_encodings(&raw, nvars, modulus, &mut budget)?;
    check_basis(&raw.basis, &mut budget)?;
    check_origin(&raw, modulus, &mut budget)?;
    check_membership(&raw, modulus, &mut budget)?;
    check_spairs(&raw, modulus, &mut budget)?;
    Ok(Accepted {
        modulus,
        nvars,
        input: raw.input,
        basis: raw.basis,
    })
}

fn check_modulus(found: u64) -> Result<u64, VerifyError> {
    if !(2..=MAX_MODULUS).contains(&found) {
        return Err(VerifyError::Modulus {
            found,
            composite: false,
        });
    }
    if !algebra::is_prime(found) {
        return Err(VerifyError::Modulus {
            found,
            composite: true,
        });
    }
    Ok(found)
}

fn check_nvars(found: u64) -> Result<usize, VerifyError> {
    if found > MAX_NVARS {
        return Err(VerifyError::Nvars {
            found,
            max: MAX_NVARS,
        });
    }
    Ok(found as usize)
}

/// Obligation 2: every polynomial is a canonical encoding.
fn check_encodings(
    raw: &RawCert,
    nvars: usize,
    modulus: u64,
    budget: &mut Budget,
) -> Result<(), VerifyError> {
    for (index, poly) in raw.input.iter().enumerate() {
        check_poly(poly, Location::Input(index), nvars, modulus, budget)?;
    }
    for (index, poly) in raw.basis.iter().enumerate() {
        check_poly(poly, Location::Basis(index), nvars, modulus, budget)?;
    }
    for (basis, entry) in raw.origin.iter().enumerate() {
        for (input, poly) in entry.iter().enumerate() {
            check_poly(
                poly,
                Location::Origin { basis, input },
                nvars,
                modulus,
                budget,
            )?;
        }
    }
    for (input, entry) in raw.membership.iter().enumerate() {
        for (basis, poly) in entry.iter().enumerate() {
            check_poly(
                poly,
                Location::Membership { input, basis },
                nvars,
                modulus,
                budget,
            )?;
        }
    }
    for entry in &raw.spairs {
        let i = usize::try_from(entry.i).unwrap_or(usize::MAX);
        let j = usize::try_from(entry.j).unwrap_or(usize::MAX);
        for (basis, poly) in entry.cofactors.iter().enumerate() {
            check_poly(
                poly,
                Location::SpairCofactor { i, j, basis },
                nvars,
                modulus,
                budget,
            )?;
        }
    }
    Ok(())
}

fn check_poly(
    poly: &Poly,
    at: Location,
    nvars: usize,
    modulus: u64,
    budget: &mut Budget,
) -> Result<(), VerifyError> {
    let fault = |fault| Err(VerifyError::Polynomial { at, fault });
    for (index, term) in poly.terms().iter().enumerate() {
        budget.step()?;
        if term.coeff() == 0 {
            return fault(PolyFault::CoefficientZero { term: index });
        }
        if term.coeff() >= modulus {
            return fault(PolyFault::CoefficientOutOfRange {
                term: index,
                coeff: term.coeff(),
            });
        }
        if term.mono().nvars() != nvars {
            return fault(PolyFault::ExponentCount {
                term: index,
                found: term.mono().nvars(),
                expected: nvars,
            });
        }
        if index > 0 {
            let previous = poly.terms()[index - 1].mono();
            match algebra::cmp_grevlex(previous, term.mono()) {
                Ordering::Greater => {}
                Ordering::Equal => {
                    return fault(PolyFault::DuplicateMonomial { term: index });
                }
                Ordering::Less => return fault(PolyFault::NotDescending { term: index }),
            }
        }
    }
    Ok(())
}

/// Collect the leading monomial of every basis element.
fn leading_monomials<'a>(
    basis: &'a [Poly],
    budget: &mut Budget,
) -> Result<Vec<&'a Mono>, VerifyError> {
    let mut lms = Vec::with_capacity(basis.len());
    for (index, element) in basis.iter().enumerate() {
        budget.step()?;
        match element.lm() {
            Some(lm) => lms.push(lm),
            None => {
                return Err(VerifyError::Basis {
                    index,
                    fault: BasisFault::Zero,
                });
            }
        }
    }
    Ok(lms)
}

/// Obligation 3: the basis is sorted, monic, interreduced, and minimal.
///
/// A strictly descending basis is duplicate-free: two equal elements share a
/// leading monomial, which the sort check rejects.
fn check_basis(basis: &[Poly], budget: &mut Budget) -> Result<(), VerifyError> {
    let lms = leading_monomials(basis, budget)?;
    for index in 1..basis.len() {
        budget.step()?;
        if algebra::cmp_grevlex(lms[index - 1], lms[index]) != Ordering::Greater {
            return Err(VerifyError::Basis {
                index,
                fault: BasisFault::NotDescending,
            });
        }
    }
    for (index, element) in basis.iter().enumerate() {
        budget.step()?;
        match element.lc() {
            Some(1) => {}
            Some(lc) => {
                return Err(VerifyError::Basis {
                    index,
                    fault: BasisFault::NotMonic { lc },
                });
            }
            None => {
                return Err(VerifyError::Basis {
                    index,
                    fault: BasisFault::Zero,
                });
            }
        }
    }
    for (index, element) in basis.iter().enumerate() {
        for (term, tail) in element.terms().iter().enumerate().skip(1) {
            for (by, lm) in lms.iter().enumerate() {
                budget.step()?;
                if by != index && lm.divides(tail.mono()) {
                    return Err(VerifyError::Basis {
                        index,
                        fault: BasisFault::TailReducible { term, by },
                    });
                }
            }
        }
    }
    for (index, lm) in lms.iter().enumerate() {
        for (by, other) in lms.iter().enumerate() {
            budget.step()?;
            if by != index && other.divides(lm) {
                return Err(VerifyError::Basis {
                    index,
                    fault: BasisFault::LeadDivisible { by },
                });
            }
        }
    }
    Ok(())
}

/// The sum of the products of the cofactors with the generators.
///
/// The sum takes one product at a time, so at most two results live at
/// once.
fn combination(
    cofactors: &[Poly],
    generators: &[Poly],
    modulus: u64,
    budget: &mut Budget,
) -> Result<Poly, VerifyError> {
    let mut sum = Poly::zero();
    for (cofactor, generator) in cofactors.iter().zip(generators) {
        budget.step()?;
        if cofactor.is_zero() || generator.is_zero() {
            continue;
        }
        let live = sum.terms().len();
        let product = algebra::mul(cofactor, generator, modulus, budget, live)?;
        sum = algebra::add(sum, product, modulus, budget, 0)?;
    }
    Ok(sum)
}

/// Check that every entry's combination of `generators` equals its target.
///
/// The outer array holds one entry per target, each entry one cofactor per
/// generator. A wrong count is a `CountMismatch` named by `names.0` (the
/// array) or `names.1` (one entry). A combination that misses its target is
/// the identity fault `fault` names.
fn check_combinations(
    entries: &[Vec<Poly>],
    generators: &[Poly],
    targets: &[Poly],
    names: (&'static str, &'static str),
    modulus: u64,
    budget: &mut Budget,
    fault: impl Fn(usize) -> VerifyError,
) -> Result<(), VerifyError> {
    if entries.len() != targets.len() {
        return Err(VerifyError::CountMismatch {
            what: names.0,
            index: None,
            found: entries.len(),
            expected: targets.len(),
        });
    }
    for (index, cofactors) in entries.iter().enumerate() {
        budget.step()?;
        if cofactors.len() != generators.len() {
            return Err(VerifyError::CountMismatch {
                what: names.1,
                index: Some(index),
                found: cofactors.len(),
                expected: generators.len(),
            });
        }
        if combination(cofactors, generators, modulus, budget)? != targets[index] {
            return Err(fault(index));
        }
    }
    Ok(())
}

/// Obligation 4: every basis element is a combination of the input.
fn check_origin(raw: &RawCert, modulus: u64, budget: &mut Budget) -> Result<(), VerifyError> {
    check_combinations(
        &raw.origin,
        &raw.input,
        &raw.basis,
        ("origin", "origin entry"),
        modulus,
        budget,
        |basis| VerifyError::OriginIdentity { basis },
    )
}

/// Obligation 5: every input polynomial is a combination of the basis.
fn check_membership(raw: &RawCert, modulus: u64, budget: &mut Budget) -> Result<(), VerifyError> {
    check_combinations(
        &raw.membership,
        &raw.basis,
        &raw.input,
        ("membership", "membership entry"),
        modulus,
        budget,
        |input| VerifyError::MembershipIdentity { input },
    )
}

/// Obligation 6: every pair is covered, and every entry holds.
fn check_spairs(raw: &RawCert, modulus: u64, budget: &mut Budget) -> Result<(), VerifyError> {
    let basis = &raw.basis;
    let lms = leading_monomials(basis, budget)?;
    let count = basis.len();
    let mut pairs: Vec<(usize, usize)> = Vec::with_capacity(raw.spairs.len());
    for (entry_index, entry) in raw.spairs.iter().enumerate() {
        budget.step()?;
        if entry.i >= count as u64 {
            return Err(VerifyError::IndexOutOfRange {
                what: "the first S-pair index",
                index: entry.i,
                bound: count,
            });
        }
        if entry.j >= count as u64 {
            return Err(VerifyError::IndexOutOfRange {
                what: "the second S-pair index",
                index: entry.j,
                bound: count,
            });
        }
        let (i, j) = (entry.i as usize, entry.j as usize);
        if i >= j {
            return Err(VerifyError::SpairIndexOrder { i, j });
        }
        if entry.cofactors.len() != count {
            return Err(VerifyError::CountMismatch {
                what: "S-pair entry",
                index: Some(entry_index),
                found: entry.cofactors.len(),
                expected: count,
            });
        }
        if let Some(&previous) = pairs.last() {
            if previous == (i, j) {
                return Err(VerifyError::SpairDuplicate { i, j });
            }
            if (i, j) < previous {
                return Err(VerifyError::SpairUnsorted { i, j });
            }
        }
        pairs.push((i, j));
    }

    let mut next = 0usize;
    for (i, left) in basis.iter().enumerate() {
        for (j, right) in basis.iter().enumerate().skip(i + 1) {
            budget.step()?;
            if pairs.get(next) == Some(&(i, j)) {
                check_spair(
                    (i, j),
                    left,
                    right,
                    &raw.spairs[next].cofactors,
                    basis,
                    modulus,
                    budget,
                )?;
                next += 1;
            } else if !lms[i].is_coprime(lms[j]) {
                return Err(VerifyError::SpairMissing { i, j });
            }
        }
    }
    Ok(())
}

/// Check one S-pair identity and its leading-monomial bounds.
///
/// The check holds one product at a time. It checks the bound on a product,
/// then adds the product to the running sum.
///
/// The zero S-polynomial has no leading monomial. It bounds no nonzero
/// summand, so every cofactor must be zero.
fn check_spair(
    pair: (usize, usize),
    left: &Poly,
    right: &Poly,
    cofactors: &[Poly],
    basis: &[Poly],
    modulus: u64,
    budget: &mut Budget,
) -> Result<(), VerifyError> {
    let (i, j) = pair;
    let target = algebra::spoly(left, right, modulus, budget, 0)?;
    let mut sum = Poly::zero();
    for (summand, (cofactor, element)) in cofactors.iter().zip(basis).enumerate() {
        budget.step()?;
        if cofactor.is_zero() || element.is_zero() {
            continue;
        }
        let live = target.terms().len() + sum.terms().len();
        let product = algebra::mul(cofactor, element, modulus, budget, live)?;
        if let Some(lm) = product.lm() {
            let above = match target.lm() {
                Some(bound) => algebra::cmp_grevlex(lm, bound) == Ordering::Greater,
                None => true,
            };
            if above {
                return Err(VerifyError::SpairBound { i, j, summand });
            }
        }
        sum = algebra::add(sum, product, modulus, budget, target.terms().len())?;
    }
    if sum != target {
        return Err(VerifyError::SpairIdentity { i, j });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verify::algebra::{Poly, Term};

    fn budget() -> Budget {
        Limits::default().budget(2)
    }

    fn poly(terms: &[(u64, &[u32])]) -> Poly {
        Poly::new(
            terms
                .iter()
                .map(|(c, e)| Term::new(*c, Mono::new(e.to_vec())))
                .collect::<Vec<_>>(),
        )
    }

    #[test]
    fn a_composite_modulus_is_named_composite() {
        assert_eq!(
            check_modulus(9),
            Err(VerifyError::Modulus {
                found: 9,
                composite: true
            })
        );
        assert_eq!(
            check_modulus(1),
            Err(VerifyError::Modulus {
                found: 1,
                composite: false
            })
        );
        assert_eq!(check_modulus(7), Ok(7));
    }

    #[test]
    fn nvars_above_the_range_is_rejected() {
        assert_eq!(check_nvars(256), Ok(256));
        assert_eq!(
            check_nvars(257),
            Err(VerifyError::Nvars {
                found: 257,
                max: 256
            })
        );
    }

    #[test]
    fn a_reducible_tail_is_rejected() {
        let good = vec![poly(&[(1, &[2, 0]), (6, &[0, 1])]), poly(&[(1, &[1, 1])])];
        assert_eq!(check_basis(&good, &mut budget()), Ok(()));
        let bad = vec![
            poly(&[(1, &[2, 0]), (6, &[0, 1])]),
            poly(&[(1, &[1, 1]), (6, &[0, 2])]),
            poly(&[(1, &[0, 2])]),
        ];
        assert_eq!(
            check_basis(&bad, &mut budget()),
            Err(VerifyError::Basis {
                index: 1,
                fault: BasisFault::TailReducible { term: 1, by: 2 }
            })
        );
    }

    #[test]
    fn a_redundant_leading_monomial_is_rejected() {
        let basis = vec![poly(&[(1, &[2, 1])]), poly(&[(1, &[1, 1])])];
        assert_eq!(
            check_basis(&basis, &mut budget()),
            Err(VerifyError::Basis {
                index: 0,
                fault: BasisFault::LeadDivisible { by: 1 }
            })
        );
    }

    #[test]
    fn an_unsorted_basis_is_rejected() {
        let basis = vec![poly(&[(1, &[1, 1])]), poly(&[(1, &[2, 0])])];
        assert_eq!(
            check_basis(&basis, &mut budget()),
            Err(VerifyError::Basis {
                index: 1,
                fault: BasisFault::NotDescending
            })
        );
    }

    #[test]
    fn a_non_monic_element_is_rejected() {
        let basis = vec![poly(&[(2, &[1, 0])])];
        assert_eq!(
            check_basis(&basis, &mut budget()),
            Err(VerifyError::Basis {
                index: 0,
                fault: BasisFault::NotMonic { lc: 2 }
            })
        );
    }

    #[test]
    fn a_duplicate_monomial_is_rejected() {
        let repeated = poly(&[(1, &[1, 0]), (2, &[1, 0])]);
        assert_eq!(
            check_poly(&repeated, Location::Input(0), 2, 7, &mut budget()),
            Err(VerifyError::Polynomial {
                at: Location::Input(0),
                fault: PolyFault::DuplicateMonomial { term: 1 }
            })
        );
    }
}
