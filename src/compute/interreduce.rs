//! Interreduction of a Gröbner basis to the reduced basis.

use super::{ComputeError, ComputeLimits, RunError};
use crate::cert::origin::{self, Origin};
use crate::poly::{ExponentOverflow, Polynomial};
use crate::ring::{PolynomialRing, PrimeOps};

/// Stop the interreduction when the budget runs out.
///
/// The final interreduction can cost more than the pair loop that feeds
/// it, so it carries the same typed partiality as the engines. Both
/// budgets off costs two comparisons per call.
fn check_limits(
    limits: &ComputeLimits,
    basis: &[Polynomial],
    origins: &[Origin],
) -> Result<(), RunError> {
    if let Some(stop) = limits.stop() {
        return Err(stop);
    }

    if let Some(limit) = limits.memory {
        let mut bytes = std::mem::size_of_val(basis);
        bytes = bytes.saturating_add(term_bytes(basis));
        bytes = bytes.saturating_add(std::mem::size_of_val(origins));
        for origin in origins {
            bytes = bytes.saturating_add(term_bytes(origin));
        }
        if bytes > limit {
            return Err(RunError::Compute(ComputeError::MemoryLimitExceeded));
        }
    }

    Ok(())
}

fn term_bytes(polys: &[Polynomial]) -> usize {
    polys.iter().map(|poly| poly.heap_bytes()).sum()
}

fn minimalize_leading_terms(mut basis: Vec<Polynomial>) -> Vec<Polynomial> {
    basis.retain(|f| !f.is_zero());
    if basis.len() <= 1 {
        return basis;
    }

    // zero polynomials were filtered, so lm exists for all entries.
    basis.sort_by(|a, b| a.lm().unwrap().cmp(b.lm().unwrap()));
    let mut minimal: Vec<Polynomial> = Vec::with_capacity(basis.len());
    for f in basis {
        let Some(lm_f) = f.lm() else { continue };
        // minimal contains only non-zero polynomials.
        if minimal.iter().any(|g| g.lm().unwrap().divides(lm_f)) {
            continue;
        }
        minimal.push(f);
    }
    minimal
}

/// Interreduce a Gröbner basis to the reduced basis, with no budget.
///
/// The input must already be a Gröbner basis: this interreduces, it does
/// not complete.
#[cfg(test)]
pub(crate) fn reduced_groebner_basis(
    ring: &PolynomialRing,
    basis: Vec<Polynomial>,
) -> Vec<Polynomial> {
    // without a budget the only error paths cannot fire.
    reduced_groebner_basis_checked(ring, basis, &ComputeLimits::default())
        .expect("interreduction without a budget cannot stop early")
}

/// Interreduce under a deadline and a memory cap.
///
/// Reports [`ComputeError::Timeout`] or
/// [`ComputeError::MemoryLimitExceeded`] when a budget runs out here rather
/// than in the pair loop.
pub(crate) fn reduced_groebner_basis_checked(
    ring: &PolynomialRing,
    mut basis: Vec<Polynomial>,
    limits: &ComputeLimits,
) -> Result<Vec<Polynomial>, RunError> {
    let ops = ring.ops();
    basis.retain(|f| !f.is_zero());
    basis = basis.into_iter().map(|f| f.make_monic(ops)).collect();

    basis = minimalize_leading_terms(basis);

    loop {
        check_limits(limits, &basis, &[])?;
        let mut next: Vec<Polynomial> = Vec::with_capacity(basis.len());
        for i in 0..basis.len() {
            check_limits(limits, &basis, &[])?;
            let reducers: Vec<&Polynomial> = basis
                .iter()
                .enumerate()
                .filter_map(|(j, g)| (i != j).then_some(g))
                .collect();
            let r = basis[i].normal_form_refs(&reducers, ops);
            if !r.is_zero() {
                next.push(r.make_monic(ops));
            }
        }

        next = minimalize_leading_terms(next);
        // zero polynomials were filtered, so lm exists for all entries.
        next.sort_by(|a, b| b.lm().unwrap().cmp(a.lm().unwrap()));

        let mut current = basis;
        // zero polynomials were filtered, so lm exists for all entries.
        current.sort_by(|a, b| b.lm().unwrap().cmp(a.lm().unwrap()));

        if next == current {
            return Ok(next);
        }
        basis = next;
    }
}

/// Interreduce and carry the origin cofactors through every step.
///
/// `basis` and `origins` line up: `origins[j]` names the combination of
/// `input` that gives `basis[j]`. The result keeps that identity for the
/// reduced basis. `input` is read only by the debug assertions.
///
/// This mirrors [`reduced_groebner_basis_checked`] step for step, so both
/// return the same basis for the same argument.
pub(crate) fn reduced_groebner_basis_tracked(
    ring: &PolynomialRing,
    basis: Vec<Polynomial>,
    origins: Vec<Origin>,
    input: &[Polynomial],
    limits: &ComputeLimits,
) -> Result<(Vec<Polynomial>, Vec<Origin>), RunError> {
    let ops = ring.ops();
    let (mut basis, mut origins) = {
        let kept: Vec<(Polynomial, Origin)> = basis
            .into_iter()
            .zip(origins)
            .filter(|(f, _)| !f.is_zero())
            .map(|(f, mut cofactors)| {
                let monic = origin::make_monic(&f, &mut cofactors, ops);
                (monic, cofactors)
            })
            .collect();
        minimalize_leading_terms_tracked(kept)
    };
    debug_assert!(
        holds_for_all(&basis, &origins, input, ops),
        "interreduction must keep the origin identity"
    );

    loop {
        check_limits(limits, &basis, &origins)?;
        let mut next: Vec<(Polynomial, Origin)> = Vec::with_capacity(basis.len());
        for i in 0..basis.len() {
            check_limits(limits, &basis, &origins)?;
            let (r, mut cofactors) = normal_form_excluding(&basis, &origins, i, ops)?;
            if !r.is_zero() {
                let monic = origin::make_monic(&r, &mut cofactors, ops);
                debug_assert!(
                    origin::holds(&cofactors, input, &monic, ops),
                    "interreduction must keep the origin identity"
                );
                next.push((monic, cofactors));
            }
        }

        let (mut next, mut next_origins) = minimalize_leading_terms_tracked(next);
        sort_descending(&mut next, &mut next_origins);
        sort_descending(&mut basis, &mut origins);

        if next == basis {
            return Ok((next, next_origins));
        }
        basis = next;
        origins = next_origins;
    }
}

/// Drop the elements whose leading monomial another element's divides,
/// keeping the cofactors with them.
///
/// This repeats the choices of [`minimalize_leading_terms`]: the same
/// stable sort by leading monomial, the same first-wins scan.
fn minimalize_leading_terms_tracked(
    mut basis: Vec<(Polynomial, Origin)>,
) -> (Vec<Polynomial>, Vec<Origin>) {
    basis.retain(|(f, _)| !f.is_zero());
    if basis.len() > 1 {
        // zero polynomials were filtered, so lm exists for all entries.
        basis.sort_by(|(a, _), (b, _)| a.lm().unwrap().cmp(b.lm().unwrap()));
        let mut minimal: Vec<(Polynomial, Origin)> = Vec::with_capacity(basis.len());
        for (f, cofactors) in basis {
            let Some(lm_f) = f.lm() else { continue };
            // minimal holds only non-zero polynomials.
            if minimal.iter().any(|(g, _)| g.lm().unwrap().divides(lm_f)) {
                continue;
            }
            minimal.push((f, cofactors));
        }
        basis = minimal;
    }
    basis.into_iter().unzip()
}

/// Sort by leading monomial, largest first, keeping the cofactors aligned.
fn sort_descending(basis: &mut Vec<Polynomial>, origins: &mut Vec<Origin>) {
    let mut pairs: Vec<(Polynomial, Origin)> = basis.drain(..).zip(origins.drain(..)).collect();
    // zero polynomials were filtered, so lm exists for all entries.
    pairs.sort_by(|(a, _), (b, _)| b.lm().unwrap().cmp(a.lm().unwrap()));
    for (poly, cofactors) in pairs {
        basis.push(poly);
        origins.push(cofactors);
    }
}

/// Reduce `basis[target]` by every other element and carry the cofactors.
///
/// This repeats the reducer choice of `Polynomial::normal_form_refs` over
/// the same elements in the same order, so it returns the same remainder.
fn normal_form_excluding(
    basis: &[Polynomial],
    origins: &[Origin],
    target: usize,
    ops: &PrimeOps,
) -> Result<(Polynomial, Origin), ExponentOverflow> {
    let mut poly = basis[target].clone();
    let mut cofactors = origins[target].clone();
    let mut remainder = poly.zero_like();

    while let Some(lt_p) = poly.lt().cloned() {
        let mut reduced = false;

        for (index, g) in basis.iter().enumerate() {
            if index == target {
                continue;
            }
            let Some(lt_g) = g.lt() else { continue };
            if lt_g.mono.divides(&lt_p.mono) {
                let m = lt_p
                    .mono
                    .quotient(&lt_g.mono)
                    // divides() implies a quotient exists.
                    .expect("divides() implies quotient()");
                let scale = lt_p.coeff.div(lt_g.coeff, ops.modulus());
                poly = poly.sub_scaled_checked(g, &scale, &m, ops)?;
                origin::sub_scaled(&mut cofactors, &origins[index], scale, &m, ops)?;
                reduced = true;
                break;
            }
        }

        if !reduced {
            let term = poly
                .pop_lt()
                // lt_p was Some, so poly is non-empty here.
                .expect("polynomial should not be empty");
            remainder.push_term(term, ops);
        }
    }

    Ok((remainder, cofactors))
}

fn holds_for_all(
    basis: &[Polynomial],
    origins: &[Origin],
    input: &[Polynomial],
    ops: &PrimeOps,
) -> bool {
    basis.len() == origins.len()
        && basis
            .iter()
            .zip(origins)
            .all(|(poly, cofactors)| origin::holds(cofactors, input, poly, ops))
}

/// Report whether every S-polynomial of the basis reduces to zero over it.
#[cfg(test)]
pub(crate) fn is_groebner_basis(basis: &[Polynomial], ops: &PrimeOps) -> bool {
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            let s = basis[i].s_polynomial(&basis[j], ops);
            if s.is_zero() {
                continue;
            }
            let reducers: Vec<&Polynomial> = basis.iter().collect();
            let r = s.normal_form_refs(&reducers, ops);
            if !r.is_zero() {
                return false;
            }
        }
    }
    true
}

/// Report whether the basis is the reduced Gröbner basis of its ideal.
#[cfg(test)]
pub(crate) fn is_reduced_basis(basis: &[Polynomial], ops: &PrimeOps) -> bool {
    for i in 0..basis.len() {
        let Some(lm_i) = basis[i].lm() else { continue };
        if basis[i].lc() != Some(&crate::ring::field::Felt::one()) {
            return false;
        }
        for j in 0..basis.len() {
            if i == j {
                continue;
            }
            let Some(lm_j) = basis[j].lm() else { continue };
            if lm_j.divides(lm_i) {
                return false;
            }
            for term in &basis[i].terms {
                if lm_j.divides(&term.mono) {
                    return false;
                }
            }
        }
    }

    is_groebner_basis(basis, ops)
}

#[cfg(test)]
mod tests {
    use super::{is_groebner_basis, is_reduced_basis, reduced_groebner_basis};
    use crate::compute::{ComputeLimits, RunError, classic};
    use crate::poly::Polynomial;
    use crate::ring::PolynomialRing;
    use std::mem::size_of;
    use std::time::{Duration, Instant};

    fn system(ring: &PolynomialRing, texts: &[&str]) -> Vec<Polynomial> {
        texts
            .iter()
            .map(|text| ring.parse_polynomial(text).expect("the text parses"))
            .collect()
    }

    #[test]
    fn interreduction_drops_a_redundant_leading_term() {
        let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
        let basis = system(&ring, &["x^2", "2*x"]);
        let reduced = reduced_groebner_basis(&ring, basis);
        assert_eq!(reduced.len(), 1);
        assert_eq!(
            reduced[0]
                .leading_term()
                .map(|(coeff, exps)| (coeff.value(), exps)),
            Some((1, [1u16].as_slice()))
        );
        assert!(is_groebner_basis(&reduced, ring.ops()));
    }

    #[test]
    fn the_classic_backend_reduces_a_small_system() {
        let ring = PolynomialRing::prime_field(32003, ["x", "y"]).expect("32003 is prime");
        let gb = classic::solve(&ring, &system(&ring, &["x*y - 1", "y^2 - y"]));
        assert!(is_reduced_basis(&gb, ring.ops()));
    }

    #[test]
    fn the_classic_backend_reduces_cyclic_3() {
        let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"]).expect("32003 is prime");
        let gb = classic::solve(
            &ring,
            &system(&ring, &["x + y + z", "x*y + y*z + z*x", "x*y*z - 1"]),
        );
        assert!(is_reduced_basis(&gb, ring.ops()));
    }

    #[test]
    fn the_classic_backend_reduces_cyclic_4() {
        let ring =
            PolynomialRing::prime_field(32003, ["x", "y", "z", "w"]).expect("32003 is prime");
        let gb = classic::solve(
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
        assert!(is_reduced_basis(&gb, ring.ops()));
    }

    #[test]
    fn the_classic_backend_reduces_katsura_4() {
        let ring =
            PolynomialRing::prime_field(32003, ["a", "b", "c", "d", "e"]).expect("32003 is prime");
        let gb = classic::solve(
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
        assert!(is_reduced_basis(&gb, ring.ops()));
    }

    #[test]
    fn an_exhausted_deadline_stops_the_interreduction() {
        let ring =
            PolynomialRing::prime_field(32003, ["x", "y", "z", "w"]).expect("32003 is prime");
        let basis = classic::solve(
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
        let expired = ComputeLimits {
            deadline: Some(Instant::now() - Duration::from_secs(1)),
            ..ComputeLimits::default()
        };

        assert_eq!(
            super::reduced_groebner_basis_checked(&ring, basis.clone(), &expired),
            Err(RunError::Compute(super::ComputeError::Timeout)),
            "the untracked interreduction must report the exhausted budget"
        );

        // Reading the basis as its own input makes the unit cofactors the
        // true origins, so the identity holds going in.
        let origins: Vec<crate::cert::origin::Origin> = (0..basis.len())
            .map(|index| crate::cert::origin::unit(&ring, index, basis.len()))
            .collect();
        assert_eq!(
            super::reduced_groebner_basis_tracked(&ring, basis.clone(), origins, &basis, &expired),
            Err(RunError::Compute(super::ComputeError::Timeout)),
            "the tracked interreduction must report the exhausted budget"
        );
    }

    /// A ring past the inline exponent width holds every exponent vector
    /// on the heap. The memory meter counts those bytes: a meter that
    /// read `size_of::<Term>()` alone missed one allocation per term.
    #[test]
    fn the_memory_meter_counts_the_spilled_exponents() {
        use crate::poly::{Term, heap_exps_bytes};

        let nvars = 16;
        let names: Vec<String> = (0..nvars).map(|index| format!("x{index}")).collect();
        let ring = PolynomialRing::prime_field(32003, names).expect("32003 is prime");
        let mut square = vec![0u16; nvars];
        square[0] = 2;
        let mut single = vec![0u16; nvars];
        single[1] = 1;
        let basis = vec![
            ring.polynomial([(1i64, square)])
                .expect("the exponents fit"),
            ring.polynomial([(1i64, single)])
                .expect("the exponents fit"),
        ];

        let terms: usize = basis.iter().map(|poly| poly.terms.len()).sum();
        let inline_only = std::mem::size_of_val(&basis[..]) + terms * size_of::<Term>();
        let spilled = terms * heap_exps_bytes(nvars);
        assert!(spilled > 0, "16 variables do not fit the inline width");

        assert_eq!(
            super::reduced_groebner_basis_checked(&ring, basis.clone(), &limit(inline_only)),
            Err(RunError::Compute(super::ComputeError::MemoryLimitExceeded)),
            "the spilled exponents do not fit a limit that counts the terms alone"
        );
        assert!(
            super::reduced_groebner_basis_checked(&ring, basis, &limit(inline_only + spilled))
                .is_ok(),
            "the same basis fits a limit that counts them"
        );
    }

    /// Limits of `bytes` memory and nothing else.
    fn limit(bytes: usize) -> ComputeLimits {
        ComputeLimits {
            memory: Some(bytes),
            ..ComputeLimits::default()
        }
    }
}
