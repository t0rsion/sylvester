//! Interreduction of a Gröbner basis to the reduced basis.

use std::time::Instant;

use super::{ComputeError, REDUCE_DEADLINE_STRIDE, poll_deadline};
use crate::cert::origin::{self, Origin};
use crate::poly::Polynomial;
use crate::ring::PolynomialRing;

/// Stop the interreduction when the budget runs out.
///
/// The final interreduction can cost more than the pair loop that feeds it,
/// so it carries the same typed partiality as the engines.
fn check_limits(
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
    basis: &[Polynomial],
    origins: &[Origin],
    nvars: usize,
) -> Result<(), ComputeError> {
    poll_deadline(deadline)?;

    if let Some(limit) = max_memory_bytes {
        let mut bytes = std::mem::size_of_val(basis);
        bytes = bytes.saturating_add(term_bytes(basis, nvars));
        bytes = bytes.saturating_add(std::mem::size_of_val(origins));
        for origin in origins {
            bytes = bytes.saturating_add(term_bytes(origin, nvars));
        }
        if bytes > limit {
            return Err(ComputeError::MemoryLimitExceeded);
        }
    }

    Ok(())
}

/// The bytes the terms of `polys` occupy over a ring of `nvars` variables.
fn term_bytes(polys: &[Polynomial], nvars: usize) -> usize {
    let per_term = super::per_term_bytes(nvars);
    polys.iter().map(|poly| poly.terms.len() * per_term).sum()
}

fn minimalize_leading_terms(mut basis: Vec<Polynomial>) -> Vec<Polynomial> {
    basis.retain(|f| !f.is_zero());
    if basis.len() <= 1 {
        return basis;
    }

    basis.sort_by(|a, b| a.lm().unwrap().cmp(b.lm().unwrap()));
    let mut minimal: Vec<Polynomial> = Vec::with_capacity(basis.len());
    for f in basis {
        let Some(lm_f) = f.lm() else { continue };
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
    reduced_groebner_basis_checked(ring, basis, None, None)
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
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
) -> Result<Vec<Polynomial>, ComputeError> {
    let modulus = ring.modulus();
    let nvars = ring.nvars();
    basis.retain(|f| !f.is_zero());
    basis = basis.into_iter().map(|f| f.make_monic(modulus)).collect();

    basis = minimalize_leading_terms(basis);

    loop {
        check_limits(deadline, max_memory_bytes, &basis, &[], nvars)?;
        let mut next: Vec<Polynomial> = Vec::with_capacity(basis.len());
        for i in 0..basis.len() {
            check_limits(deadline, max_memory_bytes, &basis, &[], nvars)?;
            let reducers: Vec<&Polynomial> = basis
                .iter()
                .enumerate()
                .filter_map(|(j, g)| (i != j).then_some(g))
                .collect();
            let r = basis[i].normal_form_refs(&reducers, modulus)?;
            if !r.is_zero() {
                next.push(r.make_monic(modulus));
            }
        }

        next = minimalize_leading_terms(next);
        next.sort_by(|a, b| b.lm().unwrap().cmp(a.lm().unwrap()));

        let mut current = basis;
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
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
) -> Result<(Vec<Polynomial>, Vec<Origin>), ComputeError> {
    let modulus = ring.modulus();
    let nvars = ring.nvars();
    let (mut basis, mut origins) = {
        let kept: Vec<(Polynomial, Origin)> = basis
            .into_iter()
            .zip(origins)
            .filter(|(f, _)| !f.is_zero())
            .map(|(f, mut cofactors)| {
                let monic = origin::make_monic(&f, &mut cofactors, modulus);
                (monic, cofactors)
            })
            .collect();
        minimalize_leading_terms_tracked(kept)
    };
    debug_assert!(
        holds_for_all(&basis, &origins, input, modulus),
        "interreduction must keep the origin identity"
    );

    loop {
        check_limits(deadline, max_memory_bytes, &basis, &origins, nvars)?;
        let mut next: Vec<(Polynomial, Origin)> = Vec::with_capacity(basis.len());
        for i in 0..basis.len() {
            check_limits(deadline, max_memory_bytes, &basis, &origins, nvars)?;
            let (r, mut cofactors) = normal_form_excluding(&basis, &origins, i, modulus, deadline)?;
            if !r.is_zero() {
                let monic = origin::make_monic(&r, &mut cofactors, modulus);
                debug_assert!(
                    origin::holds(&cofactors, input, &monic, modulus),
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

/// Drop the elements whose leading monomial another element divides,
/// keeping the cofactors with them.
///
/// This repeats the choices of [`minimalize_leading_terms`]: the same
/// stable sort by leading monomial, the same first-wins scan.
fn minimalize_leading_terms_tracked(
    mut basis: Vec<(Polynomial, Origin)>,
) -> (Vec<Polynomial>, Vec<Origin>) {
    basis.retain(|(f, _)| !f.is_zero());
    if basis.len() > 1 {
        basis.sort_by(|(a, _), (b, _)| a.lm().unwrap().cmp(b.lm().unwrap()));
        let mut minimal: Vec<(Polynomial, Origin)> = Vec::with_capacity(basis.len());
        for (f, cofactors) in basis {
            let Some(lm_f) = f.lm() else { continue };
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
    modulus: u64,
    deadline: Option<Instant>,
) -> Result<(Polynomial, Origin), ComputeError> {
    let mut poly = basis[target].clone();
    let mut cofactors = origins[target].clone();
    let mut remainder = poly.zero_like();
    let mut since_check = 0usize;

    while let Some(lt_p) = poly.lt().cloned() {
        since_check += poly.terms.len();
        if since_check >= REDUCE_DEADLINE_STRIDE {
            since_check = 0;
            poll_deadline(deadline)?;
        }
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
                    .expect("divides() implies quotient()");
                let scale = lt_p.coeff.div(lt_g.coeff, modulus);
                poly = poly.sub_scaled(g, scale, &m, modulus)?;
                origin::sub_scaled(&mut cofactors, &origins[index], scale, &m, modulus)?;
                reduced = true;
                break;
            }
        }

        if !reduced {
            let term = poly.pop_lt().expect("lt_p came from this polynomial");
            remainder.push_term(term, modulus);
        }
    }

    Ok((remainder, cofactors))
}

fn holds_for_all(
    basis: &[Polynomial],
    origins: &[Origin],
    input: &[Polynomial],
    modulus: u64,
) -> bool {
    basis.len() == origins.len()
        && basis
            .iter()
            .zip(origins)
            .all(|(poly, cofactors)| origin::holds(cofactors, input, poly, modulus))
}

/// Report whether every S-polynomial of the basis reduces to zero over it.
#[cfg(test)]
pub(crate) fn is_groebner_basis(basis: &[Polynomial], modulus: u64) -> bool {
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            let s = basis[i]
                .s_polynomial(&basis[j], modulus)
                .expect("a basis over a checked input multiplies inside the width");
            if s.is_zero() {
                continue;
            }
            let reducers: Vec<&Polynomial> = basis.iter().collect();
            let r = s
                .normal_form_refs(&reducers, modulus)
                .expect("a basis over a checked input multiplies inside the width");
            if !r.is_zero() {
                return false;
            }
        }
    }
    true
}

/// Report whether the basis is the reduced Gröbner basis of its ideal.
#[cfg(test)]
pub(crate) fn is_reduced_basis(basis: &[Polynomial], modulus: u64) -> bool {
    for i in 0..basis.len() {
        let Some(lm_i) = basis[i].lm() else { continue };
        if basis[i].lc() != Some(crate::ring::field::Felt::one()) {
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

    is_groebner_basis(basis, modulus)
}

#[cfg(test)]
mod tests {
    use super::{is_groebner_basis, is_reduced_basis, reduced_groebner_basis};
    use crate::compute::classic;
    use crate::poly::Polynomial;
    use crate::ring::PolynomialRing;
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
        assert_eq!(reduced[0].leading_term(), Some((1, [1u16].as_slice())));
        assert!(is_groebner_basis(&reduced, 7));
    }

    #[test]
    fn the_classic_backend_reduces_a_small_system() {
        let ring = PolynomialRing::prime_field(32003, ["x", "y"]).expect("32003 is prime");
        let gb = classic::solve(&ring, &system(&ring, &["x*y - 1", "y^2 - y"]));
        assert!(is_reduced_basis(&gb, 32003));
    }

    #[test]
    fn the_classic_backend_reduces_cyclic_3() {
        let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"]).expect("32003 is prime");
        let gb = classic::solve(
            &ring,
            &system(&ring, &["x + y + z", "x*y + y*z + z*x", "x*y*z - 1"]),
        );
        assert!(is_reduced_basis(&gb, 32003));
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
        assert!(is_reduced_basis(&gb, 32003));
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
        assert!(is_reduced_basis(&gb, 32003));
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
        let expired = Instant::now() - Duration::from_secs(1);

        assert_eq!(
            super::reduced_groebner_basis_checked(&ring, basis.clone(), Some(expired), None),
            Err(super::ComputeError::Timeout),
            "the untracked interreduction must report the exhausted budget"
        );

        // Reading the basis as its own input makes the unit cofactors the
        // true origins, so the identity holds going in.
        let origins: Vec<crate::cert::origin::Origin> = (0..basis.len())
            .map(|index| crate::cert::origin::unit(&ring, index, basis.len()))
            .collect();
        assert_eq!(
            super::reduced_groebner_basis_tracked(
                &ring,
                basis.clone(),
                origins,
                &basis,
                Some(expired),
                None
            ),
            Err(super::ComputeError::Timeout),
            "the tracked interreduction must report the exhausted budget"
        );
    }
}
