//! Cofactor vectors over the canonical input.
//!
//! An origin holds one cofactor per input polynomial. For a value `g` the
//! engines carry, the origin `c` states the identity
//!
//! ```text
//! g = sum_i c_i * f_i
//! ```
//!
//! over the canonicalized input F. The engines keep the identity through
//! every step they take, so [`super::assemble`] writes origins that the
//! engine derived, never origins recomputed after the fact. Every function
//! here updates an origin by the same scalar and monomial the polynomial
//! step uses.

use crate::poly::{ExponentOverflow, Monomial, Polynomial};
use crate::ring::PolynomialRing;
use crate::ring::field::Felt;

/// One cofactor per input polynomial.
pub(crate) type Origin = Vec<Polynomial>;

/// The origin of input polynomial `index` in a list of `count` inputs.
///
/// Cofactor `index` is the constant 1; the rest are zero.
pub(crate) fn unit(ring: &PolynomialRing, index: usize, count: usize) -> Origin {
    let mut origin = vec![ring.zero(); count];
    origin[index] = ring.one();
    origin
}

/// Apply `target -= coeff * mono * source`, cofactor by cofactor.
///
/// This is the origin step that matches subtracting `coeff * mono * g`
/// from a work polynomial, where `source` is the origin of `g`.
///
/// A cofactor's degree grows with every step and no polynomial lead bounds
/// it, so a product can leave the width one exponent holds. That is
/// [`ExponentOverflow`], and the engine reports it as
/// [`crate::compute::ComputeError::DegreeLimit`]. The origin is left
/// partly updated, so the caller must drop it.
pub(crate) fn sub_scaled(
    target: &mut Origin,
    source: &Origin,
    coeff: Felt,
    mono: &Monomial,
    p: u64,
) -> Result<(), ExponentOverflow> {
    for (slot, cofactor) in target.iter_mut().zip(source) {
        if cofactor.is_zero() {
            continue;
        }
        *slot = slot.sub_scaled(cofactor, coeff, mono, p)?;
    }
    Ok(())
}

/// Return `left_coeff * left_mono * left - right_coeff * right_mono * right`.
///
/// This is the origin step that matches forming an S-polynomial from two
/// term-multiplied parents. A product past the width one exponent holds is
/// [`ExponentOverflow`], as in [`sub_scaled`].
pub(crate) fn combine(
    left: &Origin,
    left_coeff: Felt,
    left_mono: &Monomial,
    right: &Origin,
    right_coeff: Felt,
    right_mono: &Monomial,
    p: u64,
) -> Result<Origin, ExponentOverflow> {
    left.iter()
        .zip(right)
        .map(|(a, b)| {
            let scaled = a.scale_monomial(left_coeff, left_mono, p)?;
            Ok(scaled.sub(&b.scale_monomial(right_coeff, right_mono, p)?, p))
        })
        .collect()
}

/// Multiply every cofactor by `coeff`.
///
/// `coeff` must be nonzero, so no term drops out and the term order holds.
pub(crate) fn scale(target: &mut Origin, coeff: Felt, p: u64) {
    debug_assert!(
        !coeff.is_zero(),
        "scaling an origin by zero loses the identity"
    );
    for cofactor in target.iter_mut() {
        for term in &mut cofactor.terms {
            term.coeff = term.coeff.mul(coeff, p);
        }
    }
}

/// Make `poly` monic and scale `target` by the same factor.
///
/// Returns the monic polynomial. The zero polynomial has no leading
/// coefficient, so it leaves the origin alone.
pub(crate) fn make_monic(poly: &Polynomial, target: &mut Origin, p: u64) -> Polynomial {
    if let Some(lc) = poly.lc() {
        scale(target, lc.inv(p), p);
    }
    poly.make_monic(p)
}

/// Report whether `poly` equals `sum_i c_i * f_i` over the input.
///
/// The engines call this from a `debug_assert`, so the check runs in debug
/// builds only.
pub(crate) fn holds(origin: &Origin, input: &[Polynomial], poly: &Polynomial, p: u64) -> bool {
    if origin.len() != input.len() {
        return false;
    }
    let mut sum = poly.zero_like();
    for (cofactor, f) in origin.iter().zip(input) {
        for term in &cofactor.terms {
            let Ok(scaled) = f.scale_monomial(term.coeff, &term.mono, p) else {
                return false;
            };
            sum = sum.add(&scaled, p);
        }
    }
    sum == *poly
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring() -> PolynomialRing {
        PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime")
    }

    fn list(ring: &PolynomialRing, texts: &[&str]) -> Vec<Polynomial> {
        texts
            .iter()
            .map(|text| ring.parse_polynomial(text).expect("the text parses"))
            .collect()
    }

    #[test]
    fn a_unit_origin_reproduces_its_input_polynomial() {
        let ring = ring();
        let input = list(&ring, &["x", "y + 3"]);
        for index in 0..input.len() {
            let origin = unit(&ring, index, input.len());
            assert!(holds(&origin, &input, &input[index], 7));
        }
    }

    #[test]
    fn a_subtraction_step_keeps_the_identity() {
        let ring = ring();
        let input = list(&ring, &["x", "y"]);
        let mono = ring
            .parse_polynomial("x*y^2")
            .expect("parses")
            .lm()
            .expect("non-zero")
            .clone();
        let coeff = Felt::new(3, 7);

        let expected = input[0]
            .sub_scaled(&input[1], coeff, &mono, 7)
            .expect("the product fits");
        let mut origin = unit(&ring, 0, 2);
        sub_scaled(&mut origin, &unit(&ring, 1, 2), coeff, &mono, 7).expect("the product fits");
        assert!(holds(&origin, &input, &expected, 7));
    }

    #[test]
    fn a_combination_step_keeps_the_identity() {
        let ring = ring();
        let input = list(&ring, &["x^2", "x*y"]);
        let left_mono = ring
            .parse_polynomial("y")
            .expect("parses")
            .lm()
            .expect("non-zero")
            .clone();
        let right_mono = ring
            .parse_polynomial("x")
            .expect("parses")
            .lm()
            .expect("non-zero")
            .clone();
        let (left_coeff, right_coeff) = (Felt::new(2, 7), Felt::new(5, 7));

        let expected = input[0]
            .scale_monomial(left_coeff, &left_mono, 7)
            .expect("the product fits")
            .sub(
                &input[1]
                    .scale_monomial(right_coeff, &right_mono, 7)
                    .expect("the product fits"),
                7,
            );
        let origin = combine(
            &unit(&ring, 0, 2),
            left_coeff,
            &left_mono,
            &unit(&ring, 1, 2),
            right_coeff,
            &right_mono,
            7,
        )
        .expect("the product fits");
        assert!(holds(&origin, &input, &expected, 7));
    }

    #[test]
    fn making_a_polynomial_monic_keeps_the_identity() {
        let ring = ring();
        let input = list(&ring, &["3*x + 2"]);
        let mut origin = unit(&ring, 0, 1);
        let monic = make_monic(&input[0], &mut origin, 7);
        assert_eq!(monic.lc(), Some(Felt::one()));
        assert!(holds(&origin, &input, &monic, 7));
    }

    #[test]
    fn a_wrong_cofactor_fails_the_identity() {
        let ring = ring();
        let input = list(&ring, &["x", "y"]);
        assert!(!holds(&unit(&ring, 0, 2), &input, &input[1], 7));
    }
}
