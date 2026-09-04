//! Rational reconstruction: one residue, one modulus, one fraction.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::compute::{ComputeLimits, RunError};

/// Lift `residue` modulo `modulus` to a rational number, or report no
/// lift.
///
/// The result is the reduced fraction `a / b` with `b > 0`,
/// `gcd(a, b) = 1`, `gcd(b, modulus) = 1`, `a = r b` modulo `modulus`,
/// `|a| <= A`, and `b <= B`, where `A = B = floor(sqrt((M - 1) / 2))`. Then
/// `2 A B < M`, and one reduced fraction at most meets every condition, so
/// the answer is unique when it exists.
///
/// The search is the extended Euclidean algorithm on `(modulus, residue)`,
/// stopped at the first remainder at or below `A`. The result is then
/// checked against every condition above, so a pair the search produces
/// and the conditions reject is reported as no lift rather than returned.
///
/// A lift that succeeds is the right answer only if `modulus` is large
/// enough for the true value. Nothing here can know that.
///
/// The Euclidean loop runs on integers the size of the modulus, so
/// `limits` stops it. No lift is `Ok(None)`. A stopped loop is `Err`.
pub(crate) fn reconstruct(
    residue: &BigInt,
    modulus: &BigInt,
    limits: &ComputeLimits,
) -> Result<Option<BigRational>, RunError> {
    debug_assert!(
        !residue.is_negative() && residue < modulus,
        "the residue is reduced"
    );
    let bound = ((modulus - 1u32) / 2u32).sqrt();
    if bound.is_zero() {
        return Ok(None);
    }

    let (numerator, denominator) = euclidean_candidate(residue, modulus, &bound, limits)?;
    Ok(
        valid_candidate(&numerator, &denominator, residue, modulus, &bound)
            .then(|| BigRational::new_raw(numerator, denominator)),
    )
}

fn euclidean_candidate(
    residue: &BigInt,
    modulus: &BigInt,
    bound: &BigInt,
    limits: &ComputeLimits,
) -> Result<(BigInt, BigInt), RunError> {
    let (mut r0, mut r1) = (modulus.clone(), residue.clone());
    let (mut t0, mut t1) = (BigInt::zero(), BigInt::one());
    let mut step = 0usize;
    while r1 > *bound {
        if let Some(stop) = limits.stop_every(step) {
            return Err(stop);
        }
        step += 1;
        let quotient = &r0 / &r1;
        let r = r0 - &quotient * &r1;
        r0 = std::mem::replace(&mut r1, r);
        let t = t0 - &quotient * &t1;
        t0 = std::mem::replace(&mut t1, t);
    }

    Ok(if t1.is_negative() {
        (-r1, -t1)
    } else {
        (r1, t1)
    })
}

fn valid_candidate(
    numerator: &BigInt,
    denominator: &BigInt,
    residue: &BigInt,
    modulus: &BigInt,
    bound: &BigInt,
) -> bool {
    denominator.is_positive()
        && denominator <= bound
        && &numerator.abs() <= bound
        && numerator.gcd(denominator).is_one()
        && denominator.gcd(modulus).is_one()
        && (numerator - residue * denominator) % modulus == BigInt::zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Budget;
    use proptest::prelude::*;
    use std::time::Duration;

    /// The lift under no limit.
    fn lift(residue: &BigInt, modulus: &BigInt) -> Option<BigRational> {
        reconstruct(residue, modulus, &ComputeLimits::default()).expect("no limit stops it")
    }

    fn integer(value: i64) -> BigInt {
        BigInt::from(value)
    }

    /// The residue of `numerator / denominator` modulo `modulus`.
    fn residue_of(numerator: i64, denominator: i64, modulus: &BigInt) -> BigInt {
        let inverse = mod_inverse(&integer(denominator), modulus);
        let value = (integer(numerator) * inverse) % modulus;
        if value.is_negative() {
            value + modulus
        } else {
            value
        }
    }

    /// The inverse of `value` modulo `modulus`, by the extended Euclidean
    /// algorithm.
    fn mod_inverse(value: &BigInt, modulus: &BigInt) -> BigInt {
        let extended = value.extended_gcd(modulus);
        assert!(extended.gcd.is_one(), "the value is invertible");
        ((extended.x % modulus) + modulus) % modulus
    }

    #[test]
    fn a_residue_of_zero_lifts_to_zero() {
        let modulus = integer(1_000_003);
        let lifted = lift(&BigInt::zero(), &modulus).expect("zero lifts");
        assert_eq!(lifted, BigRational::new_raw(BigInt::zero(), BigInt::one()));
    }

    #[test]
    fn an_integer_lifts_with_a_denominator_of_one() {
        let modulus = integer(1_000_003);
        let lifted = lift(&integer(17), &modulus).expect("17 lifts");
        assert_eq!(lifted.numer(), &integer(17));
        assert_eq!(lifted.denom(), &integer(1));
    }

    #[test]
    fn a_negative_numerator_lifts_as_itself() {
        let modulus = integer(1_000_003);
        let residue = residue_of(-5, 3, &modulus);
        let lifted = lift(&residue, &modulus).expect("-5/3 lifts");
        assert_eq!(lifted.numer(), &integer(-5));
        assert_eq!(lifted.denom(), &integer(3));
    }

    #[test]
    fn a_value_outside_the_bounds_does_not_lift_to_itself() {
        // 2 A B < M fails for this fraction, so the true value is outside
        // what the modulus supports.
        let modulus = integer(101);
        let residue = residue_of(37, 41, &modulus);
        let bound = ((&modulus - 1u32) / 2u32).sqrt();
        match lift(&residue, &modulus) {
            None => {}
            Some(lifted) => {
                assert_ne!(lifted.numer(), &integer(37));
                assert!(lifted.numer().abs() <= bound);
                assert!(lifted.denom() <= &bound);
            }
        }
    }

    /// The Euclidean loop reads the deadline, so an exhausted one stops
    /// it.
    #[test]
    fn a_lift_stops_on_an_exhausted_deadline() {
        let modulus = integer(1_000_003);
        let residue = residue_of(-5, 3, &modulus);
        let limits = ComputeLimits::of_budget(&Budget::new().timeout(Duration::ZERO));
        assert_eq!(
            reconstruct(&residue, &modulus, &limits).map(|lifted| lifted.is_some()),
            Err(RunError::Compute(crate::ComputeError::Timeout))
        );
    }

    #[test]
    fn a_modulus_below_three_lifts_nothing() {
        assert_eq!(lift(&BigInt::zero(), &BigInt::one()), None);
        assert_eq!(lift(&BigInt::zero(), &integer(2)), None);
    }

    proptest! {
        /// Whatever the input, the result meets every condition of the
        /// contract.
        #[test]
        fn a_lift_meets_the_contract(residue in 0i64..1_000_003, modulus in 3i64..1_000_003) {
            let modulus = integer(modulus);
            let residue = integer(residue) % &modulus;
            let Some(lifted) = lift(&residue, &modulus) else { return Ok(()); };
            let bound = ((&modulus - 1u32) / 2u32).sqrt();
            prop_assert!(lifted.denom().is_positive());
            prop_assert!(lifted.denom() <= &bound);
            prop_assert!(lifted.numer().abs() <= bound);
            prop_assert!(lifted.numer().gcd(lifted.denom()).is_one());
            prop_assert!(lifted.denom().gcd(&modulus).is_one());
            prop_assert!(
                (lifted.numer() - &residue * lifted.denom()) % &modulus == BigInt::zero()
            );
        }

        /// A fraction inside the bounds is the one the lift returns.
        #[test]
        fn a_fraction_inside_the_bounds_comes_back(
            numerator in -30i64..=30,
            denominator in 1i64..=30,
        ) {
            let modulus = integer(1_000_003);
            let reduced = BigRational::new(integer(numerator), integer(denominator));
            let residue = residue_of(
                numerator,
                denominator,
                &modulus,
            );
            prop_assert_eq!(lift(&residue, &modulus), Some(reduced));
        }
    }
}
