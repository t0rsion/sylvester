//! Arithmetic in the field of rational numbers.
//!
//! Every value is a [`BigRational`] in lowest terms with a positive
//! denominator, the form `num_rational` maintains. [`RationalOps`] is that
//! arithmetic in the shape the shared polynomial code calls.

use std::fmt;
use std::mem::size_of;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use super::field::{Felt, residue};
use super::{
    Coefficient, Domain, DomainOps, PolynomialRing, PrimeField, Rationals, RingError, sealed,
};
use crate::compute::{ComputeLimits, RunError};
use crate::poly::{Monomial, Polynomial, Term, heap_exps_bytes};

/// The rational operations the shared polynomial code calls.
///
/// The type stores nothing, because one field of rational numbers serves
/// every ring. The ring holds one value of it and hands it to every
/// polynomial operation, in place of the modulus [`PrimeOps`] carries.
///
/// [`PrimeOps`]: crate::PrimeOps
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RationalOps(());

impl RationalOps {
    pub(crate) fn new() -> Self {
        RationalOps(())
    }
}

impl sealed::Sealed for RationalOps {}

impl DomainOps for RationalOps {
    type Coeff = BigRational;

    fn one(&self) -> BigRational {
        BigRational::one()
    }

    fn is_one(&self, a: &BigRational) -> bool {
        a.is_one()
    }

    fn add(&self, a: &BigRational, b: &BigRational) -> Option<BigRational> {
        let sum = a + b;
        (!sum.is_zero()).then_some(sum)
    }

    fn sub(&self, a: &BigRational, b: &BigRational) -> Option<BigRational> {
        let difference = a - b;
        (!difference.is_zero()).then_some(difference)
    }

    fn neg(&self, a: &BigRational) -> BigRational {
        -a
    }

    fn mul(&self, a: &BigRational, b: &BigRational) -> BigRational {
        a * b
    }

    /// The reciprocal. The caller passes a nonzero coefficient.
    fn inv(&self, a: &BigRational) -> BigRational {
        assert!(!a.is_zero(), "zero has no inverse in a field");
        a.recip()
    }

    fn convert(&self, value: &Coefficient) -> Result<Option<BigRational>, RingError> {
        let (numerator, denominator) = value.normalized()?;
        if numerator.is_zero() {
            return Ok(None);
        }
        // normalized() reduces the fraction and makes the denominator
        // positive, which is the form BigRational holds.
        Ok(Some(BigRational::new_raw(numerator, denominator)))
    }

    /// Write the value as `a/b`, or as `a` when the denominator is 1.
    ///
    /// A negative value writes a minus sign, which the parser does not
    /// read inside a coefficient. `Display` on a polynomial writes the sign
    /// as the term separator and passes the absolute value here.
    fn write(&self, a: &BigRational, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if a.denom().is_one() {
            write!(f, "{}", a.numer())
        } else {
            write!(f, "{}/{}", a.numer(), a.denom())
        }
    }

    fn is_negative(&self, a: &BigRational) -> bool {
        a.numer().is_negative()
    }

    /// The heap bytes one coefficient holds.
    ///
    /// The number is an estimate: `num-bigint` reports the used bits of an
    /// integer and not the capacity it allocated, so this rounds the used
    /// bits of the numerator and the denominator up to whole limbs and can
    /// read below what the allocator holds.
    fn heap_bytes(&self, a: &BigRational) -> usize {
        limb_bytes(a.numer()) + limb_bytes(a.denom())
    }
}

impl Domain for Rationals {
    type Coeff = BigRational;
    type Ops = RationalOps;
    type BasisMeta = RationalMeta;
}

/// The bytes the digits of one integer hold, rounded up to whole limbs.
pub(crate) fn limb_bytes(value: &BigInt) -> usize {
    const LIMB: usize = size_of::<u64>();
    let bits = usize::try_from(value.bits()).unwrap_or(usize::MAX);
    bits.div_ceil(LIMB * 8) * LIMB
}

/// A rational polynomial with its denominators cleared.
///
/// Every coefficient is an integer, the greatest common divisor of the
/// coefficients is 1, and the terms run ascending under grevlex, so the
/// last term is the leading term. Scaling a generator by a nonzero
/// rational leaves the ideal it generates alone, so the modular driver
/// keeps no record of the factor it divided out.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ClearedPolynomial {
    terms: Vec<(BigInt, Monomial)>,
}

impl ClearedPolynomial {
    /// Clear the denominators of `f` and divide out the content.
    ///
    /// The multiplier is the least common multiple of the denominators,
    /// and the divisor is the greatest common divisor of the integer
    /// coefficients that follow. The zero polynomial holds no term and
    /// clears to no term.
    ///
    /// Every loop here runs once per term of `f`, over exact integers that
    /// grow with the input, so `limits` stops each of them.
    pub(crate) fn of(f: &Polynomial<Rationals>, limits: &ComputeLimits) -> Result<Self, RunError> {
        let multiplier = denominator_lcm(f, limits)?;
        let mut terms = clear_terms(f, &multiplier, limits)?;
        let content = coefficient_gcd(&terms, limits)?;
        divide_content(&mut terms, &content, limits)?;
        Ok(ClearedPolynomial { terms })
    }

    /// The leading coefficient, or `None` for the zero polynomial.
    pub(crate) fn leading_coefficient(&self) -> Option<&BigInt> {
        self.terms.last().map(|(coeff, _)| coeff)
    }

    /// The image of the polynomial in `ring`.
    ///
    /// Every coefficient is reduced modulo the prime of `ring`, and a term
    /// the prime divides drops out. The content is 1, so the prime cannot
    /// divide every coefficient and the image is never zero.
    pub(crate) fn image(&self, ring: &PolynomialRing<PrimeField>) -> Polynomial<PrimeField> {
        let p = ring.modulus();
        let terms: Vec<Term<PrimeField>> = self
            .terms
            .iter()
            .filter_map(|(coeff, mono)| {
                let value = residue(coeff, p);
                (value != 0).then(|| Term {
                    coeff: Felt::from_residue(value),
                    mono: mono.clone(),
                })
            })
            .collect();
        Polynomial::from_sorted_terms(ring.clone(), terms)
    }

    /// The heap bytes the polynomial holds.
    pub(crate) fn heap_bytes(&self, nvars: usize) -> usize {
        let per_term = size_of::<(BigInt, Monomial)>() + heap_exps_bytes(nvars);
        self.terms.iter().fold(0usize, |bytes, (coeff, _)| {
            bytes
                .saturating_add(per_term)
                .saturating_add(limb_bytes(coeff))
        })
    }
}

fn denominator_lcm(
    polynomial: &Polynomial<Rationals>,
    limits: &ComputeLimits,
) -> Result<BigInt, RunError> {
    let mut multiplier = BigInt::one();
    for (index, term) in polynomial.terms.iter().enumerate() {
        if let Some(stop) = limits.stop_every(index) {
            return Err(stop);
        }
        multiplier = multiplier.lcm(term.coeff.denom());
    }
    Ok(multiplier)
}

fn clear_terms(
    polynomial: &Polynomial<Rationals>,
    multiplier: &BigInt,
    limits: &ComputeLimits,
) -> Result<Vec<(BigInt, Monomial)>, RunError> {
    let mut terms = Vec::with_capacity(polynomial.terms.len());
    for (index, term) in polynomial.terms.iter().enumerate() {
        if let Some(stop) = limits.stop_every(index) {
            return Err(stop);
        }
        let scaled = term.coeff.numer() * (multiplier / term.coeff.denom());
        terms.push((scaled, term.mono.clone()));
    }
    Ok(terms)
}

fn coefficient_gcd(
    terms: &[(BigInt, Monomial)],
    limits: &ComputeLimits,
) -> Result<BigInt, RunError> {
    let mut content = BigInt::zero();
    for (index, (coefficient, _)) in terms.iter().enumerate() {
        if let Some(stop) = limits.stop_every(index) {
            return Err(stop);
        }
        content = content.gcd(coefficient);
    }
    Ok(content)
}

fn divide_content(
    terms: &mut [(BigInt, Monomial)],
    content: &BigInt,
    limits: &ComputeLimits,
) -> Result<(), RunError> {
    if content.is_zero() || content.is_one() {
        return Ok(());
    }
    for (index, (coefficient, _)) in terms.iter_mut().enumerate() {
        if let Some(stop) = limits.stop_every(index) {
            return Err(stop);
        }
        *coefficient /= content;
    }
    Ok(())
}

/// What a basis over the rationals records about its own origin.
///
/// It is the basis metadata of [`Rationals`], read through
/// `GroebnerBasis::lift`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RationalMeta {
    /// The multimodular driver produced the basis.
    Lifted(ModularLift),
    /// The caller supplied the basis and the checked constructor accepted
    /// it. No modular run produced it.
    Checked,
}

/// What one multimodular run did, and what its result establishes.
///
/// The counters are exact. They count consumed runs, so
/// `primes_folded + primes_discarded` equals `primes_consumed`.
/// `primes_skipped` counts primes rejected before any run and stays
/// outside that sum.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModularLift {
    /// The runs that finished and were consumed.
    pub primes_consumed: usize,
    /// The primes that divided the leading coefficient of a cleared
    /// generator. No run started for them.
    pub primes_skipped: usize,
    /// The consumed runs in the final accumulator.
    pub primes_folded: usize,
    /// The consumed runs outside the final accumulator.
    pub primes_discarded: usize,
    /// The runs that left the lift unchanged at the stop.
    pub confirming_primes: usize,
    /// The bit length of the product of the folded primes.
    pub modulus_bits: u64,
    /// What the run establishes about the basis.
    pub established: Established,
}

/// What a rational run establishes about the basis it returns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Established {
    /// The lift did not change over the confirming primes. Nothing about
    /// the ideal follows.
    Unchanged,
    /// The lift did not change. Every generator reduces to zero modulo
    /// the basis, and so does every S-polynomial of the basis. The basis is
    /// then the reduced Gröbner basis of an ideal that contains the input
    /// ideal. Equality is not established: the basis `{1}` passes both
    /// tests.
    ContainsInput,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::PolynomialRing;

    fn ops() -> RationalOps {
        RationalOps::new()
    }

    fn rational(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    /// The unreduced fraction `numerator / denominator`.
    fn fraction(numerator: i64, denominator: i64) -> Coefficient {
        Coefficient::Fraction {
            numerator: BigInt::from(numerator),
            denominator: BigInt::from(denominator),
        }
    }

    #[test]
    fn arithmetic_is_exact() {
        let ops = ops();
        let a = rational(1, 2);
        let b = rational(1, 3);
        assert_eq!(ops.add(&a, &b), Some(rational(5, 6)));
        assert_eq!(ops.sub(&a, &b), Some(rational(1, 6)));
        assert_eq!(ops.mul(&a, &b), rational(1, 6));
        assert_eq!(ops.neg(&a), rational(-1, 2));
        assert_eq!(ops.inv(&rational(-2, 3)), rational(-3, 2));
        assert!(ops.is_one(&ops.mul(&a, &ops.inv(&a))));
    }

    #[test]
    fn a_result_of_zero_is_reported_as_none() {
        let ops = ops();
        let a = rational(3, 4);
        assert_eq!(ops.add(&a, &ops.neg(&a)), None);
        assert_eq!(ops.sub(&a, &a), None);
    }

    #[test]
    fn conversion_reduces_and_puts_the_sign_in_the_numerator() {
        let ops = ops();
        assert_eq!(
            ops.convert(&fraction(3, 6)),
            Ok(Some(BigRational::new(BigInt::from(1), BigInt::from(2))))
        );
        let converted = ops
            .convert(&fraction(1, -2))
            .expect("the denominator is not zero")
            .expect("the value is not zero");
        assert_eq!(converted, rational(-1, 2));
        assert!(converted.denom().is_positive());
        assert_eq!(ops.convert(&fraction(0, 3)), Ok(None));
        assert_eq!(
            ops.convert(&fraction(1, 0)),
            Err(RingError::ZeroDenominator)
        );
    }

    #[test]
    fn every_integer_width_is_a_coefficient() {
        let ring = PolynomialRing::rationals(["x"]).expect("the name is a variable name");
        let expected = ring
            .polynomial([(Coefficient::Integer(BigInt::from(u64::MAX)), [0])])
            .expect("the exponent vector matches the ring");
        assert_eq!(
            ring.polynomial([(u64::MAX, [0])])
                .expect("the exponent vector matches the ring"),
            expected
        );
        assert_eq!(expected.to_string(), format!("{}", u64::MAX));
    }

    #[test]
    fn a_negative_coefficient_writes_the_sign_once() {
        let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
        let f = ring
            .polynomial([
                (fraction(-1, 2), [1, 0]),
                (Coefficient::Small(-1), [0, 1]),
                (Coefficient::Small(3), [0, 0]),
            ])
            .expect("the exponent vectors match the ring");
        assert_eq!(f.to_string(), "-1/2*x - y + 3");
    }

    #[test]
    fn clearing_multiplies_by_the_denominators_and_divides_the_content() {
        let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
        let f = ring
            .polynomial([(fraction(1, 2), [2, 0]), (fraction(1, 3), [0, 1])])
            .expect("the exponent vectors match the ring");
        let cleared =
            ClearedPolynomial::of(&f, &ComputeLimits::default()).expect("no limit stops it");
        // 6*f is 3*x^2 + 2*y, whose content is 1.
        assert_eq!(cleared.leading_coefficient(), Some(&BigInt::from(3)));
        let image = cleared
            .image(&PolynomialRing::prime_field(7, ["x", "y"]).expect("the modulus is prime"));
        assert_eq!(image.to_string(), "3*x^2 + 2*y");
    }

    #[test]
    fn a_content_of_more_than_one_divides_out() {
        let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
        let f = ring
            .polynomial([
                (Coefficient::Small(6), [2, 0]),
                (Coefficient::Small(-4), [0, 1]),
            ])
            .expect("the exponent vectors match the ring");
        let cleared =
            ClearedPolynomial::of(&f, &ComputeLimits::default()).expect("no limit stops it");
        assert_eq!(cleared.leading_coefficient(), Some(&BigInt::from(3)));
        let image = cleared
            .image(&PolynomialRing::prime_field(11, ["x", "y"]).expect("the modulus is prime"));
        assert_eq!(image.to_string(), "3*x^2 + 9*y");
    }

    #[test]
    fn a_coefficient_the_prime_divides_leaves_the_image() {
        let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
        let f = ring
            .polynomial([
                (Coefficient::Small(1), [2, 0]),
                (Coefficient::Small(7), [0, 1]),
            ])
            .expect("the exponent vectors match the ring");
        let image = ClearedPolynomial::of(&f, &ComputeLimits::default())
            .expect("no limit stops it")
            .image(&PolynomialRing::prime_field(7, ["x", "y"]).expect("the modulus is prime"));
        assert_eq!(image.to_string(), "x^2");
    }

    /// The clearing loops read the deadline, so an exhausted one stops
    /// them.
    #[test]
    fn clearing_stops_on_an_exhausted_deadline() {
        let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
        let f = ring
            .polynomial([(fraction(1, 2), [2, 0]), (fraction(1, 3), [0, 1])])
            .expect("the exponent vectors match the ring");
        let limits =
            ComputeLimits::of_budget(&crate::Budget::new().timeout(std::time::Duration::ZERO));
        assert_eq!(
            ClearedPolynomial::of(&f, &limits).map(|_| ()),
            Err(RunError::Compute(crate::ComputeError::Timeout))
        );
    }

    #[test]
    fn the_zero_polynomial_clears_to_no_term() {
        let ring = PolynomialRing::rationals(["x"]).expect("the name is a variable name");
        let cleared = ClearedPolynomial::of(&ring.zero(), &ComputeLimits::default())
            .expect("no limit stops it");
        assert_eq!(cleared.leading_coefficient(), None);
        assert_eq!(cleared.heap_bytes(1), 0);
    }

    #[test]
    fn the_heap_estimate_grows_with_the_value() {
        let ops = ops();
        // One limb each for 1 and 2.
        assert_eq!(ops.heap_bytes(&rational(1, 2)), 16);
        // 2^500 holds 501 bits, which is eight limbs, and 3 holds one.
        let large = BigRational::new(BigInt::from(1) << 500, BigInt::from(3));
        assert_eq!(ops.heap_bytes(&large), 72);
        // A numerator of zero holds no digits.
        let zero = BigRational::new_raw(BigInt::from(0), BigInt::from(1));
        assert_eq!(ops.heap_bytes(&zero), 8);
    }
}
