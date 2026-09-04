//! Property tests for the ring constructors and the text syntax.
//!
//! `Display` writes what `parse_polynomial` reads, so every polynomial must
//! survive the trip through its own text. Both domains are covered, because
//! one reader and one writer serve both.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use proptest::prelude::*;
use sylvester::{Coefficient, Domain, Polynomial, PolynomialRing, Rationals};

const PRIMES: [u64; 4] = [2, 3, 7, 32003];

/// The variable names of a ring with `nvars` variables.
fn names(nvars: usize) -> Vec<String> {
    ["x", "y", "z", "w"][..nvars]
        .iter()
        .map(|name| name.to_string())
        .collect()
}

/// A ring with one to four variables over one of the test primes.
fn any_prime_ring() -> impl Strategy<Value = PolynomialRing> {
    (0..PRIMES.len(), 1usize..=4).prop_map(|(index, nvars)| {
        PolynomialRing::prime_field(PRIMES[index], names(nvars)).expect("the modulus is prime")
    })
}

/// A ring with one to four variables over the rationals.
fn any_rational_ring() -> impl Strategy<Value = PolynomialRing<Rationals>> {
    (1usize..=4).prop_map(|nvars| {
        PolynomialRing::rationals(names(nvars)).expect("the names are variable names")
    })
}

/// A polynomial of a random prime-field ring, with up to six terms of
/// degree at most nine in each variable.
fn any_prime_polynomial() -> impl Strategy<Value = Polynomial> {
    any_prime_ring().prop_flat_map(|ring| {
        let nvars = ring.nvars();
        let terms = proptest::collection::vec(
            (
                -1000i64..1000,
                proptest::collection::vec(0u16..10, nvars..=nvars),
            ),
            0..6,
        );
        terms.prop_map(move |terms| {
            ring.polynomial(terms)
                .expect("the exponent vectors match the ring")
        })
    })
}

/// A rational coefficient: a small integer, an integer past the width of an
/// `i64`, or a fraction with a small or a large denominator.
fn any_coefficient() -> impl Strategy<Value = Coefficient> {
    prop_oneof![
        (-1000i64..1000).prop_map(Coefficient::Small),
        any::<i64>().prop_map(|value| Coefficient::Integer(BigInt::from(value) << 70)),
        (any::<i64>(), 1i64..1000).prop_map(|(numerator, denominator)| Coefficient::Fraction {
            numerator: BigInt::from(numerator),
            denominator: BigInt::from(denominator),
        }),
        (any::<i64>(), 1i64..1000).prop_map(|(numerator, denominator)| Coefficient::Fraction {
            numerator: BigInt::from(numerator) << 70,
            denominator: (BigInt::from(denominator) << 70) + 1,
        }),
    ]
}

/// A polynomial of a random rational ring, with up to six terms of degree
/// at most nine in each variable.
fn any_rational_polynomial() -> impl Strategy<Value = Polynomial<Rationals>> {
    any_rational_ring().prop_flat_map(|ring| {
        let nvars = ring.nvars();
        let terms = proptest::collection::vec(
            (
                any_coefficient(),
                proptest::collection::vec(0u16..10, nvars..=nvars),
            ),
            0..6,
        );
        terms.prop_map(move |terms| {
            ring.polynomial(terms)
                .expect("the exponent vectors match the ring")
        })
    })
}

/// The text parses back to the same polynomial, and writes the same text
/// again.
fn round_trips<D: Domain>(poly: &Polynomial<D>) -> Result<(), TestCaseError> {
    let text = poly.to_string();
    let read = poly
        .ring()
        .parse_polynomial(&text)
        .unwrap_or_else(|error| panic!("\"{text}\" must parse: {error}"));
    prop_assert_eq!(&read, poly);
    prop_assert_eq!(read.to_string(), text);
    Ok(())
}

/// Whitespace around every operator leaves the value alone.
///
/// The fraction bar takes no whitespace, so it is not one of the operators
/// this pads.
fn survives_padding<D: Domain>(poly: &Polynomial<D>) -> Result<(), TestCaseError> {
    let text = poly.to_string();
    let mut spaced = String::new();
    for c in text.chars() {
        if matches!(c, '+' | '-' | '*' | '^') {
            spaced.push_str("  ");
            spaced.push(c);
            spaced.push_str("  ");
        } else {
            spaced.push(c);
        }
    }
    let read = poly
        .ring()
        .parse_polynomial(&spaced)
        .unwrap_or_else(|error| panic!("\"{spaced}\" must parse: {error}"));
    prop_assert_eq!(&read, poly);
    Ok(())
}

/// The monomials of the polynomial, largest first.
fn monomials<D: Domain>(poly: &Polynomial<D>) -> Vec<Vec<u16>> {
    poly.terms().map(|(_, exps)| exps.to_vec()).collect()
}

/// The terms run strictly descending under grevlex, so no monomial repeats.
fn descends(monomials: &[Vec<u16>]) -> Result<(), TestCaseError> {
    for pair in monomials.windows(2) {
        prop_assert!(
            grevlex_cmp(&pair[0], &pair[1]) == std::cmp::Ordering::Greater,
            "terms must run strictly descending"
        );
    }
    Ok(())
}

proptest! {
    #[test]
    fn text_round_trips_over_a_prime_field(poly in any_prime_polynomial()) {
        round_trips(&poly)?;
    }

    #[test]
    fn text_round_trips_over_the_rationals(poly in any_rational_polynomial()) {
        round_trips(&poly)?;
    }

    #[test]
    fn whitespace_never_changes_a_prime_field_value(poly in any_prime_polynomial()) {
        survives_padding(&poly)?;
    }

    #[test]
    fn whitespace_never_changes_a_rational_value(poly in any_rational_polynomial()) {
        survives_padding(&poly)?;
    }

    #[test]
    fn a_built_polynomial_is_reduced_and_ordered(poly in any_prime_polynomial()) {
        let modulus = poly.ring().modulus();
        for (coeff, exps) in poly.terms() {
            let coeff = coeff.value();
            prop_assert!(coeff >= 1 && coeff < modulus, "coefficient {coeff} is outside [1, p)");
            prop_assert_eq!(exps.len(), poly.ring().nvars());
        }
        descends(&monomials(&poly))?;
        prop_assert_eq!(poly.is_zero(), poly.terms().len() == 0);
    }

    #[test]
    fn a_built_rational_polynomial_is_in_lowest_terms(poly in any_rational_polynomial()) {
        for (coeff, exps) in poly.terms() {
            prop_assert!(!coeff.numer().is_zero(), "a stored coefficient is never zero");
            prop_assert!(coeff.denom().is_positive(), "the denominator carries no sign");
            prop_assert!(
                coeff.numer().gcd(coeff.denom()).is_one(),
                "{coeff} is not in lowest terms"
            );
            prop_assert_eq!(exps.len(), poly.ring().nvars());
        }
        descends(&monomials(&poly))?;
        prop_assert_eq!(poly.is_zero(), poly.terms().len() == 0);
    }
}

#[test]
fn a_coefficient_past_the_width_of_an_i64_round_trips() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
    let numerator = BigInt::from(u64::MAX) * BigInt::from(u64::MAX);
    let text = format!("{numerator}/7*x^2 - {numerator}*y");
    let poly = ring.parse_polynomial(&text).expect("the text parses");
    assert_eq!(poly.to_string(), text);
    let (coeff, _) = poly.leading_term().expect("the polynomial is not zero");
    assert_eq!(coeff.numer(), &numerator);
    assert_eq!(coeff.denom(), &BigInt::from(7));
}

/// Grevlex from the definition, written here so the test does not read the
/// engine's comparison.
fn grevlex_cmp(a: &[u16], b: &[u16]) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let da: u64 = a.iter().map(|&e| e as u64).sum();
    let db: u64 = b.iter().map(|&e| e as u64).sum();
    match da.cmp(&db) {
        Ordering::Equal => {}
        unequal => return unequal,
    }
    for (ea, eb) in a.iter().zip(b).rev() {
        match ea.cmp(eb) {
            Ordering::Equal => continue,
            ord => return ord.reverse(),
        }
    }
    Ordering::Equal
}
