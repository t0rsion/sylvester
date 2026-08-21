//! Property tests for the ring constructors and the text syntax.
//!
//! `Display` writes what `parse_polynomial` reads, so every polynomial must
//! survive the trip through its own text.

use proptest::prelude::*;
use sylvester::{Polynomial, PolynomialRing};

const PRIMES: [u64; 4] = [2, 3, 7, 32003];

/// A ring with one to four variables over one of the test primes.
fn any_ring() -> impl Strategy<Value = PolynomialRing> {
    (0..PRIMES.len(), 1usize..=4).prop_map(|(index, nvars)| {
        let names: Vec<String> = ["x", "y", "z", "w"][..nvars]
            .iter()
            .map(|name| name.to_string())
            .collect();
        PolynomialRing::prime_field(PRIMES[index], names).expect("the modulus is prime")
    })
}

/// A polynomial of a random ring, with up to six terms of degree at most
/// nine in each variable.
fn any_polynomial() -> impl Strategy<Value = Polynomial> {
    any_ring().prop_flat_map(|ring| {
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

proptest! {
    #[test]
    fn text_round_trips_through_display_and_parse(poly in any_polynomial()) {
        let text = poly.to_string();
        let read = poly
            .ring()
            .parse_polynomial(&text)
            .unwrap_or_else(|error| panic!("\"{text}\" must parse: {error}"));
        prop_assert_eq!(read, poly);
    }

    #[test]
    fn display_is_a_function_of_the_polynomial(poly in any_polynomial()) {
        let once = poly.to_string();
        let twice = poly
            .ring()
            .parse_polynomial(&once)
            .expect("the text parses")
            .to_string();
        prop_assert_eq!(once, twice);
    }

    #[test]
    fn a_built_polynomial_is_reduced_and_ordered(poly in any_polynomial()) {
        let modulus = poly.ring().modulus();
        let mut monomials: Vec<Vec<u16>> = Vec::new();
        for (coeff, exps) in poly.terms() {
            prop_assert!(coeff >= 1 && coeff < modulus, "coefficient {coeff} is outside [1, p)");
            prop_assert_eq!(exps.len(), poly.ring().nvars());
            monomials.push(exps.to_vec());
        }

        for pair in monomials.windows(2) {
            prop_assert!(
                grevlex_cmp(&pair[0], &pair[1]) == std::cmp::Ordering::Greater,
                "terms must run strictly descending"
            );
        }
        prop_assert_eq!(poly.is_zero(), poly.terms().len() == 0);
    }

    #[test]
    fn whitespace_between_tokens_never_changes_the_value(poly in any_polynomial()) {
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
        prop_assert_eq!(
            poly.ring().parse_polynomial(&spaced).expect("the text parses"),
            poly
        );
    }
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
