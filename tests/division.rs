//! Quotients and remainders from multivariate division.

use std::time::Duration;

use num_rational::BigRational;
use sylvester::{
    BasisError, Budget, GroebnerBasis, NormalFormError, Polynomial, PolynomialRing, Rationals,
};

const MODULUS: u64 = 32003;

fn prime_ring() -> PolynomialRing {
    PolynomialRing::prime_field(MODULUS, ["x", "y"])
        .expect("the modulus is prime and the names are distinct")
}

fn rational_ring() -> PolynomialRing<Rationals> {
    PolynomialRing::rationals(["x", "y"]).expect("the names are distinct")
}

fn prime_basis(ring: &PolynomialRing) -> GroebnerBasis {
    GroebnerBasis::<sylvester::PrimeField>::from_polynomials(
        ring,
        vec![
            ring.parse_polynomial("y^2").expect("the syntax holds"),
            ring.parse_polynomial("x").expect("the syntax holds"),
        ],
        Budget::new(),
    )
    .expect("the list is a reduced basis")
}

fn rational_basis(ring: &PolynomialRing<Rationals>) -> GroebnerBasis<Rationals> {
    GroebnerBasis::<Rationals>::from_polynomials(
        ring,
        vec![
            ring.parse_polynomial("y^2").expect("the syntax holds"),
            ring.parse_polynomial("x").expect("the syntax holds"),
        ],
        Budget::new(),
    )
    .expect("the list is a reduced basis")
}

fn prime_product(ring: &PolynomialRing, left: &Polynomial, right: &Polynomial) -> Polynomial {
    let mut terms = Vec::new();
    for (left_term, left_exps) in left.terms() {
        for (right_term, right_exps) in right.terms() {
            let exps = left_exps
                .iter()
                .zip(right_exps)
                .map(|(a, b)| a + b)
                .collect::<Vec<_>>();
            let coeff = (left_term.value() * right_term.value()) % MODULUS;
            terms.push((coeff, exps));
        }
    }
    ring.polynomial(terms).expect("the exponent vectors fit")
}

fn rational_product(
    ring: &PolynomialRing<Rationals>,
    left: &Polynomial<Rationals>,
    right: &Polynomial<Rationals>,
) -> Polynomial<Rationals> {
    let mut terms = Vec::new();
    for (left_term, left_exps) in left.terms() {
        for (right_term, right_exps) in right.terms() {
            let exps = left_exps
                .iter()
                .zip(right_exps)
                .map(|(a, b)| a + b)
                .collect::<Vec<_>>();
            let coeff: BigRational = left_term * right_term;
            terms.push((coeff, exps));
        }
    }
    ring.polynomial(terms).expect("the exponent vectors fit")
}

fn prime_sum(ring: &PolynomialRing, polynomials: &[Polynomial]) -> Polynomial {
    let terms = polynomials
        .iter()
        .flat_map(|poly| {
            poly.terms()
                .map(|(coeff, exps)| (coeff.value(), exps.to_vec()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    ring.polynomial(terms).expect("the exponent vectors fit")
}

fn rational_sum(
    ring: &PolynomialRing<Rationals>,
    polynomials: &[Polynomial<Rationals>],
) -> Polynomial<Rationals> {
    let terms = polynomials
        .iter()
        .flat_map(|poly| {
            poly.terms()
                .map(|(coeff, exps)| (coeff.clone(), exps.to_vec()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    ring.polynomial(terms).expect("the exponent vectors fit")
}

#[test]
fn prime_division_rebuilds_the_input_and_aligns_quotients() {
    let ring = prime_ring();
    let basis = prime_basis(&ring);
    let input = ring
        .parse_polynomial("x*y + 2*x + 3*y + 7")
        .expect("the syntax holds");
    let result = basis
        .divide(&input, Budget::new())
        .expect("the division fits the unlimited budget");

    assert_eq!(result.quotients().len(), basis.len());
    assert!(result.quotients()[0].is_zero());
    assert_eq!(
        result.quotients()[1],
        ring.parse_polynomial("y + 2").unwrap()
    );
    assert_eq!(
        result.remainder(),
        &ring.parse_polynomial("3*y + 7").unwrap()
    );
    assert_eq!(
        prime_sum(
            &ring,
            &[
                prime_product(&ring, &result.quotients()[0], &basis[0]),
                prime_product(&ring, &result.quotients()[1], &basis[1]),
                result.remainder().clone(),
            ],
        ),
        input
    );
}

#[test]
fn rational_division_rebuilds_the_input_and_aligns_quotients() {
    let ring = rational_ring();
    let basis = rational_basis(&ring);
    let input = ring
        .parse_polynomial("1/2*x*y + 2*x + 3*y + 7")
        .expect("the syntax holds");
    let result = basis
        .divide(&input, Budget::new())
        .expect("the division fits the unlimited budget");

    assert_eq!(result.quotients().len(), basis.len());
    assert!(result.quotients()[0].is_zero());
    assert_eq!(
        result.quotients()[1],
        ring.parse_polynomial("1/2*y + 2").unwrap()
    );
    assert_eq!(
        result.remainder(),
        &ring.parse_polynomial("3*y + 7").unwrap()
    );
    assert_eq!(
        rational_sum(
            &ring,
            &[
                rational_product(&ring, &result.quotients()[0], &basis[0]),
                rational_product(&ring, &result.quotients()[1], &basis[1]),
                result.remainder().clone(),
            ],
        ),
        input
    );
}

#[test]
fn prime_zero_unit_and_empty_bases_have_the_expected_parts() {
    let prime = prime_ring();
    let input = prime.parse_polynomial("x^2 + y + 1").unwrap();
    let zero = prime.zero();
    let empty =
        GroebnerBasis::<sylvester::PrimeField>::from_polynomials(&prime, Vec::new(), Budget::new())
            .unwrap();
    let empty_result = empty.divide(&input, Budget::new()).unwrap();
    assert!(empty_result.quotients().is_empty());
    assert_eq!(empty_result.remainder(), &input);
    let zero_result = empty.divide(&zero, Budget::new()).unwrap();
    assert!(zero_result.quotients().is_empty());
    assert!(zero_result.remainder().is_zero());

    let unit = GroebnerBasis::<sylvester::PrimeField>::from_polynomials(
        &prime,
        vec![prime.one()],
        Budget::new(),
    )
    .unwrap();
    let unit_result = unit.divide(&input, Budget::new()).unwrap();
    assert_eq!(unit_result.quotients(), std::slice::from_ref(&input));
    assert!(unit_result.remainder().is_zero());
    let unit_zero = unit.divide(&zero, Budget::new()).unwrap();
    assert!(unit_zero.quotients()[0].is_zero());
    assert!(unit_zero.remainder().is_zero());
}

#[test]
fn rational_zero_unit_and_empty_bases_have_the_expected_parts() {
    let rational = rational_ring();
    let rational_input = rational.parse_polynomial("x + 1").unwrap();
    let rational_zero = rational.zero();
    let rational_empty =
        GroebnerBasis::<Rationals>::from_polynomials(&rational, Vec::new(), Budget::new()).unwrap();
    let rational_result = rational_empty
        .divide(&rational_input, Budget::new())
        .unwrap();
    assert!(rational_result.quotients().is_empty());
    assert_eq!(rational_result.remainder(), &rational_input);
    let rational_zero_result = rational_empty
        .divide(&rational_zero, Budget::new())
        .unwrap();
    assert!(rational_zero_result.quotients().is_empty());
    assert!(rational_zero_result.remainder().is_zero());

    let rational_unit = GroebnerBasis::<Rationals>::from_polynomials(
        &rational,
        vec![rational.one()],
        Budget::new(),
    )
    .unwrap();
    let rational_unit_result = rational_unit
        .divide(&rational_input, Budget::new())
        .unwrap();
    assert_eq!(rational_unit_result.quotients(), &[rational_input]);
    assert!(rational_unit_result.remainder().is_zero());
}

#[test]
fn division_reports_overflow_timeout_and_memory_limits() {
    let ring = prime_ring();
    let divisor = ring
        .polynomial([(1u64, [65535, 0]), (1u64, [0, 65535])])
        .unwrap();
    let overflow_basis = GroebnerBasis::<sylvester::PrimeField>::from_polynomials(
        &ring,
        vec![divisor],
        Budget::new(),
    )
    .expect("the divisor is a reduced basis");
    let overflow_input = ring.polynomial([(1u64, [65535, 1])]).unwrap();
    assert_eq!(
        overflow_basis.divide(&overflow_input, Budget::new()),
        Err(NormalFormError::ExponentLimit { limit: 65535 })
    );

    let basis = prime_basis(&ring);
    assert_eq!(
        basis.divide(
            &ring.parse_polynomial("x + y").unwrap(),
            Budget::new().timeout(Duration::ZERO),
        ),
        Err(NormalFormError::Timeout)
    );
    assert_eq!(
        basis.divide(
            &ring.parse_polynomial("x + y").unwrap(),
            Budget::new().memory_limit(0),
        ),
        Err(NormalFormError::MemoryLimitExceeded)
    );
}

#[test]
fn empty_inputs_still_poll_an_expired_budget() {
    let ring = prime_ring();
    let basis =
        GroebnerBasis::<sylvester::PrimeField>::from_polynomials(&ring, Vec::new(), Budget::new())
            .unwrap();
    let zero = ring.zero();
    let expired = Budget::new().timeout(Duration::ZERO);

    assert_eq!(
        basis.divide(&zero, expired.clone()),
        Err(NormalFormError::Timeout)
    );
    assert_eq!(
        basis.normal_form(&zero, expired.clone()),
        Err(NormalFormError::Timeout)
    );
    assert_eq!(
        GroebnerBasis::<sylvester::PrimeField>::from_polynomials(&ring, Vec::new(), expired),
        Err(BasisError::Timeout)
    );
}

#[test]
fn division_rejects_a_polynomial_from_another_ring() {
    let ring = prime_ring();
    let basis = prime_basis(&ring);
    let other = PolynomialRing::prime_field(7, ["x", "y"]).unwrap();
    let input = other.parse_polynomial("x + 1").unwrap();

    assert_eq!(
        basis.divide(&input, Budget::new()),
        Err(NormalFormError::RingMismatch)
    );
}
