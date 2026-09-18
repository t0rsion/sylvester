use std::time::Duration;

use sylvester::{
    ArithmeticError, Budget, CancellationToken, Coefficient, PolynomialRing, RingError,
};

#[test]
fn prime_results_match_independent_term_construction() {
    let ring = PolynomialRing::prime_field(11, ["x", "y"]).expect("11 is prime");
    let left = ring
        .polynomial([(1, [1, 0]), (1, [0, 1]), (1, [0, 0])])
        .expect("left terms fit the ring");
    let right = ring
        .polynomial([(2, [1, 0]), (3, [0, 1]), (4, [0, 0])])
        .expect("right terms fit the ring");
    let expected_sum = ring
        .polynomial([(3, [1, 0]), (4, [0, 1]), (5, [0, 0])])
        .expect("sum terms fit the ring");
    let expected_difference = ring
        .polynomial([(10, [1, 0]), (9, [0, 1]), (8, [0, 0])])
        .expect("difference terms fit the ring");
    let expected_product = ring
        .polynomial([
            (2, [2, 0]),
            (5, [1, 1]),
            (6, [1, 0]),
            (3, [0, 2]),
            (7, [0, 1]),
            (4, [0, 0]),
        ])
        .expect("product terms fit the ring");

    assert_eq!(
        left.try_add(&right, Budget::new())
            .expect("unlimited addition"),
        expected_sum
    );
    assert_eq!(
        left.try_sub(&right, Budget::new())
            .expect("unlimited subtraction"),
        expected_difference
    );
    assert_eq!(
        left.try_mul(&right, Budget::new())
            .expect("unlimited multiplication"),
        expected_product
    );
}

#[test]
fn prime_field_arithmetic_obeys_ring_identities() {
    let ring = PolynomialRing::prime_field(101, ["x", "y"]).expect("101 is prime");
    let f = ring.parse_polynomial("x^2 + 3*x*y - 2").expect("parses");
    let g = ring.parse_polynomial("x - y + 4").expect("parses");
    let sum = f.try_add(&g, Budget::new()).expect("unlimited addition");

    assert_eq!(
        sum.try_sub(&g, Budget::new())
            .expect("unlimited subtraction"),
        f
    );
    assert_eq!(
        f.try_mul(&g, Budget::new())
            .expect("unlimited multiplication"),
        g.try_mul(&f, Budget::new())
            .expect("unlimited multiplication")
    );
    assert_eq!(
        f.try_pow(2, Budget::new()).expect("unlimited power"),
        f.try_mul(&f, Budget::new())
            .expect("unlimited multiplication")
    );
    assert_eq!(
        f.try_neg(Budget::new())
            .expect("unlimited negation")
            .try_neg(Budget::new())
            .expect("unlimited negation"),
        f
    );
    assert_eq!(f.try_scale(1, Budget::new()).expect("unlimited scaling"), f);
    assert_eq!(
        f.try_pow(0, Budget::new()).expect("unlimited power"),
        ring.one()
    );
}

#[test]
fn rational_arithmetic_keeps_exact_coefficients() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("valid variable names");
    let f = ring.parse_polynomial("1/2*x + 1/3*y").expect("parses");
    let g = ring.parse_polynomial("2/3*x - 1/5").expect("parses");

    assert_eq!(
        f.try_add(&g, Budget::new())
            .expect("unlimited addition")
            .to_string(),
        "7/6*x + 1/3*y - 1/5"
    );
    assert_eq!(
        f.try_scale(
            Coefficient::Fraction {
                numerator: 3.into(),
                denominator: 2.into(),
            },
            Budget::new(),
        )
        .expect("unlimited scaling")
        .to_string(),
        "3/4*x + 1/2*y"
    );
    assert_eq!(
        f.try_pow(2, Budget::new())
            .expect("unlimited power")
            .to_string(),
        "1/4*x^2 + 1/3*x*y + 1/9*y^2"
    );
}

#[test]
fn rational_memory_budget_charges_coefficient_workspace() {
    let ring = PolynomialRing::rationals(["x"]).expect("valid variable names");
    let value = ring.parse_polynomial("1/2*x").expect("parses");

    assert_eq!(
        value.try_mul(&value, Budget::new().memory_limit(480)),
        Err(ArithmeticError::MemoryLimitExceeded)
    );
}

#[test]
fn arithmetic_reports_ring_and_coefficient_errors() {
    let left_ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    let right_ring = PolynomialRing::prime_field(11, ["x"]).expect("11 is prime");
    let left = left_ring.parse_polynomial("x").expect("parses");
    let right = right_ring.parse_polynomial("x").expect("parses");

    assert_eq!(
        left.try_add(&right, Budget::new()),
        Err(ArithmeticError::RingMismatch)
    );
    assert_eq!(
        left.try_scale(
            Coefficient::Fraction {
                numerator: 1.into(),
                denominator: 7.into(),
            },
            Budget::new(),
        ),
        Err(ArithmeticError::CoefficientConversion(
            RingError::CoefficientNotInvertible {
                denominator: 7.into(),
            }
        ))
    );
}

#[test]
fn arithmetic_reports_exponent_overflow_before_product_allocation() {
    let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    let high = ring.polynomial([(1, [u16::MAX])]).expect("one exponent");
    let x = ring.parse_polynomial("x").expect("parses");

    assert_eq!(
        high.try_mul(&x, Budget::new()),
        Err(ArithmeticError::ExponentLimit {
            limit: u16::MAX as u32
        })
    );
    assert_eq!(
        high.try_pow(2, Budget::new()),
        Err(ArithmeticError::ExponentLimit {
            limit: u16::MAX as u32
        })
    );
}

#[test]
fn exhausted_budgets_apply_to_zero_work_and_results() {
    let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    let zero = ring.zero();
    let x = ring.parse_polynomial("x").expect("parses");
    let expired = Budget::new().timeout(Duration::ZERO);
    let cancelled = CancellationToken::new();
    cancelled.cancel();

    assert_eq!(zero.try_add(&zero, expired), Err(ArithmeticError::Timeout));
    assert_eq!(
        zero.try_add(&zero, Budget::new().cancellation(cancelled)),
        Err(ArithmeticError::Timeout)
    );
    assert_eq!(
        x.try_add(&x, Budget::new().memory_limit(0)),
        Err(ArithmeticError::MemoryLimitExceeded)
    );
    assert_eq!(
        zero.try_pow(0, Budget::new().memory_limit(0)),
        Err(ArithmeticError::MemoryLimitExceeded)
    );
}

#[test]
fn arithmetic_memory_budget_charges_retained_term_capacity() {
    let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    let cancelled = ring
        .polynomial([(1, [0]), (-1, [0])])
        .expect("terms fit the ring");

    assert!(cancelled.is_zero());
    assert_eq!(
        cancelled.try_add(&cancelled, Budget::new().memory_limit(0)),
        Err(ArithmeticError::MemoryLimitExceeded)
    );
}
