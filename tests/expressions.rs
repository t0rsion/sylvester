use std::time::Duration;

use sylvester::{Budget, ExpressionError, ParseError, PolynomialRing};

#[test]
fn generators_follow_the_variable_order() {
    let ring = PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime");
    let generators = ring.generators();
    assert_eq!(generators.len(), 2);
    assert_eq!(generators[0].to_string(), "x");
    assert_eq!(generators[1].to_string(), "y");
    assert_eq!(ring.generator(0), Some(generators[0].clone()));
    assert!(ring.generator(2).is_none());
}

#[test]
fn flat_input_has_the_same_value_in_both_parsers() {
    let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"]).expect("32003 is prime");
    let text = "x^2*y - 3*z + 1";
    assert_eq!(
        ring.parse_polynomial_with_budget(text, Budget::new())
            .expect("the expression parses"),
        ring.parse_polynomial(text).expect("the flat syntax parses")
    );
}

#[test]
fn flat_fraction_has_the_same_field_value_in_both_parsers() {
    let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    let text = "3/6*x - 1/2";
    assert_eq!(
        ring.parse_polynomial_with_budget(text, Budget::new())
            .expect("the expression parses"),
        ring.parse_polynomial(text).expect("the flat syntax parses")
    );
}

#[test]
fn parentheses_and_precedence_expand_products() {
    let ring = PolynomialRing::prime_field(32003, ["x", "y"]).expect("32003 is prime");
    let parsed = ring
        .parse_polynomial_with_budget("(x + y) * (x - y)", Budget::new())
        .expect("the expression parses");
    let expected = ring
        .polynomial([(1, [2, 0]), (-1, [0, 2])])
        .expect("the terms fit");
    assert_eq!(parsed, expected);
    assert_eq!(
        ring.parse_polynomial_with_budget("x + y * x", Budget::new())
            .expect("the expression parses")
            .to_string(),
        "x*y + x"
    );
}

#[test]
fn unary_signs_and_fractions_stay_exact() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are valid");
    let parsed = ring
        .parse_polynomial_with_budget("-(1 / 2) * (x + -y)", Budget::new())
        .expect("the expression parses");
    assert_eq!(parsed.to_string(), "-1/2*x + 1/2*y");
    assert_eq!(
        ring.parse_polynomial_with_budget("(1/2 + 1/2) * x", Budget::new())
            .expect("the expression parses")
            .to_string(),
        "x"
    );
}

#[test]
fn double_star_is_a_power_operator() {
    let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    let parsed = ring
        .parse_polynomial_with_budget("(x + 1) ** 2", Budget::new())
        .expect("the expression parses");
    assert_eq!(parsed.to_string(), "x^2 + 2*x + 1");
}

#[test]
fn division_between_polynomials_is_rejected() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are valid");
    assert_eq!(
        ring.parse_polynomial_with_budget("x / y", Budget::new()),
        Err(ExpressionError::Parse(ParseError::UnexpectedCharacter {
            character: '/',
            position: 2,
        }))
    );
}

#[test]
fn expansion_budget_is_call_wide() {
    let ring = PolynomialRing::prime_field(32003, ["x", "y"]).expect("32003 is prime");
    assert_eq!(
        ring.parse_polynomial_with_budget("(x + y)^8", Budget::new().memory_limit(1),),
        Err(ExpressionError::MemoryLimitExceeded)
    );
    assert!(
        ring.parse_polynomial_with_budget("(x + y)^8", Budget::new())
            .is_ok()
    );
}

#[test]
fn timeout_is_reported_before_expansion() {
    let ring = PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime");
    assert_eq!(
        ring.parse_polynomial_with_budget("(x + y)^8", Budget::new().timeout(Duration::ZERO),),
        Err(ExpressionError::Timeout)
    );
}

#[test]
fn exponent_width_is_a_typed_expression_error() {
    let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    assert_eq!(
        ring.parse_polynomial_with_budget("x^65536", Budget::new()),
        Err(ExpressionError::ExponentLimit { limit: 65535 })
    );
}

#[test]
fn an_exponent_past_u64_is_rejected_without_saturation() {
    let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    assert_eq!(
        ring.parse_polynomial_with_budget("(-1)^18446744073709551616", Budget::new(),),
        Err(ExpressionError::ExponentLimit { limit: 65535 })
    );
}

#[test]
fn zero_to_zero_checks_the_one_result_against_memory() {
    let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    assert_eq!(
        ring.parse_polynomial_with_budget("0^0", Budget::new().memory_limit(0)),
        Err(ExpressionError::MemoryLimitExceeded)
    );
}

#[test]
fn nesting_limit_prevents_stack_overflow() {
    let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
    let text = format!("{}x{}", "(".repeat(300), ")".repeat(300));
    assert!(matches!(
        ring.parse_polynomial_with_budget(&text, Budget::new()),
        Err(ExpressionError::NestingLimit { .. })
    ));
}
