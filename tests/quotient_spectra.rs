//! Independent checks for finite quotient matrices and spectra.

use sylvester::{
    Budget, ComputeOptions, GroebnerBasis, PolynomialRing, PrimeField, RationalOptions,
    UnivariatePolynomial,
};

#[test]
fn characteristic_polynomial_handles_degree_above_the_characteristic() {
    let ring = PolynomialRing::prime_field(2, ["x"]).expect("2 is prime");
    let x = ring.parse_polynomial("x").expect("x parses");
    let ideal = ring
        .ideal([ring
            .parse_polynomial("x^3 + x")
            .expect("the generator parses")])
        .expect("the generator belongs to the ring");
    let basis = ideal
        .groebner_basis(ComputeOptions::new())
        .expect("the basis computes");
    let quotient = basis
        .finite_quotient(sylvester::Budget::new())
        .expect("the quotient is finite");

    let matrix = quotient
        .multiplication_matrix(&x, sylvester::Budget::new())
        .expect("the quotient matrix computes");
    assert_eq!(matrix.dimension(), 3);
    assert_eq!(matrix.entries().len(), 9);

    let characteristic = quotient
        .characteristic_polynomial(&x, sylvester::Budget::new())
        .expect("the characteristic polynomial computes");
    let minimal = quotient
        .minimal_polynomial(&x, sylvester::Budget::new())
        .expect("the minimal polynomial computes");
    assert_eq!(characteristic.to_string(), "t^3 + t");
    assert_eq!(minimal.to_string(), "t^3 + t");
    assert_eq!(characteristic.degree(), Some(3));
}

#[test]
fn characteristic_polynomial_handles_degree_above_three() {
    let ring = PolynomialRing::prime_field(3, ["x"]).expect("3 is prime");
    let x = ring.parse_polynomial("x").expect("x parses");
    let ideal = ring
        .ideal([ring
            .parse_polynomial("x^4 - 1")
            .expect("the generator parses")])
        .expect("the generator belongs to the ring");
    let basis = ideal
        .groebner_basis(ComputeOptions::new())
        .expect("the basis computes");
    let quotient = basis
        .finite_quotient(sylvester::Budget::new())
        .expect("the quotient is finite");

    let characteristic = quotient
        .characteristic_polynomial(&x, sylvester::Budget::new())
        .expect("the characteristic polynomial computes");
    let minimal = quotient
        .minimal_polynomial(&x, sylvester::Budget::new())
        .expect("the minimal polynomial computes");
    assert_eq!(characteristic.to_string(), "t^4 + 2");
    assert_eq!(minimal.to_string(), "t^4 + 2");
}

#[test]
fn rational_spectra_are_exact() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the ring builds");
    let x = ring.parse_polynomial("x").expect("x parses");
    let ideal = ring
        .ideal([
            ring.parse_polynomial("x^2 + y^2 - 1")
                .expect("the first generator parses"),
            ring.parse_polynomial("4*x*y - 1")
                .expect("the second generator parses"),
        ])
        .expect("the generators belong to the ring");
    let basis = ideal
        .groebner_basis(RationalOptions::new())
        .expect("the rational basis computes");
    let quotient = basis
        .finite_quotient(sylvester::Budget::new())
        .expect("the quotient is finite");

    let characteristic = quotient
        .characteristic_polynomial(&x, sylvester::Budget::new())
        .expect("the characteristic polynomial computes");
    let minimal = quotient
        .minimal_polynomial(&x, sylvester::Budget::new())
        .expect("the minimal polynomial computes");
    assert_eq!(characteristic.to_string(), "t^4 - t^2 + 1/16");
    assert_eq!(minimal.to_string(), "t^4 - t^2 + 1/16");
}

#[test]
fn minimal_polynomial_records_the_first_krylov_dependence() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the ring builds");
    let x = ring.parse_polynomial("x").expect("x parses");
    let ideal = ring
        .ideal([
            ring.parse_polynomial("x^2").expect("x^2 parses"),
            ring.parse_polynomial("y^2").expect("y^2 parses"),
        ])
        .expect("the generators belong to the ring");
    let basis = ideal
        .groebner_basis(RationalOptions::new())
        .expect("the rational basis computes");
    let quotient = basis
        .finite_quotient(sylvester::Budget::new())
        .expect("the quotient is finite");

    let x_coordinates = quotient
        .coordinates(&x, sylvester::Budget::new())
        .expect("the x residue coordinates read");
    assert!(x_coordinates[0].is_none(), "x is not a scalar residue");
    assert!(x_coordinates.iter().skip(1).any(Option::is_some));
    assert!(
        quotient
            .multiply(&x, &x, sylvester::Budget::new())
            .expect("x squared reduces")
            .is_zero()
    );

    let characteristic = quotient
        .characteristic_polynomial(&x, sylvester::Budget::new())
        .expect("the characteristic polynomial computes");
    let minimal = quotient
        .minimal_polynomial(&x, sylvester::Budget::new())
        .expect("the minimal polynomial computes");
    assert_eq!(characteristic.to_string(), "t^4");
    assert_eq!(minimal.to_string(), "t^2");
    assert!(minimal.degree() < characteristic.degree());
}

#[test]
fn zero_quotient_has_unit_spectra() {
    let ring = PolynomialRing::prime_field(5, ["x"]).expect("5 is prime");
    let x = ring.parse_polynomial("x").expect("x parses");
    let ideal = ring.ideal([ring.one()]).expect("one belongs to the ring");
    let basis = ideal
        .groebner_basis(ComputeOptions::new())
        .expect("the basis computes");
    let quotient = basis
        .finite_quotient(sylvester::Budget::new())
        .expect("the quotient is finite");

    let matrix = quotient
        .multiplication_matrix(&x, sylvester::Budget::new())
        .expect("the zero quotient has a matrix");
    assert_eq!(matrix.dimension(), 0);
    assert_eq!(matrix.entries(), &[]);
    let characteristic: UnivariatePolynomial = quotient
        .characteristic_polynomial(&x, sylvester::Budget::new())
        .expect("the characteristic polynomial computes");
    let minimal = quotient
        .minimal_polynomial(&x, sylvester::Budget::new())
        .expect("the minimal polynomial computes");
    assert_eq!(characteristic.to_string(), "1");
    assert_eq!(minimal.to_string(), "1");
}

#[test]
fn constant_ring_has_linear_spectra() {
    let ring = PolynomialRing::prime_field(5, std::iter::empty::<&str>()).expect("the ring builds");
    let basis = GroebnerBasis::<PrimeField>::from_polynomials(&ring, Vec::new(), Budget::new())
        .expect("the zero basis computes");
    let quotient = basis
        .finite_quotient(Budget::new())
        .expect("the quotient is finite");
    let one = ring.one();

    let matrix = quotient
        .multiplication_matrix(&one, Budget::new())
        .expect("the matrix computes");
    assert_eq!(matrix.dimension(), 1);
    assert_eq!(matrix.to_string(), "[[1]]");
    assert_eq!(
        quotient
            .characteristic_polynomial(&one, Budget::new())
            .expect("the characteristic polynomial computes")
            .to_string(),
        "t + 4"
    );
    assert_eq!(
        quotient
            .minimal_polynomial(&one, Budget::new())
            .expect("the minimal polynomial computes")
            .to_string(),
        "t + 4"
    );
}
