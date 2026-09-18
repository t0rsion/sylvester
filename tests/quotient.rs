//! Checks for finite quotient construction and residue arithmetic.

use sylvester::{
    Budget, ComputeOptions, Felt, FiniteQuotient, GroebnerBasis, Polynomial, PolynomialRing,
    QuotientError,
};

fn checked_basis(ring: &PolynomialRing, polynomials: &[&str]) -> GroebnerBasis {
    let polynomials = polynomials
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the polynomial parses"))
        .collect();
    GroebnerBasis::<sylvester::PrimeField>::from_polynomials(ring, polynomials, Budget::new())
        .expect("the list is a Groebner basis")
}

fn values(coordinates: Vec<Option<Felt>>) -> Vec<Option<u64>> {
    coordinates
        .into_iter()
        .map(|coefficient| coefficient.map(Felt::value))
        .collect()
}

fn first_reduce_memory_limit(quotient: &FiniteQuotient, polynomial: &Polynomial) -> usize {
    let mut high = 1usize;
    while quotient
        .reduce(polynomial, Budget::new().memory_limit(high))
        .is_err()
    {
        high = high
            .checked_mul(2)
            .expect("the operation fits a memory limit");
    }
    let mut low = 0usize;
    while low < high {
        let middle = low + (high - low) / 2;
        if quotient
            .reduce(polynomial, Budget::new().memory_limit(middle))
            .is_ok()
        {
            high = middle;
        } else {
            low = middle + 1;
        }
    }
    low
}

#[test]
fn staircase_coordinates_and_residue_arithmetic_follow_grevlex() {
    let ring = PolynomialRing::prime_field(101, ["x", "y"]).expect("the ring builds");
    let basis = checked_basis(&ring, &["y^3", "x^2"]);
    assert!(basis.is_zero_dimensional());
    let quotient = basis
        .finite_quotient(Budget::new())
        .expect("the quotient is finite");

    assert_eq!(
        quotient.standard_monomials(),
        &[
            vec![0, 0],
            vec![0, 1],
            vec![1, 0],
            vec![0, 2],
            vec![1, 1],
            vec![1, 2],
        ]
    );
    assert_eq!(quotient.vector_space_dimension(), 6);
    assert_eq!(quotient.source_basis(), &basis);

    let residue = ring
        .parse_polynomial("x^2 + x*y^2 + 2*y")
        .expect("the polynomial parses");
    assert_eq!(
        quotient
            .reduce(&residue, Budget::new())
            .unwrap()
            .to_string(),
        "x*y^2 + 2*y"
    );
    assert_eq!(
        values(quotient.coordinates(&residue, Budget::new()).unwrap()),
        vec![None, Some(2), None, None, None, Some(1)]
    );

    let x = ring.parse_polynomial("x").expect("x parses");
    let y = ring.parse_polynomial("y").expect("y parses");
    assert_eq!(
        quotient
            .multiply(&x, &y, Budget::new())
            .unwrap()
            .to_string(),
        "x*y"
    );
    assert_eq!(
        quotient.multiply(&y, &x, Budget::new()).unwrap(),
        quotient.multiply(&x, &y, Budget::new()).unwrap()
    );
    assert!(quotient.pow(&x, 3, Budget::new()).unwrap().is_zero());
}

#[test]
fn multiplication_matrix_columns_are_normal_form_coordinates() {
    let ring = PolynomialRing::prime_field(101, ["x", "y"]).expect("the ring builds");
    let basis = checked_basis(&ring, &["x^2", "y^2"]);
    let quotient = basis
        .finite_quotient(Budget::new())
        .expect("the quotient is finite");
    let polynomial = ring
        .parse_polynomial("x + y")
        .expect("the polynomial parses");
    let matrix = quotient
        .multiplication_matrix(&polynomial, Budget::new())
        .expect("the matrix builds");

    assert_eq!(matrix.dimension(), quotient.vector_space_dimension());
    for (column, exponents) in quotient.standard_monomials().iter().enumerate() {
        let basis_monomial = ring
            .polynomial([(1, exponents.clone())])
            .expect("the basis monomial belongs to the ring");
        let expected = quotient
            .coordinates(
                &quotient
                    .multiply(&polynomial, &basis_monomial, Budget::new())
                    .expect("the product reduces"),
                Budget::new(),
            )
            .expect("the coordinates read");
        for (row, coefficient) in expected.into_iter().enumerate() {
            assert_eq!(matrix.entry(row, column).copied(), coefficient);
        }
    }
}

#[test]
fn unit_and_zero_variable_quotients_have_their_edge_dimensions() {
    let ring = PolynomialRing::prime_field(101, ["x"]).expect("the ring builds");
    let unit_basis = ring
        .ideal([ring.one()])
        .expect("the generator belongs to the ring")
        .groebner_basis(ComputeOptions::new())
        .expect("the basis computes");
    assert!(unit_basis.is_zero_dimensional());
    let unit = unit_basis
        .finite_quotient(Budget::new())
        .expect("the unit quotient is finite");
    assert_eq!(unit.vector_space_dimension(), 0);
    assert!(unit.standard_monomials().is_empty());

    let constants = PolynomialRing::prime_field(101, std::iter::empty::<&str>())
        .expect("the constant ring builds");
    let zero_basis = checked_basis(&constants, &[]);
    assert!(zero_basis.is_zero_dimensional());
    let zero = zero_basis
        .finite_quotient(Budget::new())
        .expect("the zero-variable quotient is finite");
    assert_eq!(zero.vector_space_dimension(), 1);
    assert_eq!(zero.standard_monomials(), &[Vec::<u16>::new()]);
    assert_eq!(
        values(zero.coordinates(&constants.one(), Budget::new()).unwrap()),
        vec![Some(1)]
    );
}

#[test]
fn positive_dimension_ring_mismatch_and_budget_are_typed() {
    let ring = PolynomialRing::prime_field(101, ["x", "y"]).expect("the ring builds");
    let basis = checked_basis(&ring, &["x"]);
    assert!(!basis.is_zero_dimensional());
    assert!(matches!(
        basis.finite_quotient(Budget::new()),
        Err(QuotientError::NotFinite)
    ));

    let finite_basis = checked_basis(&ring, &["x^2", "y^2"]);
    assert!(matches!(
        finite_basis.finite_quotient(Budget::new().memory_limit(0)),
        Err(QuotientError::MemoryLimitExceeded)
    ));
    let quotient = finite_basis
        .finite_quotient(Budget::new())
        .expect("the quotient is finite");
    let other_ring = PolynomialRing::prime_field(101, ["x", "z"]).expect("the other ring builds");
    let other_polynomial = other_ring.parse_polynomial("x").expect("x parses");
    assert_eq!(
        quotient.reduce(&other_polynomial, Budget::new()),
        Err(QuotientError::RingMismatch)
    );
}

#[test]
fn quotient_budget_counts_spare_polynomial_capacity() {
    let ring = PolynomialRing::prime_field(101, ["x"]).expect("the ring builds");
    let basis = checked_basis(&ring, &["x^2"]);
    let quotient = basis
        .finite_quotient(Budget::new())
        .expect("the quotient is finite");
    let tail = ring
        .polynomial((1..=64).map(|exponent| (1i64, vec![exponent])))
        .expect("the cancellation tail builds");
    let left = ring
        .one()
        .try_add(&tail, Budget::new())
        .expect("the sum builds");
    let spare = left
        .try_sub(&tail, Budget::new())
        .expect("the cancellation builds");
    let compact = ring.one();
    assert_eq!(spare, compact);
    assert!(spare.estimated_heap_bytes() > compact.estimated_heap_bytes());

    let compact_limit = first_reduce_memory_limit(&quotient, &compact);
    assert_eq!(
        quotient.reduce(&spare, Budget::new().memory_limit(compact_limit)),
        Err(QuotientError::MemoryLimitExceeded)
    );
}
