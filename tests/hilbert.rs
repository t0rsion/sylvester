//! The Hilbert series, the dimension, and the multiplicity.
//!
//! Every reference value in this file comes from Singular 4.4.1, one call
//! per system:
//!
//! ```text
//! ring r = 32003, (x1,x2,x3,x4), dp;
//! ideal I = <generators>;
//! ideal G = std(I);
//! hilb(G, 1); dim(G); mult(G); vdim(G);
//! ```
//!
//! `hilb(G, 1)` prints the numerator over `(1 - t)^n`, lowest degree
//! first, and closes the list with one 0 that is not a coefficient. The
//! constants below are that output with the closing 0 dropped. `dim` is
//! the Krull dimension of the quotient, `mult` its multiplicity, and
//! `vdim` the dimension of the quotient over `k`, which is -1 when the
//! quotient is not finite over `k`.
//!
//! The hand-checked cases carry their derivation next to the test.

use std::time::Duration;

use num_bigint::BigInt;
use proptest::prelude::*;
use sylvester::{
    Budget, ComputeOptions, Domain, GroebnerBasis, HilbertError, HilbertSeries, PolynomialRing,
    PrimeField, RationalOptions, Rationals,
};

const MODULUS: u64 = 32003;

/// cyclic-4 in four variables, the generators of `benchmarks/gb-comparison`.
const CYCLIC_4: [&str; 4] = [
    "x1 + x2 + x3 + x4",
    "x1*x2 + x2*x3 + x1*x4 + x3*x4",
    "x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4",
    "x1*x2*x3*x4 - 1",
];

/// The numerator of cyclic-4 over `(1 - t)^4`.
///
/// Singular: `1,-1,-1,0,1,0,1,0,-2,1,0`, `dim(G) = 1`, `mult(G) = 4`,
/// `vdim(G) = -1`. The same values come out over `Q` and over `F_32003`.
const CYCLIC_4_NUMERATOR: [i64; 10] = [1, -1, -1, 0, 1, 0, 1, 0, -2, 1];

/// katsura-4 in five variables, the POSSO definition.
const KATSURA_4: [&str; 5] = [
    "x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 - 1",
    "x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 - x1",
    "2*x1*x2 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 - x2",
    "x2^2 + 2*x1*x3 + 2*x2*x4 + 2*x3*x5 - x3",
    "2*x2*x3 + 2*x1*x4 + 2*x2*x5 - x4",
];

/// The numerator of katsura-4 over `(1 - t)^5`.
///
/// Singular: `1,-1,-4,4,6,-6,-4,4,1,-1,0`, `dim(G) = 0`,
/// `mult(G) = vdim(G) = 16`. The numerator is `(1 - t)^5 * (1 + t)^4`, so
/// the graded pieces of the quotient have dimensions 1, 4, 6, 4, 1 and
/// then 0.
const KATSURA_4_NUMERATOR: [i64; 10] = [1, -1, -4, 4, 6, -6, -4, 4, 1, -1];

/// cyclic-6 in six variables.
const CYCLIC_6: [&str; 6] = [
    "x1 + x2 + x3 + x4 + x5 + x6",
    "x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x6 + x6*x1",
    "x1*x2*x3 + x2*x3*x4 + x3*x4*x5 + x4*x5*x6 + x5*x6*x1 + x6*x1*x2",
    "x1*x2*x3*x4 + x2*x3*x4*x5 + x3*x4*x5*x6 + x4*x5*x6*x1 + x5*x6*x1*x2 + x6*x1*x2*x3",
    "x1*x2*x3*x4*x5 + x2*x3*x4*x5*x6 + x3*x4*x5*x6*x1 + x4*x5*x6*x1*x2 + x5*x6*x1*x2*x3 + x6*x1*x2*x3*x4",
    "x1*x2*x3*x4*x5*x6 - 1",
];

/// The numerator of cyclic-6 over `(1 - t)^6`.
///
/// Singular: `1,-1,-1,-4,1,34,-60,37,-5,3,-1,-34,59,-36,6,1,0`,
/// `dim(G) = 0`, `mult(G) = vdim(G) = 156`.
const CYCLIC_6_NUMERATOR: [i64; 16] =
    [1, -1, -1, -4, 1, 34, -60, 37, -5, 3, -1, -34, 59, -36, 6, 1];

/// The reduced basis of cyclic-3 over `Q`, largest leading monomial first.
///
/// The leading monomials are `z^3`, `y^2`, and `x`, one pure power per
/// variable, so the numerator is `(1 - t)*(1 - t^2)*(1 - t^3)`.
/// Singular: `1,-1,-1,0,1,1,-1,0`, `dim(G) = 0`,
/// `mult(G) = vdim(G) = 6`.
const CYCLIC_3_BASIS: [&str; 3] = ["z^3 - 1", "y^2 + y*z + z^2", "x + y + z"];

const CYCLIC_3_NUMERATOR: [i64; 7] = [1, -1, -1, 0, 1, 1, -1];

fn prime_ring(variables: &[&str]) -> PolynomialRing {
    PolynomialRing::prime_field(MODULUS, variables.iter().copied()).expect("the modulus is prime")
}

fn parse_all<D: Domain>(ring: &PolynomialRing<D>, texts: &[&str]) -> Vec<sylvester::Polynomial<D>> {
    texts
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the syntax holds"))
        .collect()
}

/// The basis the F4 engine computes over `F_32003`.
fn prime_basis(ring: &PolynomialRing, generators: &[&str]) -> GroebnerBasis {
    ring.ideal(parse_all(ring, generators))
        .expect("the generators share the ring")
        .groebner_basis(ComputeOptions::new())
        .expect("the system fits every budget")
}

fn series(basis: &GroebnerBasis) -> HilbertSeries {
    basis
        .hilbert_series(Budget::new())
        .expect("the recursion fits every budget")
}

fn numerator(coefficients: &[i64]) -> Vec<BigInt> {
    coefficients.iter().copied().map(BigInt::from).collect()
}

/// The number of monomials of each degree up to `max_degree` that no
/// leading monomial of `basis` divides.
///
/// These are the standard monomials, and they are a basis over `k` of the
/// quotient by the leading monomial ideal. The count is the independent
/// value the series is checked against.
fn standard_monomials<D: Domain>(basis: &GroebnerBasis<D>, max_degree: u32) -> Vec<u64> {
    let leading: Vec<Vec<u16>> = basis
        .iter()
        .filter_map(|f| f.leading_term().map(|(_, exps)| exps.to_vec()))
        .collect();
    let nvars = basis.ring().nvars();
    let mut counts = vec![0u64; max_degree as usize + 1];
    let mut exps = vec![0u16; nvars];
    walk(&mut exps, 0, max_degree, &leading, &mut counts);
    counts
}

/// Walk every exponent vector of degree at most `left` and count the ones
/// no leading monomial divides.
fn walk(exps: &mut Vec<u16>, index: usize, left: u32, leading: &[Vec<u16>], counts: &mut Vec<u64>) {
    if index == exps.len() {
        let degree: u32 = exps.iter().map(|&e| u32::from(e)).sum();
        let divisible = leading
            .iter()
            .any(|lead| lead.iter().zip(exps.iter()).all(|(a, b)| a <= b));
        if !divisible {
            counts[degree as usize] += 1;
        }
        return;
    }
    for exp in 0..=left {
        exps[index] = exp as u16;
        walk(exps, index + 1, left - exp, leading, counts);
    }
    exps[index] = 0;
}

/// `(x^2, x*y)` in `k[x, y]`.
///
/// By hand: the standard monomials are 1, then `x` and `y`, then `y^d`
/// alone in every degree above 1, so the series is
/// `1 + 2*t + t^2 + t^3 + ...`. Over `(1 - t)^2` the numerator is
/// `1 - 2*t^2 + t^3`. The recursion reaches the same value: the pivot is
/// `x`, `L + (x) = (x)` gives `1 - t`, `L : x = (x, y)` gives
/// `(1 - t)^2`, and `(1 - t) + t*(1 - t)^2 = 1 - 2*t^2 + t^3`.
/// Singular: `1,0,-2,1,0`, `dim(G) = 1`, `mult(G) = 1`.
#[test]
fn the_series_of_the_hand_checked_monomial_ideal_holds() {
    let ring = prime_ring(&["x", "y"]);
    let basis = prime_basis(&ring, &["x^2", "x*y"]);
    let series = series(&basis);

    assert_eq!(series.numerator(), numerator(&[1, 0, -2, 1]));
    assert_eq!(series.denominator_power(), 2);
    assert_eq!(series.dimension(), Some(1));
    assert_eq!(series.multiplicity(), Some(BigInt::from(1)));
    assert_eq!(series.to_string(), "(1 - 2*t^2 + t^3)/(1 - t)^2");
    assert_eq!(
        basis.krull_dimension(Budget::new()).expect("it fits"),
        Some(1)
    );
}

#[test]
fn the_series_of_the_maximal_ideal_holds() {
    let plane = prime_ring(&["x", "y"]);
    let maximal = series(&prime_basis(&plane, &["x", "y"]));
    assert_eq!(maximal.numerator(), numerator(&[1, -2, 1]));
    assert_eq!(maximal.dimension(), Some(0));
    assert_eq!(maximal.multiplicity(), Some(BigInt::from(1)));
}

#[test]
fn the_series_of_a_monomial_hypersurface_holds() {
    let space = prime_ring(&["x", "y", "z"]);
    let surface = series(&prime_basis(&space, &["x*y"]));
    assert_eq!(surface.numerator(), numerator(&[1, 0, -1]));
    assert_eq!(surface.dimension(), Some(2));
    assert_eq!(surface.multiplicity(), Some(BigInt::from(2)));
    assert_eq!(surface.to_string(), "(1 - t^2)/(1 - t)^3");
}

/// Without the pivot restriction in design section 6.2, `L + (x) = L`.
#[test]
fn the_series_of_x_and_yz_holds() {
    let space = prime_ring(&["x", "y", "z"]);
    let curve = series(&prime_basis(&space, &["x", "y*z"]));
    assert_eq!(curve.numerator(), numerator(&[1, -1, -1, 1]));
    assert_eq!(curve.dimension(), Some(1));
    assert_eq!(curve.multiplicity(), Some(BigInt::from(2)));
}

/// The principal ideal `(x^5)` in `k[x, y]`.
///
/// Singular: `1,0,0,0,0,-1,0`, `dim(G) = 1`, `mult(G) = 5`.
#[test]
fn the_series_of_a_principal_ideal_is_one_minus_the_generator() {
    let ring = prime_ring(&["x", "y"]);
    let series = series(&prime_basis(&ring, &["x^5"]));

    assert_eq!(series.numerator(), numerator(&[1, 0, 0, 0, 0, -1]));
    assert_eq!(series.dimension(), Some(1));
    assert_eq!(series.multiplicity(), Some(BigInt::from(5)));
    assert_eq!(series.to_string(), "(1 - t^5)/(1 - t)^2");
}

#[test]
fn the_zero_ideal_holds_its_conventions() {
    let ring = prime_ring(&["x", "y"]);
    let zero = GroebnerBasis::<PrimeField>::from_polynomials(&ring, Vec::new(), Budget::new())
        .expect("the empty list is the basis of the zero ideal");
    let zero = series(&zero);
    assert_eq!(zero.numerator(), numerator(&[1]));
    assert_eq!(zero.denominator_power(), 2);
    assert_eq!(zero.dimension(), Some(2));
    assert_eq!(zero.multiplicity(), Some(BigInt::from(1)));
    assert_eq!(zero.to_string(), "(1)/(1 - t)^2");
    for degree in 0..6u32 {
        assert_eq!(zero.coefficient(degree), BigInt::from(degree + 1));
    }
}

#[test]
fn the_unit_ideal_holds_its_conventions() {
    let ring = prime_ring(&["x", "y"]);
    let one = GroebnerBasis::<PrimeField>::from_polynomials(
        &ring,
        parse_all(&ring, &["1"]),
        Budget::new(),
    )
    .expect("the list holding 1 is the basis of the unit ideal");
    let unit = series(&one);
    assert!(unit.numerator().is_empty());
    assert_eq!(unit.dimension(), None);
    assert_eq!(unit.multiplicity(), None);
    assert_eq!(unit.to_string(), "(0)/(1 - t)^2");
    assert_eq!(unit.coefficient(0), BigInt::from(0));
    assert_eq!(one.krull_dimension(Budget::new()).expect("it fits"), None);
}

/// cyclic-4 over `F_32003` against the Singular numerator.
#[test]
fn the_series_of_cyclic_4_is_the_singular_numerator() {
    let ring = prime_ring(&["x1", "x2", "x3", "x4"]);
    let basis = prime_basis(&ring, &CYCLIC_4);
    let series = series(&basis);

    assert_eq!(series.numerator(), numerator(&CYCLIC_4_NUMERATOR));
    assert_eq!(series.denominator_power(), 4);
    assert_eq!(series.dimension(), Some(1));
    assert_eq!(series.multiplicity(), Some(BigInt::from(4)));
}

/// katsura-4 over `F_32003` against the Singular numerator.
///
/// The quotient is finite over `k`: `vdim(G) = 16` is the sum of the
/// graded dimensions 1, 4, 6, 4, 1, which is also `mult(G)`.
#[test]
fn the_series_of_katsura_4_is_the_singular_numerator() {
    let ring = prime_ring(&["x1", "x2", "x3", "x4", "x5"]);
    let basis = prime_basis(&ring, &KATSURA_4);
    let series = series(&basis);

    assert_eq!(series.numerator(), numerator(&KATSURA_4_NUMERATOR));
    assert_eq!(series.dimension(), Some(0));
    assert_eq!(series.multiplicity(), Some(BigInt::from(16)));

    let graded: Vec<BigInt> = (0..6).map(|degree| series.coefficient(degree)).collect();
    assert_eq!(graded, numerator(&[1, 4, 6, 4, 1, 0]));
    let total: BigInt = graded.into_iter().sum();
    assert_eq!(total, BigInt::from(16));
}

/// cyclic-6 over `F_32003` against the Singular numerator.
///
/// The system is the largest one in this file. Its quotient is finite over
/// `k`, and the multiplicity 156 is the number of solutions with
/// multiplicity.
#[test]
fn the_series_of_cyclic_6_is_the_singular_numerator() {
    let ring = prime_ring(&["x1", "x2", "x3", "x4", "x5", "x6"]);
    let series = series(&prime_basis(&ring, &CYCLIC_6));

    assert_eq!(series.numerator(), numerator(&CYCLIC_6_NUMERATOR));
    assert_eq!(series.dimension(), Some(0));
    assert_eq!(series.multiplicity(), Some(BigInt::from(156)));
}

/// The series of one ideal is the same over `F_p` and over `Q`.
///
/// The computation reads leading monomials, and cyclic-4 has the same
/// leading monomials over both domains, so the two series agree. The
/// rational basis comes from the multimodular engine, the prime-field one
/// from F4.
#[test]
fn the_series_does_not_depend_on_the_domain() {
    let prime = series(&prime_basis(
        &prime_ring(&["x1", "x2", "x3", "x4"]),
        &CYCLIC_4,
    ));

    let ring = PolynomialRing::rationals(["x1", "x2", "x3", "x4"]).expect("the names hold");
    let rational = ring
        .ideal(parse_all(&ring, &CYCLIC_4))
        .expect("the generators share the ring")
        .groebner_basis(RationalOptions::new())
        .expect("cyclic-4 lifts")
        .hilbert_series(Budget::new())
        .expect("the recursion fits every budget");

    assert_eq!(rational.numerator(), numerator(&CYCLIC_4_NUMERATOR));
    assert_eq!(rational, prime);
}

/// A checked basis over `Q` carries the series of its own leading ideal.
#[test]
fn a_checked_rational_basis_has_the_series_of_cyclic_3() {
    let ring = PolynomialRing::rationals(["x", "y", "z"]).expect("the names hold");
    let basis = GroebnerBasis::<Rationals>::from_polynomials(
        &ring,
        parse_all(&ring, &CYCLIC_3_BASIS),
        Budget::new(),
    )
    .expect("the list is the reduced basis of cyclic-3");
    let series = basis
        .hilbert_series(Budget::new())
        .expect("the recursion fits every budget");

    assert_eq!(series.numerator(), numerator(&CYCLIC_3_NUMERATOR));
    assert_eq!(series.dimension(), Some(0));
    assert_eq!(series.multiplicity(), Some(BigInt::from(6)));
    assert_eq!(
        basis.krull_dimension(Budget::new()).expect("it fits"),
        Some(0)
    );
}

/// A homogeneous ideal, where the series is the ideal's own Hilbert
/// series.
///
/// `I = (x^2 - y*z, x*y - z^2)` is homogeneous, so its reduced basis is,
/// and the graded dimensions of `k[x,y,z] / I` are the graded dimensions
/// of the quotient by the leading ideal. The test counts the standard
/// monomials of each degree and compares them with the coefficients.
/// Singular: `1,0,-2,0,1,0`, `dim(G) = 1`, `mult(G) = 4`.
#[test]
fn a_homogeneous_ideal_carries_its_own_hilbert_series() {
    let ring = prime_ring(&["x", "y", "z"]);
    let ideal = ring
        .ideal(parse_all(&ring, &["x^2 - y*z", "x*y - z^2"]))
        .expect("the generators share the ring");
    assert!(ideal.has_homogeneous_generators());
    let basis = ideal
        .groebner_basis(ComputeOptions::new())
        .expect("the system fits every budget");
    assert!(basis.is_homogeneous());

    let series = series(&basis);
    assert_eq!(series.numerator(), numerator(&[1, 0, -2, 0, 1]));
    assert_eq!(series.dimension(), Some(1));
    assert_eq!(series.multiplicity(), Some(BigInt::from(4)));

    let counts = standard_monomials(&basis, 8);
    for (degree, count) in counts.iter().enumerate() {
        assert_eq!(
            series.coefficient(degree as u32),
            BigInt::from(*count),
            "degree {degree}"
        );
    }
}

/// An inhomogeneous ideal, where the coefficients are the first difference
/// of the affine Hilbert function and not that function.
///
/// katsura-4 has inhomogeneous generators and an inhomogeneous basis. The
/// running sum of the coefficients up to degree `d` is the affine Hilbert
/// function `H(d)`, which the count of standard monomials of degree at
/// most `d` gives independently.
#[test]
fn the_coefficients_are_the_first_difference_of_the_affine_function() {
    let ring = prime_ring(&["x1", "x2", "x3", "x4", "x5"]);
    let ideal = ring
        .ideal(parse_all(&ring, &KATSURA_4))
        .expect("the generators share the ring");
    assert!(!ideal.has_homogeneous_generators());
    let basis = ideal
        .groebner_basis(ComputeOptions::new())
        .expect("the system fits every budget");
    assert!(!basis.is_homogeneous());

    let series = series(&basis);
    let counts = standard_monomials(&basis, 8);
    let mut running = BigInt::from(0);
    for (degree, count) in counts.iter().enumerate() {
        running += series.coefficient(degree as u32);
        assert_eq!(
            running,
            BigInt::from(counts[..=degree].iter().sum::<u64>()),
            "degree {degree}"
        );
        assert_eq!(series.coefficient(degree as u32), BigInt::from(*count));
    }
}

/// The zero ideal in one variable: the coefficients are 1 in every degree
/// while the affine Hilbert function is `d + 1`. Section 6.1 of
/// `docs/rational-design.md` names this case.
#[test]
fn the_coefficients_are_not_the_affine_hilbert_function() {
    let ring = prime_ring(&["x"]);
    let basis = GroebnerBasis::<PrimeField>::from_polynomials(&ring, Vec::new(), Budget::new())
        .expect("the empty list is the basis of the zero ideal");
    let series = series(&basis);

    let mut running = BigInt::from(0);
    for degree in 0..6u32 {
        assert_eq!(series.coefficient(degree), BigInt::from(1));
        running += series.coefficient(degree);
        assert_eq!(running, BigInt::from(degree + 1));
    }
}

/// An exhausted budget is typed, and it says nothing about the ideal.
#[test]
fn an_exhausted_budget_is_typed() {
    let ring = prime_ring(&["x1", "x2", "x3", "x4"]);
    let basis = prime_basis(&ring, &CYCLIC_4);

    assert_eq!(
        basis.hilbert_series(Budget::new().timeout(Duration::ZERO)),
        Err(HilbertError::Timeout)
    );
    assert_eq!(
        basis.hilbert_series(Budget::new().memory_limit(0)),
        Err(HilbertError::MemoryLimitExceeded)
    );
    assert_eq!(
        basis.krull_dimension(Budget::new().memory_limit(0)),
        Err(HilbertError::MemoryLimitExceeded)
    );
}

/// A base case charges its vector before it allocates it.
///
/// The single generator `x^65535` gives the numerator `1 - t^65535`, a
/// dense vector of 65,536 coefficients. The limit below is smaller than
/// that vector, so the call reports the limit instead of allocating.
#[test]
fn a_base_case_charges_the_vector_it_allocates() {
    let ring = prime_ring(&["x", "y"]);
    let basis = GroebnerBasis::<PrimeField>::from_polynomials(
        &ring,
        vec![
            ring.polynomial([(1, [u16::MAX, 0].as_slice())])
                .expect("the width holds"),
        ],
        Budget::new(),
    )
    .expect("one monic monomial is a reduced basis");

    assert_eq!(
        basis.hilbert_series(Budget::new().memory_limit(1_000)),
        Err(HilbertError::MemoryLimitExceeded)
    );
}

/// A product of pure powers charges every factor and every product.
///
/// The four generators are pure powers of distinct variables, so the
/// numerator is the product of the four factors `1 - t^16384`. The
/// product holds 65,533 coefficients, past the limit below.
#[test]
fn a_product_of_pure_powers_charges_its_factors() {
    const DEGREE: u16 = 16_384;
    let ring = prime_ring(&["x1", "x2", "x3", "x4"]);
    let mut polynomials = Vec::new();
    for index in 0..4 {
        let mut exps = [0u16; 4];
        exps[index] = DEGREE;
        polynomials.push(
            ring.polynomial([(1, exps.as_slice())])
                .expect("the width holds"),
        );
    }
    let basis = GroebnerBasis::<PrimeField>::from_polynomials(&ring, polynomials, Budget::new())
        .expect("pure powers of distinct variables are a reduced basis");

    assert_eq!(
        basis.hilbert_series(Budget::new().memory_limit(1_000)),
        Err(HilbertError::MemoryLimitExceeded)
    );
    assert_eq!(
        basis.hilbert_series(Budget::new().timeout(Duration::ZERO)),
        Err(HilbertError::Timeout)
    );
}

/// Two calls on one basis return one value.
#[test]
fn the_series_is_deterministic() {
    let ring = prime_ring(&["x1", "x2", "x3", "x4", "x5"]);
    let basis = prime_basis(&ring, &KATSURA_4);
    let first = series(&basis);
    let second = series(&basis);

    assert_eq!(first, second);
    assert_eq!(first.to_string(), second.to_string());
}

/// A monomial ideal from `exponents`, as a basis over `F_32003`.
fn monomial_ideal(variables: &[&str], exponents: &[Vec<u16>]) -> GroebnerBasis {
    let ring = prime_ring(variables);
    let generators = exponents.iter().map(|exps| {
        ring.polynomial([(1, exps.as_slice())])
            .expect("the width holds")
    });
    ring.ideal(generators)
        .expect("the generators share the ring")
        .groebner_basis(ComputeOptions::new())
        .expect("a monomial ideal fits every budget")
}

proptest! {
    /// The coefficients count the standard monomials of random monomial
    /// ideals in three and four variables, up to degree 8.
    #[test]
    fn the_coefficients_count_the_standard_monomials_in_three_variables(
        exponents in prop::collection::vec(prop::collection::vec(0u16..4, 3), 1..6)
    ) {
        let basis = monomial_ideal(&["x", "y", "z"], &exponents);
        let series = basis.hilbert_series(Budget::new()).expect("it fits");
        for (degree, count) in standard_monomials(&basis, 8).iter().enumerate() {
            prop_assert_eq!(series.coefficient(degree as u32), BigInt::from(*count));
        }
    }

    #[test]
    fn the_coefficients_count_the_standard_monomials_in_four_variables(
        exponents in prop::collection::vec(prop::collection::vec(0u16..3, 4), 1..6)
    ) {
        let basis = monomial_ideal(&["w", "x", "y", "z"], &exponents);
        let series = basis.hilbert_series(Budget::new()).expect("it fits");
        for (degree, count) in standard_monomials(&basis, 8).iter().enumerate() {
            prop_assert_eq!(series.coefficient(degree as u32), BigInt::from(*count));
        }
    }
}

/// The colon branch of a chain returns the closed form.
///
/// `L = (x^d*y, z)` in `k[x, y, z]` gives `R / L = k[x, y] / (x^d*y)`,
/// whose series is `(1 - t^(d+1)) / (1 - t)^2`. Over `(1 - t)^3` the
/// numerator is `(1 - t^(d+1))*(1 - t) = 1 - t - t^(d+1) + t^(d+2)`. The
/// recursion reaches it one colon at a time: the pivot is `x`, and
/// `L : x = (x^(d-1)*y, z)`, so the chain is `d` nodes deep.
#[test]
fn a_colon_chain_returns_the_closed_form() {
    const DEGREE: usize = 300;
    let basis = chain_basis(DEGREE as u16);
    let series = basis
        .hilbert_series(Budget::new())
        .expect("300 nodes fit every budget");

    let mut expected = vec![BigInt::from(0); DEGREE + 3];
    expected[0] = BigInt::from(1);
    expected[1] = BigInt::from(-1);
    expected[DEGREE + 1] = BigInt::from(-1);
    expected[DEGREE + 2] = BigInt::from(1);
    assert_eq!(series.numerator(), expected.as_slice());
    assert_eq!(series.denominator_power(), 3);
    assert_eq!(series.dimension(), Some(1));
    assert_eq!(series.multiplicity(), Some(BigInt::from(DEGREE + 1)));
}

/// A colon chain of 65,535 nodes reports its budget and does not abort.
///
/// The chain is deeper than the calling stack holds, so the worklist runs
/// on the heap. The descent charges about 128 bytes of generator sets and
/// one frame per node, which stays under the limit below. The numerators
/// the chain builds on the way back hold 65,538 coefficients each, so the
/// memo passes the limit after a few of them and the call reports it.
#[test]
fn a_deep_colon_chain_reports_its_budget() {
    let basis = chain_basis(u16::MAX);
    assert_eq!(
        basis.hilbert_series(Budget::new().memory_limit(64 << 20)),
        Err(HilbertError::MemoryLimitExceeded)
    );
}

/// The basis `[x^degree*y, z]`, whose colon chain is `degree` nodes deep.
fn chain_basis(degree: u16) -> GroebnerBasis {
    let ring = prime_ring(&["x", "y", "z"]);
    let polynomials = vec![
        ring.polynomial([(1, [degree, 1, 0].as_slice())])
            .expect("the width holds"),
        ring.polynomial([(1, [0u16, 0, 1].as_slice())])
            .expect("the width holds"),
    ];
    GroebnerBasis::<PrimeField>::from_polynomials(&ring, polynomials, Budget::new())
        .expect("two monomials with a coprime pair are a reduced basis")
}
