//! Normal form, membership, and the checked basis constructor.
//!
//! A basis here comes from one of two places: the F4 engine over `F_p`,
//! or the checked constructor over a list this file states. The list is
//! cyclic-3, verified by hand. Substituting `x = -(y + z)` into
//! `x*y + y*z + z*x` gives `-(y^2 + y*z + z^2)`, and reducing `x*y*z - 1`
//! by both gives `z^3 - 1`. Those three, largest leading monomial first,
//! are the reduced basis under grevlex.

use std::time::Duration;

use num_rational::BigRational;
use proptest::prelude::*;
use sylvester::{
    BasisError, Budget, ComputeOptions, GroebnerBasis, NormalFormError, Polynomial, PolynomialRing,
    PrimeField, Rationals,
};

const MODULUS: u64 = 32003;

/// The generators of cyclic-3 as text, in three variables.
const CYCLIC_3: [&str; 3] = ["x + y + z", "x*y + y*z + z*x", "x*y*z - 1"];

/// The reduced Gröbner basis of cyclic-3 under grevlex, largest leading
/// monomial first.
const CYCLIC_3_BASIS: [&str; 3] = ["z^3 - 1", "y^2 + y*z + z^2", "x + y + z"];

fn prime_ring(variables: &[&str]) -> PolynomialRing {
    PolynomialRing::prime_field(MODULUS, variables.iter().copied()).expect("the modulus is prime")
}

fn rational_ring(variables: &[&str]) -> PolynomialRing<Rationals> {
    PolynomialRing::rationals(variables.iter().copied()).expect("the names are variable names")
}

fn parse_all<D: sylvester::Domain>(ring: &PolynomialRing<D>, texts: &[&str]) -> Vec<Polynomial<D>> {
    texts
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the syntax holds"))
        .collect()
}

/// The basis of cyclic-3 over `F_p`, from the F4 engine.
fn cyclic_3_prime_basis(ring: &PolynomialRing) -> GroebnerBasis {
    ring.ideal(parse_all(ring, &CYCLIC_3))
        .expect("the generators share the ring")
        .groebner_basis(ComputeOptions::new())
        .expect("cyclic-3 fits every budget")
}

/// The basis of cyclic-3 over `Q`, through the checked constructor.
fn cyclic_3_rational_basis(ring: &PolynomialRing<Rationals>) -> GroebnerBasis<Rationals> {
    GroebnerBasis::<Rationals>::from_polynomials(
        ring,
        parse_all(ring, &CYCLIC_3_BASIS),
        Budget::new(),
    )
    .expect("the list is the reduced basis of cyclic-3")
}

/// Multiply two polynomials over `F_p`.
///
/// The crate offers no public multiplication, and the membership tests
/// need one to build a combination of the basis.
fn prime_product(ring: &PolynomialRing, f: &Polynomial, g: &Polynomial) -> Polynomial {
    let mut terms = Vec::new();
    for (left, left_exps) in f.terms() {
        for (right, right_exps) in g.terms() {
            let coeff = (left.value() * right.value()) % MODULUS;
            let exps: Vec<u16> = left_exps
                .iter()
                .zip(right_exps)
                .map(|(a, b)| a + b)
                .collect();
            terms.push((coeff, exps));
        }
    }
    ring.polynomial(terms).expect("the widths match the ring")
}

/// Multiply two polynomials over `Q`.
fn rational_product(
    ring: &PolynomialRing<Rationals>,
    f: &Polynomial<Rationals>,
    g: &Polynomial<Rationals>,
) -> Polynomial<Rationals> {
    let mut terms = Vec::new();
    for (left, left_exps) in f.terms() {
        for (right, right_exps) in g.terms() {
            let coeff: BigRational = left * right;
            let exps: Vec<u16> = left_exps
                .iter()
                .zip(right_exps)
                .map(|(a, b)| a + b)
                .collect();
            terms.push((coeff, exps));
        }
    }
    ring.polynomial(terms).expect("the widths match the ring")
}

/// Report whether no leading monomial of the basis divides a monomial of
/// `f`.
fn is_reduced<D: sylvester::Domain>(basis: &GroebnerBasis<D>, f: &Polynomial<D>) -> bool {
    let leads: Vec<&[u16]> = basis
        .iter()
        .map(|g| g.leading_term().expect("a basis element is nonzero").1)
        .collect();
    f.terms().all(|(_, exps)| {
        !leads
            .iter()
            .any(|lead| lead.iter().zip(exps).all(|(a, b)| a <= b))
    })
}

#[test]
fn every_generator_reduces_to_zero() {
    let ring = prime_ring(&["x", "y", "z"]);
    let basis = cyclic_3_prime_basis(&ring);
    for generator in parse_all(&ring, &CYCLIC_3) {
        assert!(
            basis
                .normal_form(&generator, Budget::new())
                .expect("the division fits every budget")
                .is_zero()
        );
        assert!(
            basis
                .contains(&generator, Budget::new())
                .expect("the division fits every budget")
        );
    }

    let rationals = rational_ring(&["x", "y", "z"]);
    let rational_basis = cyclic_3_rational_basis(&rationals);
    for generator in parse_all(&rationals, &CYCLIC_3) {
        assert!(
            rational_basis
                .contains(&generator, Budget::new())
                .expect("the division fits every budget")
        );
    }
}

#[test]
fn every_basis_element_reduces_to_zero() {
    let ring = prime_ring(&["x", "y", "z"]);
    let basis = cyclic_3_prime_basis(&ring);
    for element in basis.iter() {
        assert!(
            basis
                .normal_form(element, Budget::new())
                .expect("the division fits every budget")
                .is_zero()
        );
    }

    let rationals = rational_ring(&["x", "y", "z"]);
    let rational_basis = cyclic_3_rational_basis(&rationals);
    for element in rational_basis.iter() {
        assert!(
            rational_basis
                .normal_form(element, Budget::new())
                .expect("the division fits every budget")
                .is_zero()
        );
    }
}

#[test]
fn a_remainder_is_reduced_and_the_division_is_idempotent() {
    let ring = prime_ring(&["x", "y", "z"]);
    let basis = cyclic_3_prime_basis(&ring);
    let f = ring
        .parse_polynomial("x^3*y^2 + 7*x*y*z^4 - z^5 + 11")
        .expect("the syntax holds");

    let remainder = basis
        .normal_form(&f, Budget::new())
        .expect("the division fits every budget");
    assert!(!remainder.is_zero());
    assert!(is_reduced(&basis, &remainder));

    let again = basis
        .normal_form(&remainder, Budget::new())
        .expect("the division fits every budget");
    assert_eq!(again, remainder);
}

/// The reducer choice is the first divisor in basis order, so one input
/// has one remainder. Two calls, and two computations of the basis, agree.
#[test]
fn the_division_is_deterministic() {
    let ring = prime_ring(&["x", "y", "z"]);
    let f = ring
        .parse_polynomial("x^4 - 5*y^3*z + 2*z^2 + 9")
        .expect("the syntax holds");

    let first = cyclic_3_prime_basis(&ring);
    let second = cyclic_3_prime_basis(&ring);
    assert_eq!(first, second);

    let expected = first
        .normal_form(&f, Budget::new())
        .expect("the division fits every budget");
    for _ in 0..4 {
        let repeat = second
            .normal_form(&f, Budget::new())
            .expect("the division fits every budget");
        assert_eq!(repeat, expected);
    }
}

/// A polynomial outside the ideal keeps a nonzero remainder. Over
/// cyclic-3 the variable `x` reduces to `-y - z`, which no leading
/// monomial of the basis divides.
#[test]
fn a_polynomial_outside_the_ideal_has_a_nonzero_normal_form() {
    let ring = prime_ring(&["x", "y", "z"]);
    let basis = cyclic_3_prime_basis(&ring);
    let x = ring.parse_polynomial("x").expect("the syntax holds");
    let remainder = basis
        .normal_form(&x, Budget::new())
        .expect("the division fits every budget");
    assert_eq!(
        remainder,
        ring.parse_polynomial("-y - z").expect("the syntax holds")
    );
    assert!(
        !basis
            .contains(&x, Budget::new())
            .expect("the division fits every budget")
    );

    let rationals = rational_ring(&["x", "y", "z"]);
    let rational_basis = cyclic_3_rational_basis(&rationals);
    let x = rationals.parse_polynomial("x").expect("the syntax holds");
    assert_eq!(
        rational_basis
            .normal_form(&x, Budget::new())
            .expect("the division fits every budget"),
        rationals
            .parse_polynomial("-y - z")
            .expect("the syntax holds")
    );
    assert!(
        !rational_basis
            .contains(&x, Budget::new())
            .expect("the division fits every budget")
    );
}

/// Every combination of the basis elements is a member, whatever the
/// cofactors are.
#[test]
fn a_combination_over_the_rationals_is_a_member() {
    let ring = rational_ring(&["x", "y", "z"]);
    let basis = cyclic_3_rational_basis(&ring);
    let cofactors = parse_all(&ring, &["1/2*x^2 - 3", "y*z + 1/3", "-2/5*z^3 + x"]);

    let mut sum = ring.zero();
    for (cofactor, element) in cofactors.iter().zip(basis.iter()) {
        let product = rational_product(&ring, cofactor, element);
        let mut terms: Vec<(BigRational, Vec<u16>)> = Vec::new();
        for (coeff, exps) in sum.terms().chain(product.terms()) {
            terms.push((coeff.clone(), exps.to_vec()));
        }
        sum = ring.polynomial(terms).expect("the widths match the ring");
    }

    assert!(!sum.is_zero());
    assert!(
        basis
            .contains(&sum, Budget::new())
            .expect("the division fits every budget")
    );

    let stray = rational_product(
        &ring,
        &ring.parse_polynomial("1/7*x").expect("the syntax holds"),
        &ring.parse_polynomial("y + 1").expect("the syntax holds"),
    );
    assert!(
        !basis
            .contains(&stray, Budget::new())
            .expect("the division fits every budget")
    );
}

proptest! {
    /// A combination of the basis elements over `F_p` is a member. The
    /// cofactors are arbitrary, so this is the membership direction that
    /// needs no oracle.
    #[test]
    fn a_combination_over_a_prime_field_is_a_member(
        cofactors in proptest::collection::vec(
            proptest::collection::vec(
                (1u64..MODULUS, proptest::collection::vec(0u16..4, 3..=3)),
                0..4,
            ),
            3..=3,
        ),
    ) {
        let ring = prime_ring(&["x", "y", "z"]);
        let basis = cyclic_3_prime_basis(&ring);

        let mut terms: Vec<(u64, Vec<u16>)> = Vec::new();
        for (cofactor, element) in cofactors.iter().zip(basis.iter()) {
            let cofactor = ring
                .polynomial(cofactor.clone())
                .expect("the widths match the ring");
            let product = prime_product(&ring, &cofactor, element);
            for (coeff, exps) in product.terms() {
                terms.push((coeff.value(), exps.to_vec()));
            }
        }
        let combination = ring.polynomial(terms).expect("the widths match the ring");

        prop_assert!(basis
            .contains(&combination, Budget::new())
            .expect("the division fits every budget"));
    }
}

/// Dividing `x^65535*y` by `x^65535 + y^65535` needs `y^65536`. The
/// leading monomials and the quotient monomial all fit; the tail multiple
/// does not.
#[test]
fn an_exponent_past_the_width_is_an_error_and_not_a_panic() {
    let ring = prime_ring(&["x", "y"]);
    let divisor = ring
        .polynomial([(1u64, [65535u16, 0]), (1u64, [0, 65535])])
        .expect("the widths match the ring");
    let basis = GroebnerBasis::<PrimeField>::from_polynomials(&ring, vec![divisor], Budget::new())
        .expect("one monic element is a reduced basis");
    let f = ring
        .polynomial([(1u64, [65535u16, 1])])
        .expect("the widths match the ring");

    assert_eq!(
        basis.normal_form(&f, Budget::new()),
        Err(NormalFormError::ExponentLimit { limit: 65535 })
    );
}

#[test]
fn an_exhausted_budget_is_typed() {
    let ring = prime_ring(&["x", "y", "z"]);
    let basis = cyclic_3_prime_basis(&ring);
    let f = ring
        .parse_polynomial("x^3*y^2 + z^4")
        .expect("the syntax holds");

    assert_eq!(
        basis.normal_form(&f, Budget::new().timeout(Duration::ZERO)),
        Err(NormalFormError::Timeout)
    );
    assert_eq!(
        basis.normal_form(&f, Budget::new().memory_limit(0)),
        Err(NormalFormError::MemoryLimitExceeded)
    );

    let polynomials = basis.clone().into_polynomials();
    assert_eq!(
        GroebnerBasis::<PrimeField>::from_polynomials(
            &ring,
            polynomials,
            Budget::new().timeout(Duration::ZERO)
        ),
        Err(BasisError::Timeout)
    );
}

/// The bytes one term costs the budget of a division.
///
/// The meter is `Polynomial::heap_bytes`, which is not public, so the
/// figure is the smallest memory limit under which a division that holds
/// one term runs. `f` must hold one term no basis element divides.
fn one_term_bytes(basis: &GroebnerBasis, f: &Polynomial) -> usize {
    let fits = |bytes: usize| {
        basis
            .normal_form(f, Budget::new().memory_limit(bytes))
            .is_ok()
    };
    let mut low = 0usize;
    let mut high = 4096usize;
    assert!(fits(high), "the search starts above what one term needs");
    while low + 1 < high {
        let middle = (low + high) / 2;
        if fits(middle) {
            high = middle;
        } else {
            low = middle;
        }
    }
    high
}

/// A reduction step charges the polynomial it is about to build.
///
/// The step holds the working value and allocates a replacement of at
/// most the terms of both operands, so it charges three times what the
/// working value holds. A budget that covers the working value and not
/// the replacement reports the limit instead of allocating.
#[test]
fn a_step_charges_the_replacement_before_it_allocates_it() {
    const TAIL: u16 = 200;
    let ring = prime_ring(&["x", "y"]);
    let mut terms = vec![(1i64, [TAIL, 0])];
    for exponent in 0..TAIL {
        terms.push((1, [0, exponent]));
    }
    // Every tail term has a degree below TAIL, so the leading monomial is
    // x^TAIL and one monic element is the reduced basis of its own ideal.
    let g = ring
        .polynomial(terms)
        .expect("the exponent vectors match the ring");
    let basis =
        GroebnerBasis::<PrimeField>::from_polynomials(&ring, vec![g.clone()], Budget::new())
            .expect("one monic polynomial is a reduced basis");
    let outside = ring.parse_polynomial("y^300").expect("the syntax holds");
    let one_term = one_term_bytes(&basis, &outside);
    let held = usize::from(TAIL + 1) * one_term;

    // Twice what g holds covers g and not the replacement, which is
    // three times it.
    assert_eq!(
        basis.normal_form(&g, Budget::new().memory_limit(2 * held)),
        Err(NormalFormError::MemoryLimitExceeded)
    );
    assert!(
        basis
            .normal_form(&g, Budget::new().memory_limit(4 * held))
            .is_ok_and(|remainder| remainder.is_zero())
    );
}

#[test]
fn a_polynomial_of_another_ring_is_a_ring_mismatch() {
    let ring = prime_ring(&["x", "y", "z"]);
    let basis = cyclic_3_prime_basis(&ring);
    let other = PolynomialRing::prime_field(7, ["x", "y", "z"]).expect("7 is prime");
    let f = other.parse_polynomial("x + 1").expect("the syntax holds");

    assert_eq!(
        basis.normal_form(&f, Budget::new()),
        Err(NormalFormError::RingMismatch)
    );
    assert_eq!(
        GroebnerBasis::<PrimeField>::from_polynomials(&ring, vec![f], Budget::new()),
        Err(BasisError::RingMismatch)
    );
}

/// The constructor accepts what the engine computed, and the value it
/// builds is the value the engine returned.
#[test]
fn the_constructor_accepts_a_computed_basis() {
    let ring = prime_ring(&["x", "y", "z"]);
    let computed = cyclic_3_prime_basis(&ring);
    let checked = GroebnerBasis::<PrimeField>::from_polynomials(
        &ring,
        computed.clone().into_polynomials(),
        Budget::new(),
    )
    .expect("the engine returns a reduced basis");
    assert_eq!(checked, computed);

    let rationals = rational_ring(&["x", "y", "z"]);
    let basis = cyclic_3_rational_basis(&rationals);
    assert!(basis.lift().is_none(), "no modular run produced it");
}

#[test]
fn the_constructor_names_the_property_the_list_does_not_have() {
    let ring = prime_ring(&["x", "y", "z"]);
    let budget = Budget::new();

    assert_eq!(
        GroebnerBasis::<PrimeField>::from_polynomials(&ring, vec![ring.zero()], budget),
        Err(BasisError::ZeroPolynomial { index: 0 })
    );

    let not_monic = parse_all(&ring, &["2*x + y"]);
    assert_eq!(
        GroebnerBasis::<PrimeField>::from_polynomials(&ring, not_monic, budget),
        Err(BasisError::NotMonic { index: 0 })
    );

    let mut unsorted = parse_all(&ring, &CYCLIC_3_BASIS);
    unsorted.reverse();
    assert_eq!(
        GroebnerBasis::<PrimeField>::from_polynomials(&ring, unsorted, budget),
        Err(BasisError::NotSorted { index: 1 })
    );

    let twice = parse_all(&ring, &["x + y + z", "x + y + z"]);
    assert_eq!(
        GroebnerBasis::<PrimeField>::from_polynomials(&ring, twice, budget),
        Err(BasisError::NotSorted { index: 1 })
    );

    let plane = prime_ring(&["x", "y"]);
    let not_interreduced = parse_all(&plane, &["x^2", "x"]);
    assert_eq!(
        GroebnerBasis::<PrimeField>::from_polynomials(&plane, not_interreduced, budget),
        Err(BasisError::NotInterreduced { index: 0 })
    );

    let not_groebner = parse_all(&plane, &["x^2 - y", "x*y - 1"]);
    assert_eq!(
        GroebnerBasis::<PrimeField>::from_polynomials(&plane, not_groebner, budget),
        Err(BasisError::NotGroebner { left: 0, right: 1 })
    );
}

/// Over `Q` the same checks run on exact coefficients.
#[test]
fn the_constructor_checks_the_rational_domain_too() {
    let ring = rational_ring(&["x", "y"]);
    let budget = Budget::new();

    let not_monic = parse_all(&ring, &["1/2*x + y"]);
    assert_eq!(
        GroebnerBasis::<Rationals>::from_polynomials(&ring, not_monic, budget),
        Err(BasisError::NotMonic { index: 0 })
    );

    let not_groebner = parse_all(&ring, &["x^2 - y", "x*y - 1/3"]);
    assert_eq!(
        GroebnerBasis::<Rationals>::from_polynomials(&ring, not_groebner, budget),
        Err(BasisError::NotGroebner { left: 0, right: 1 })
    );

    let not_interreduced = parse_all(&ring, &["x^2 + y", "x"]);
    assert_eq!(
        GroebnerBasis::<Rationals>::from_polynomials(&ring, not_interreduced, budget),
        Err(BasisError::NotInterreduced { index: 0 })
    );
}
