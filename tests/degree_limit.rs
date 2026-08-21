//! Inputs past the degree an exponent's width supports.
//!
//! Every run here must end in a typed error. A panic is the defect these
//! tests exist to block: the engines used to multiply monomials with an
//! `expect` and abort on the product.

use sylvester::{Backend, ComputeError, ComputeOptions, PolynomialRing};

const LIMIT: ComputeError = ComputeError::DegreeLimit { limit: 65535 };

fn backends() -> [Backend; 2] {
    [Backend::Classic, Backend::Matrix]
}

/// A generator of total degree 131070, which no product of the engines
/// could hold. The ring accepts it because each exponent fits a `u16` on
/// its own.
#[test]
fn an_input_generator_past_the_limit_is_a_typed_error() {
    let ring = PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime");
    let f = ring
        .polynomial([(1, [65535, 0]), (1, [0, 65535])])
        .expect("each exponent fits");
    let g = ring
        .polynomial([(1, [65535, 65535])])
        .expect("each exponent fits");
    assert_eq!(g.degree(), Some(131070));

    let ideal = ring.ideal([f, g]).expect("one ring");
    for backend in backends() {
        assert_eq!(
            ideal.groebner_basis(ComputeOptions::new().backend(backend)),
            Err(LIMIT),
            "backend {backend:?}"
        );
    }
    assert_eq!(
        ideal.groebner_basis_certified(ComputeOptions::new()),
        Err(sylvester::CertifyError::Engine(LIMIT))
    );
}

/// Every generator is inside the limit and every critical pair's least
/// common multiple is too, but a signature product leaves the width one
/// exponent holds.
#[test]
fn a_signature_product_past_the_limit_is_a_typed_error() {
    let ring = PolynomialRing::prime_field(2, ["x"]).expect("2 is prime");
    let f0 = ring
        .polynomial([(1, [65535]), (1, [40000]), (1, [0])])
        .expect("each exponent fits");
    let f1 = ring
        .polynomial([(1, [21845]), (1, [3]), (1, [1])])
        .expect("each exponent fits");
    let f2 = ring.polynomial([(1, [21845])]).expect("each exponent fits");
    for f in [&f0, &f1, &f2] {
        assert!(f.degree().is_some_and(|deg| deg <= 65535));
    }

    let ideal = ring.ideal([f0, f1, f2]).expect("one ring");
    for backend in backends() {
        assert_eq!(
            ideal.groebner_basis(ComputeOptions::new().backend(backend)),
            Err(LIMIT),
            "backend {backend:?}"
        );
    }
    assert_eq!(
        ideal.groebner_basis_certified(ComputeOptions::new()),
        Err(sylvester::CertifyError::Engine(LIMIT))
    );
}
