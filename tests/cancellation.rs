//! Cancellation stops prime-field and multimodular computations.

use std::sync::atomic::Ordering;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

use sylvester::{
    Backend, Budget, CancellationToken, ComputeError, ComputeOptions, Ideal, Polynomial,
    PolynomialRing, RationalOptions, Rationals,
};

fn prime_ideal() -> Ideal {
    let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"])
        .expect("the modulus and variables are valid");
    ring.ideal([
        ring.parse_polynomial("x^2 + y^2 + z^2 - 1")
            .expect("the polynomial parses"),
        ring.parse_polynomial("x*y - z")
            .expect("the polynomial parses"),
        ring.parse_polynomial("y*z - x")
            .expect("the polynomial parses"),
    ])
    .expect("the generators belong to one ring")
}

fn cyclic(ring: &PolynomialRing<Rationals>, n: usize) -> Vec<Polynomial<Rationals>> {
    let mut system = Vec::with_capacity(n);
    for degree in 1..n {
        let mut terms = Vec::with_capacity(n);
        for offset in 0..n {
            let mut exponents = vec![0u16; n];
            for step in 0..degree {
                exponents[(offset + step) % n] += 1;
            }
            terms.push((1, exponents));
        }
        system.push(ring.polynomial(terms).expect("the exponents fit"));
    }
    system.push(
        ring.polynomial(vec![(1, vec![1u16; n]), (-1, vec![0u16; n])])
            .expect("the exponents fit"),
    );
    system
}

#[test]
fn an_already_cancelled_token_stops_both_prime_field_backends() {
    for backend in [Backend::F4, Backend::Classic] {
        let token = CancellationToken::new();
        token.cancel();
        let options = ComputeOptions::new()
            .backend(backend)
            .budget(Budget::new().cancellation(token.clone()));
        assert_eq!(
            prime_ideal().groebner_basis(options),
            Err(ComputeError::Timeout)
        );
        assert!(token.is_cancelled());
    }
}

#[test]
fn an_already_cancelled_token_stops_the_multimodular_driver() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the variables are valid");
    let generator = ring
        .parse_polynomial("x^2 - y")
        .expect("the polynomial parses");
    let ideal = ring
        .ideal([generator])
        .expect("the generator belongs to the ring");
    let token = CancellationToken::new();
    token.cancel();
    let options = RationalOptions::new().compute(
        ComputeOptions::new()
            .threads(8)
            .budget(Budget::new().cancellation(token.clone())),
    );
    assert_eq!(ideal.groebner_basis(options), Err(ComputeError::Timeout));
    assert!(token.is_cancelled());
}

#[test]
fn cancellation_during_a_modular_wave_joins_every_worker() {
    let ring = PolynomialRing::rationals(["x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7"])
        .expect("the variables are valid");
    let ideal = ring
        .ideal(cyclic(&ring, 8))
        .expect("the generators belong to one ring");
    let token = CancellationToken::new();
    let cancel = token.clone();
    let (started_tx, started_rx) = mpsc::sync_channel(1);
    let worker = thread::spawn(move || {
        started_tx.send(()).expect("the caller is waiting");
        ideal.groebner_basis(
            RationalOptions::new().compute(
                ComputeOptions::new()
                    .threads(8)
                    .budget(Budget::new().cancellation(token)),
            ),
        )
    });
    started_rx.recv().expect("the computation thread started");
    thread::sleep(Duration::from_millis(20));
    cancel.cancel();
    assert_eq!(
        worker.join().expect("all modular workers joined"),
        Err(ComputeError::Timeout)
    );
}

#[test]
fn token_equality_uses_identity_after_cancellation() {
    let first = CancellationToken::new();
    let clone = first.clone();
    let shared_flag = first.shared_flag();
    let second = CancellationToken::new();
    assert_eq!(first, clone);
    assert_ne!(first, second);
    first.cancel();
    second.cancel();
    assert!(shared_flag.load(Ordering::Relaxed));
    assert_eq!(first, clone);
    assert_ne!(first, second);
}
