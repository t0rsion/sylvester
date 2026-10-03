//! Measurements for the v2 verifier on fixed certificate bytes.
//!
//! Run this ignored test in release mode with the benchmark cores pinned.
//! It reports the median of repeated checks, so an arithmetic change can be
//! compared on the same generated certificate.

use std::time::Instant;

use sylvester::verify::{Limits, verify_with_limits};
use sylvester::{ComputeOptions, PolynomialRing};

const P: u64 = 1_073_741_827;
const RUNS: usize = 11;

fn cyclic_six_certificate() -> Vec<u8> {
    let names: Vec<String> = (0..6).map(|index| format!("x{index}")).collect();
    let ring = PolynomialRing::prime_field(P, names).expect("the modulus is prime");
    let texts = [
        "x0 + x1 + x2 + x3 + x4 + x5",
        "x0*x1 + x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x0",
        "x0*x1*x2 + x1*x2*x3 + x2*x3*x4 + x3*x4*x5 + x4*x5*x0 + x5*x0*x1",
        "x0*x1*x2*x3 + x1*x2*x3*x4 + x2*x3*x4*x5 + x3*x4*x5*x0 + x4*x5*x0*x1 + x5*x0*x1*x2",
        "x0*x1*x2*x3*x4 + x1*x2*x3*x4*x5 + x2*x3*x4*x5*x0 + x3*x4*x5*x0*x1 + x4*x5*x0*x1*x2 + x5*x0*x1*x2*x3",
        "x0*x1*x2*x3*x4*x5 - 1",
    ];
    let generators: Vec<_> = texts
        .into_iter()
        .map(|text| ring.parse_polynomial(text).expect("the polynomial parses"))
        .collect();
    ring.ideal(generators)
        .expect("the generators share a ring")
        .groebner_basis_certified(ComputeOptions::new().threads(1))
        .expect("the certificate holds")
        .certificate()
        .to_vec()
}

fn median_check(bytes: &[u8]) -> f64 {
    let limits = Limits {
        max_work_units: u64::MAX,
        ..Limits::default()
    };
    let mut samples = Vec::with_capacity(RUNS);
    for _ in 0..RUNS {
        let start = Instant::now();
        verify_with_limits(bytes, &limits).expect("the certificate holds");
        samples.push(start.elapsed().as_secs_f64());
    }
    samples.sort_by(f64::total_cmp);
    samples[RUNS / 2]
}

#[test]
#[ignore = "a measurement, not a check"]
fn v2_verifier_median_on_fixed_certificates() {
    let generated = cyclic_six_certificate();
    let cases = [
        (
            "tiny",
            include_bytes!("fixtures/cert-v2/tiny.cert").as_slice(),
        ),
        (
            "square",
            include_bytes!("fixtures/cert-v2/square.cert").as_slice(),
        ),
        ("cyclic-6", generated.as_slice()),
    ];
    for (name, bytes) in cases {
        println!(
            "{name}: {} bytes, median {:.6} s",
            bytes.len(),
            median_check(bytes)
        );
    }
}
