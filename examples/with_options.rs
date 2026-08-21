//! Compute a certified basis under a deadline and a memory limit.

use std::time::Duration;
use sylvester::{Backend, ComputeOptions, PolynomialRing};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ring = PolynomialRing::prime_field(32003, ["x", "y"])?;
    let ideal = ring.ideal([
        ring.parse_polynomial("x^2 - 1")?,
        ring.parse_polynomial("x*y - 1")?,
    ])?;

    let options = ComputeOptions::new()
        .backend(Backend::Classic)
        .timeout(Duration::from_secs(2))
        .memory_limit(512 << 20);

    // `groebner_basis_certified` always runs the classic backend and does not
    // read `backend`; the option here picks the backend for `groebner_basis`
    // alone.
    let basis = ideal.groebner_basis(options.clone())?;
    println!("the basis holds {} polynomials", basis.len());

    let certified = ideal.groebner_basis_certified(options)?;
    println!(
        "the verifier accepted {} certificate bytes",
        certified.certificate().len()
    );
    Ok(())
}
