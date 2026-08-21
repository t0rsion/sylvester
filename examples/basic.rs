//! Compute the Gröbner basis of cyclic-3 over F_32003.

use sylvester::{ComputeOptions, PolynomialRing};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"])?;
    let ideal = ring.ideal([
        ring.parse_polynomial("x + y + z")?,
        ring.parse_polynomial("x*y + y*z + z*x")?,
        ring.parse_polynomial("x*y*z - 1")?,
    ])?;

    let basis = ideal.groebner_basis(ComputeOptions::new())?;
    println!("the basis holds {} polynomials", basis.len());
    for polynomial in &basis {
        println!("  {polynomial}");
    }
    Ok(())
}
