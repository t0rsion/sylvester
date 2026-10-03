//! A certified input ideal and computations in its finite quotient.

use sylvester::{Budget, ComputeOptions, FiniteQuotient, PolynomialRing};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let quotient = nilpotent_quotient()?;
    let x = quotient.ring().generator(0).expect("the ring has x");
    assert_eq!(quotient.vector_space_dimension(), 4);
    let square = quotient.pow(&x, 2, Budget::new())?;
    assert!(square.is_zero());
    let minimal = quotient.minimal_polynomial(&x, Budget::new())?;
    let characteristic = quotient.characteristic_polynomial(&x, Budget::new())?;
    assert_eq!(minimal.to_string(), "t^2");
    assert_eq!(characteristic.to_string(), "t^4");
    println!("dimension: {}", quotient.vector_space_dimension());
    println!("minimal polynomial of x: {minimal}");
    println!("characteristic polynomial of multiplication by x: {characteristic}");
    Ok(())
}

fn nilpotent_quotient() -> Result<FiniteQuotient, Box<dyn std::error::Error>> {
    let ring = PolynomialRing::prime_field(7, ["x", "y"])?;
    let input = ring.ideal([ring.parse_polynomial("x^2")?, ring.parse_polynomial("y^2")?])?;
    let certified = input.groebner_basis_certified(ComputeOptions::new())?;
    Ok(certified.basis().finite_quotient(Budget::new())?)
}
