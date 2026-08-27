//! The public half of the monomial tests of design section 8.2.
//!
//! `MonomialTable` and every type around it are `pub(crate)`, so an
//! integration test cannot reach them. The reference proptest against
//! `poly::Monomial`, the packing edges, the mask soundness sweep, the
//! cross-store cases, and the overflow cases are unit tests in
//! `src/compute/f4/monomial.rs`. What stays here is the one monomial fact
//! the public API states: the exponent contract the table holds to at every
//! lane width.

use sylvester::{ParseError, PolynomialRing, RingError};

#[test]
fn the_public_exponent_contract_stops_at_65535() {
    let ring = PolynomialRing::prime_field(32003, ["x", "y"]).unwrap();
    assert!(ring.parse_polynomial("x^65535").is_ok());
    assert!(matches!(
        ring.parse_polynomial("x^65536"),
        Err(RingError::Parse(ParseError::ExponentTooLarge { .. }))
    ));
    assert!(ring.polynomial([(1i64, [65535u16, 0])]).is_ok());
}
