//! Gröbner bases over finite prime fields, with a default F4 engine, a
//! classic F5 oracle, and independent certificate verifiers.
//!
//! Every value is built through a [`PolynomialRing`], which fixes the prime,
//! the variable names, and the variable order. The monomial order is
//! grevlex over that variable order, and certificates name it
//! `grevlex-v1`.
//!
//! ```
//! use sylvester::{ComputeOptions, PolynomialRing};
//!
//! let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"])?;
//! let ideal = ring.ideal([
//!     ring.parse_polynomial("x^2*y - 3*z + 1")?,
//!     ring.parse_polynomial("x*z - y")?,
//! ])?;
//!
//! let basis = ideal.groebner_basis(ComputeOptions::new())?;
//! for polynomial in &basis {
//!     println!("{polynomial}");
//! }
//!
//! let certified = ideal.groebner_basis_certified(ComputeOptions::new())?;
//! assert_eq!(certified.basis().len(), basis.len());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Correctness status: repaired, checked against oracles, not proven
//!
//! Counterexamples proved the originally extracted engines incorrect.
//! This crate carries the repair. The tests in `tests/known_defects.rs`
//! assert the hand-verified correct bases through an independent checker.
//! A randomized differential suite checks F4 and classic against a
//! self-contained Buchberger oracle in `tests/differential.rs`.
//! Termination has a pen-and-paper proof by Dickson's lemma. No proof
//! here is machine-checked. The verifier is trusted, unproven code.
//!
//! [`Ideal::groebner_basis`] is the unproven path. It reports an exhausted
//! budget, or a monomial past the degree an exponent's width supports: an
//! input generator, a critical pair's least common multiple, or a signature
//! or cofactor product. It says nothing else about the basis it returns.
//! [`Ideal::groebner_basis_certified`] runs the selected backend and writes
//! its certificate. It returns a value only after the independent verifier
//! in [`verify`] accepts the bytes, and the basis it returns is decoded
//! from those bytes. Classic writes `sylv-gb-cert-v1`; F4 writes
//! `sylv-gb-cert-v2`.
//!
//! # Backends
//!
//! [`Backend::F4`] is the default. It batches critical pairs by degree,
//! builds sparse matrices over interned monomials, and uses Gebauer-Moller
//! pair management. [`Backend::Classic`] processes pairs one at a time in
//! signature order. It is the differential oracle and the v1 path.
//!
//! [`ComputeOptions::threads`] controls F4's rayon parallelism. Classic and
//! certified runs stay on one thread. The crate has no feature flags.

#![deny(unsafe_code)]
#![warn(missing_docs)]
#![warn(unreachable_pub)]

mod cert;
mod certificate;
mod compute;
mod ideal;
mod poly;
mod ring;
pub mod verify;

pub use certificate::{
    CertificateCap, CertifiedGroebnerBasis, CertifyError, EmitterFault, Place, TraceFault,
};
pub use compute::{Backend, ComputeError, ComputeOptions, ComputeReport, F4Counters};
pub use ideal::{GroebnerBasis, Ideal};
pub use poly::Polynomial;
pub use ring::{ParseError, PolynomialRing, RingError};
