//! Signature-based Gröbner bases over finite prime fields.
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
//! Two counterexamples proved both backends incorrect as originally
//! extracted, and this crate carries the repair. The counterexample tests
//! in `tests/known_defects.rs` assert correct behavior through an
//! independent checker, and a randomized differential suite checks both
//! backends against a self-contained Buchberger oracle in
//! `tests/differential.rs`. Termination has a pen-and-paper proof by
//! Dickson's lemma. No proof here is machine-checked, and the verifier
//! itself is a trusted, unproven base.
//!
//! [`Ideal::groebner_basis`] is the unproven path. It reports an exhausted
//! budget, or a monomial past the degree an exponent's width supports (an
//! input generator, a critical pair's least common multiple, or a signature
//! or cofactor product), and nothing else about the basis it returns.
//! [`Ideal::groebner_basis_certified`] runs the classic backend with
//! cofactor tracking and writes a certificate. It returns a value only
//! after the independent verifier in [`verify`] accepts the bytes, and the
//! basis it returns is decoded from those bytes. Certification is
//! classic-only in this release.
//!
//! # Backends
//!
//! [`Backend::Classic`] processes critical pairs one at a time in signature
//! order. [`Backend::Matrix`] batches pairs by degree and reduces them in
//! sparse Macaulay matrices with F4-style elimination. It applies the F5
//! syzygy criterion to skip rows, and it replaces a rewritable
//! S-polynomial row with the multiple its canonical rewriter names.
//!
//! # Feature flags
//!
//! `parallel` adds rayon row reduction to the matrix backend. Rows that
//! share one signature reduce against the same pivots, so such a group
//! runs in parallel once it reaches an internal threshold. There is no
//! option for it.

#![deny(unsafe_code)]
#![warn(missing_docs)]

mod cert;
mod certificate;
mod compute;
mod ideal;
mod poly;
mod ring;
pub mod verify;

pub use certificate::{CertifiedGroebnerBasis, CertifyError, EmitterFault, Place};
pub use compute::{Backend, ComputeError, ComputeOptions};
pub use ideal::{GroebnerBasis, Ideal};
pub use poly::Polynomial;
pub use ring::{ParseError, PolynomialRing, RingError};
