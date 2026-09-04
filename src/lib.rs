//! Gröbner bases over prime fields and the rationals.
//!
//! The crate provides an F4 engine, a multimodular rational driver, a
//! classic F5 oracle, and independent certificate verifiers.
//!
//! Every value is built through a [`PolynomialRing`], which fixes the
//! coefficient domain, the variable names, and the variable order. The
//! domain is a sealed type parameter. [`PrimeField`] is the default, so
//! `PolynomialRing` alone names a ring over `F_p`, and [`Rationals`] is
//! the exact domain `Q`. The monomial order is grevlex over that variable
//! order, and certificates name it `grevlex-v1`. [`Ideal::groebner_basis`]
//! returns the basis the engine computed.
//! [`Ideal::groebner_basis_certified`] returns one a verifier accepted.
//!
//! Over `Q` the engine is multimodular: it computes a basis modulo many
//! primes and lifts the residues to rational numbers. That path is a
//! heuristic and has no certified form. `GroebnerBasis::lift` carries the
//! counters of the run and names what its stopping rule observed.
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
//! # Correctness status: checked, not proven
//!
//! Two counterexamples proved the two engines of the original extraction
//! incorrect. This crate carries the repair. `KNOWN_ISSUES.md` records the
//! systems, the wrong outputs, the witnesses, and the root causes.
//!
//! The counterexample tests in `tests/known_defects.rs` assert correct
//! behavior through an independent checker and pass here. So does a
//! randomized differential suite against a self-contained Buchberger
//! oracle, in `tests/differential.rs`. The F4 backend is checked three
//! ways: against the classic backend, against that Buchberger oracle, and
//! against msolve on the comparison record under
//! `benchmarks/gb-comparison`. Termination of the classic backend has a
//! pen-and-paper proof by Dickson's lemma, also in `KNOWN_ISSUES.md`. No
//! proof here is machine-checked, and neither backend is proven correct.
//!
//! [`Ideal::groebner_basis`] is the unproven path. It reports an exhausted
//! budget, or a monomial past the width of one exponent, and nothing else
//! about the basis it returns.
//! [`Ideal::groebner_basis_certified`] runs the backend the options name
//! and writes a certificate for it. It returns a value only after the
//! independent verifier in [`verify`] accepts the bytes, and the basis it
//! returns is decoded from those bytes. The format follows the backend:
//! the classic backend writes `sylv-gb-cert-v1` from the cofactors it
//! tracks, and F4 writes `sylv-gb-cert-v2` from the trace it records.
//!
//! # Backends
//!
//! [`Backend::F4`] is the default. It batches critical pairs by degree,
//! builds one sparse matrix per batch over interned monomials, and reduces
//! it. Pair management is Gebauer-Moller. There are no signatures in it,
//! and it writes `sylv-gb-cert-v2`.
//!
//! [`Backend::Classic`] processes critical pairs one at a time in
//! signature order, as classic F5 does. It is much slower than F4 on the
//! comparison record. It stays as the oracle the differential suite
//! compares against and as the `sylv-gb-cert-v1` path.
//!
//! # Threads
//!
//! [`ComputeOptions::threads`] sets the thread count of one computation.
//! The default is the rayon global pool, which reads `RAYON_NUM_THREADS`.
//! The F4 kernel reduces the rows of one batch against the frozen pivots
//! in parallel once the batch passes an internal work threshold. Below the
//! threshold, and at a count of 1, every row reduces on the calling
//! thread. The classic backend computes on one thread. The thread count
//! changes no byte of a basis or of a certificate. The crate has no
//! feature flags.

#![deny(unsafe_code)]
#![warn(missing_docs)]

mod cert;
mod certificate;
mod compute;
mod hilbert;
mod ideal;
mod normal_form;
mod poly;
mod ring;
pub mod verify;

pub use certificate::{
    CertificateCap, CertifiedGroebnerBasis, CertifyError, EmitterFault, Place, TraceFault,
};
pub use compute::{
    Backend, Budget, ComputeError, ComputeOptions, ComputeReport, F4Counters, RationalOptions,
    RationalStop,
};
pub use hilbert::{HilbertError, HilbertSeries};
pub use ideal::{GroebnerBasis, Ideal};
pub use normal_form::{BasisError, NormalFormError};
pub use poly::Polynomial;
pub use ring::{
    Coefficient, Domain, DomainOps, Established, Felt, ModularLift, ParseError, PolynomialRing,
    PrimeField, PrimeOps, RationalMeta, RationalOps, Rationals, RingError,
};
