//! Ideals and the bases computed from them.

use std::ops::Deref;

use crate::certificate::{CertifiedGroebnerBasis, CertifyError};
use crate::compute::{self, ComputeError, ComputeOptions};
use crate::poly::Polynomial;
use crate::ring::PolynomialRing;

/// The ideal a list of polynomials generates.
///
/// Build one with [`PolynomialRing::ideal`].
///
/// ```
/// use sylvester::{ComputeOptions, PolynomialRing};
///
/// let ring = PolynomialRing::prime_field(32003, ["x", "y"])?;
/// let ideal = ring.ideal([
///     ring.parse_polynomial("x^2 - 1")?,
///     ring.parse_polynomial("x*y - 1")?,
/// ])?;
/// let basis = ideal.groebner_basis(ComputeOptions::new())?;
/// assert_eq!(basis.len(), 2);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ideal {
    ring: PolynomialRing,
    generators: Vec<Polynomial>,
}

impl Ideal {
    pub(crate) fn new(ring: PolynomialRing, generators: Vec<Polynomial>) -> Self {
        Ideal { ring, generators }
    }

    /// The ring the ideal lives in.
    pub fn ring(&self) -> &PolynomialRing {
        &self.ring
    }

    /// The generators, in the order the caller gave them.
    pub fn generators(&self) -> &[Polynomial] {
        &self.generators
    }

    /// Compute the reduced Gröbner basis under grevlex.
    ///
    /// The options pick the backend and the resource budget. Nothing else
    /// is checked. The result carries no proof: it is the value the engine
    /// returned, unchecked. Use [`Ideal::groebner_basis_certified`] for a
    /// basis an independent verifier has accepted.
    pub fn groebner_basis(&self, options: ComputeOptions) -> Result<GroebnerBasis, ComputeError> {
        let polynomials = compute::groebner_basis(&self.ring, &self.generators, &options)?;
        Ok(GroebnerBasis {
            ring: self.ring.clone(),
            polynomials,
        })
    }

    /// Compute the reduced Gröbner basis and the certificate that proves
    /// it.
    ///
    /// The value exists only after the independent verifier in
    /// [`crate::verify`] accepts the certificate bytes. The basis it
    /// carries is decoded from those bytes, not taken from the engine.
    ///
    /// This release certifies the classic backend alone. The backend in
    /// `options` is not read. The deadline and the memory limit are, and
    /// they cover the whole run: the engine, the certificate, and the
    /// verifier. An exhausted budget is [`CertifyError::Engine`] up to the
    /// verifier and [`CertifyError::VerifierExhausted`] inside it. Neither
    /// says anything about the basis.
    ///
    /// [`CertifyError::Rejected`] means the engine wrote a certificate that
    /// does not hold. That is an engine defect, reported rather than
    /// panicked.
    pub fn groebner_basis_certified(
        &self,
        options: ComputeOptions,
    ) -> Result<CertifiedGroebnerBasis, CertifyError> {
        compute::groebner_basis_certified(&self.ring, &self.generators, &options)
    }
}

/// The basis a computation returned.
///
/// The elements are monic and run strictly descending by leading monomial.
/// The value derefs to `&[Polynomial]`, so slice methods and iteration work
/// on it directly.
///
/// A basis from [`Ideal::groebner_basis`] is the engine's claim. A basis
/// from [`Ideal::groebner_basis_certified`] passed the verifier.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GroebnerBasis {
    ring: PolynomialRing,
    polynomials: Vec<Polynomial>,
}

impl GroebnerBasis {
    pub(crate) fn new(ring: PolynomialRing, polynomials: Vec<Polynomial>) -> Self {
        GroebnerBasis { ring, polynomials }
    }

    /// The ring the basis lives in.
    pub fn ring(&self) -> &PolynomialRing {
        &self.ring
    }

    /// Take the owned polynomials.
    pub fn into_polynomials(self) -> Vec<Polynomial> {
        self.polynomials
    }
}

impl Deref for GroebnerBasis {
    type Target = [Polynomial];

    fn deref(&self) -> &[Polynomial] {
        &self.polynomials
    }
}

impl IntoIterator for GroebnerBasis {
    type Item = Polynomial;
    type IntoIter = std::vec::IntoIter<Polynomial>;

    fn into_iter(self) -> Self::IntoIter {
        self.polynomials.into_iter()
    }
}

impl<'a> IntoIterator for &'a GroebnerBasis {
    type Item = &'a Polynomial;
    type IntoIter = std::slice::Iter<'a, Polynomial>;

    fn into_iter(self) -> Self::IntoIter {
        self.polynomials.iter()
    }
}
