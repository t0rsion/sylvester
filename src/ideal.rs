//! Ideals and the bases computed from them.

use std::ops::Deref;

use crate::certificate::{CertifiedGroebnerBasis, CertifyError};
use crate::compute::{self, ComputeError, ComputeOptions, ComputeReport};
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

    /// Compute the reduced Gröbner basis.
    ///
    /// The options pick the backend and the resource budget. The result is
    /// the value the engine returned, unchecked.
    ///
    /// [`Ideal::groebner_basis_certified`] returns a basis an independent
    /// verifier has accepted.
    pub fn groebner_basis(&self, options: ComputeOptions) -> Result<GroebnerBasis, ComputeError> {
        let polynomials = compute::groebner_basis(&self.ring, &self.generators, &options)?;
        Ok(GroebnerBasis {
            ring: self.ring.clone(),
            polynomials,
        })
    }

    /// Compute the reduced Gröbner basis and report what the run did.
    ///
    /// The basis is the one [`Ideal::groebner_basis`] returns for the same
    /// options. The report adds the backend, the counters of the run, and
    /// the wall time of the engine call.
    pub fn groebner_basis_with_report(
        &self,
        options: ComputeOptions,
    ) -> Result<(GroebnerBasis, ComputeReport), ComputeError> {
        let (polynomials, report) =
            compute::groebner_basis_with_report(&self.ring, &self.generators, &options)?;
        let basis = GroebnerBasis {
            ring: self.ring.clone(),
            polynomials,
        };
        Ok((basis, report))
    }

    /// Compute the reduced Gröbner basis and the certificate the verifier
    /// accepted for it.
    ///
    /// The value exists only after the independent verifier in
    /// [`crate::verify`] accepts the certificate bytes. The basis it
    /// carries is decoded from those bytes.
    ///
    /// The certificate format follows the backend. Classic writes
    /// `sylv-gb-cert-v1` from the cofactors it tracks. F4 writes
    /// `sylv-gb-cert-v2` from the trace it records. There is no format
    /// option.
    ///
    /// The deadline and the memory limit cover the whole run: the engine,
    /// the certificate, and the verifier. An exhausted budget is
    /// [`CertifyError::Engine`] up to the verifier and
    /// [`CertifyError::VerifierExhausted`] inside it. Neither says
    /// anything about the basis.
    ///
    /// A run that records its trace stays on one thread, so the thread
    /// count changes no byte of a `sylv-gb-cert-v2` certificate.
    ///
    /// [`CertifyError::Rejected`] means the engine wrote a certificate that
    /// does not hold. That is an engine defect.
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
/// The value derefs to `&[Polynomial]`.
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

    /// The polynomials, consuming the basis.
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
