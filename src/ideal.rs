//! Ideals and the bases computed from them.

use std::ops::Deref;

use crate::certificate::{CertifiedGroebnerBasis, CertifyError};
use crate::compute::{
    self, Budget, ComputeError, ComputeLimits, ComputeOptions, ComputeReport, RationalOptions,
};
use crate::hilbert::{self, HilbertError, HilbertSeries};
use crate::normal_form::{self, BasisError, NormalFormError};
use crate::poly::Polynomial;
use crate::ring::{Domain, ModularLift, PolynomialRing, PrimeField, RationalMeta, Rationals};

/// The ideal a list of polynomials generates.
///
/// Build one with [`PolynomialRing::ideal`]. The domain parameter defaults
/// to [`PrimeField`].
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
pub struct Ideal<D: Domain = PrimeField> {
    ring: PolynomialRing<D>,
    generators: Vec<Polynomial<D>>,
}

impl<D: Domain> Ideal<D> {
    pub(crate) fn new(ring: PolynomialRing<D>, generators: Vec<Polynomial<D>>) -> Self {
        Ideal { ring, generators }
    }

    /// The ring the ideal lives in.
    pub fn ring(&self) -> &PolynomialRing<D> {
        &self.ring
    }

    /// The generators, in the order the caller gave them.
    pub fn generators(&self) -> &[Polynomial<D>] {
        &self.generators
    }

    /// Report whether every generator is homogeneous.
    ///
    /// A homogeneous polynomial holds terms of one total degree, and the
    /// zero polynomial counts as one. The test is sufficient and not
    /// necessary: a homogeneous ideal can be given by inhomogeneous
    /// generators, and then this returns false.
    /// [`GroebnerBasis::is_homogeneous`] decides the ideal.
    pub fn has_homogeneous_generators(&self) -> bool {
        hilbert::all_homogeneous(&self.generators)
    }
}

impl Ideal<PrimeField> {
    /// Compute the reduced Gröbner basis under grevlex.
    ///
    /// The options pick the backend and the resource budget. Nothing else
    /// is checked. The result carries no proof: it is the value the engine
    /// returned, unchecked. Use [`Ideal::groebner_basis_certified`] for a
    /// basis an independent verifier has accepted.
    pub fn groebner_basis(&self, options: ComputeOptions) -> Result<GroebnerBasis, ComputeError> {
        let polynomials = compute::groebner_basis(&self.ring, &self.generators, &options)?;
        Ok(GroebnerBasis::new(self.ring.clone(), polynomials, ()))
    }

    /// Compute the reduced Gröbner basis and report what the run did.
    ///
    /// The basis is the one [`Ideal::groebner_basis`] returns for the same
    /// options, and it carries no proof either. The report adds the
    /// backend, the counters of the run, and the wall time of the engine
    /// call.
    pub fn groebner_basis_with_report(
        &self,
        options: ComputeOptions,
    ) -> Result<(GroebnerBasis, ComputeReport), ComputeError> {
        let (polynomials, report) =
            compute::groebner_basis_with_report(&self.ring, &self.generators, &options)?;
        Ok((
            GroebnerBasis::new(self.ring.clone(), polynomials, ()),
            report,
        ))
    }

    /// Compute the reduced Gröbner basis and the certificate that proves
    /// it.
    ///
    /// The value exists only after the independent verifier in
    /// [`crate::verify`] accepts the certificate bytes. The basis it
    /// carries is decoded from those bytes, not taken from the engine.
    ///
    /// The certificate format follows the backend. The classic backend
    /// writes `sylv-gb-cert-v1` from the cofactors it tracks, and the F4
    /// backend, which is the default, writes `sylv-gb-cert-v2` from the
    /// trace it records. There is no format option.
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
    /// does not hold. That is an engine defect, reported rather than
    /// panicked.
    pub fn groebner_basis_certified(
        &self,
        options: ComputeOptions,
    ) -> Result<CertifiedGroebnerBasis, CertifyError> {
        compute::groebner_basis_certified(&self.ring, &self.generators, &options)
    }
}

impl Ideal<Rationals> {
    /// Check that a candidate is the reduced basis of this input ideal.
    ///
    /// The check tracks exact origin cofactors over `Q` and tests both
    /// ideal inclusions. It retains the input, basis, and origin witnesses.
    /// This is a library check, not an independent certificate. It can be
    /// much more expensive than the multimodular computation; a finite
    /// budget is appropriate for inputs of unknown cost.
    pub fn check_basis_equality(
        &self,
        basis: &GroebnerBasis<Rationals>,
        budget: Budget,
    ) -> Result<crate::RationalEqualityCheck, crate::EqualityCheckError> {
        crate::rational_check::check_basis(&self.ring, &self.generators, basis, budget)
    }

    /// Compute the reduced Gröbner basis under grevlex.
    ///
    /// The engine is multimodular: it clears the denominators of the
    /// generators, computes a basis modulo many 31-bit primes with the
    /// backend the options name, combines the runs by the Chinese
    /// remainder theorem, and lifts each residue to a rational number.
    ///
    /// The result is a heuristic. No isolated verifier checks the vote
    /// over the leading monomials, the combination, or the lift, and there
    /// is no certified path over `Q`.
    /// [`GroebnerBasis::lift`] carries the counters of the run and
    /// [`Established`], which names what the stopping rule observed.
    /// [`RationalStop::ContainsInput`] is the stronger of the two rules,
    /// and it still does not establish that the basis generates the input
    /// ideal.
    ///
    /// ```
    /// use sylvester::{Established, PolynomialRing, RationalOptions};
    ///
    /// let ring = PolynomialRing::rationals(["x", "y", "z"])?;
    /// let ideal = ring.ideal([
    ///     ring.parse_polynomial("x + y + z")?,
    ///     ring.parse_polynomial("x*y + y*z + z*x")?,
    ///     ring.parse_polynomial("x*y*z - 1")?,
    /// ])?;
    /// let basis = ideal.groebner_basis(RationalOptions::new())?;
    /// assert_eq!(basis[0].to_string(), "z^3 - 1");
    /// assert_eq!(
    ///     basis.lift().map(|lift| lift.established),
    ///     Some(Established::Unchanged)
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// [`Established`]: crate::Established
    /// [`RationalStop::ContainsInput`]: crate::RationalStop::ContainsInput
    pub fn groebner_basis(
        &self,
        options: RationalOptions,
    ) -> Result<GroebnerBasis<Rationals>, ComputeError> {
        let (polynomials, lift) =
            compute::rational_groebner_basis(&self.ring, &self.generators, &options)?;
        Ok(GroebnerBasis::new(
            self.ring.clone(),
            polynomials,
            RationalMeta::Lifted(lift),
        ))
    }

    /// Compute the reduced Gröbner basis and report what the run did.
    ///
    /// The basis is the one [`Ideal::<Rationals>::groebner_basis`] returns
    /// for the same options. The report carries no F4 counters: they
    /// describe one prime run, a rational run has many, and summing them
    /// would name a run that never happened.
    /// [`ComputeReport::modular`] carries the record of the run and
    /// [`ComputeReport::modular_concurrency`] the number of prime runs the
    /// driver ran at once.
    ///
    /// [`Ideal::<Rationals>::groebner_basis`]: Ideal::groebner_basis
    pub fn groebner_basis_with_report(
        &self,
        options: RationalOptions,
    ) -> Result<(GroebnerBasis<Rationals>, ComputeReport), ComputeError> {
        let (polynomials, lift, report) =
            compute::rational_groebner_basis_with_report(&self.ring, &self.generators, &options)?;
        Ok((
            GroebnerBasis::new(self.ring.clone(), polynomials, RationalMeta::Lifted(lift)),
            report,
        ))
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
///
/// The domain parameter defaults to [`PrimeField`]. Two bases are equal
/// when they hold the same ring and the same polynomials in the same
/// order. What the domain records about the origin of the basis is not
/// part of that, so two computations of one basis stay equal.
#[derive(Clone, Debug)]
pub struct GroebnerBasis<D: Domain = PrimeField> {
    ring: PolynomialRing<D>,
    polynomials: Vec<Polynomial<D>>,
    /// What the domain records about where the basis came from. It is `()`
    /// over [`PrimeField`] and [`RationalMeta`] over [`Rationals`].
    meta: D::BasisMeta,
}

impl<D: Domain> GroebnerBasis<D> {
    pub(crate) fn new(
        ring: PolynomialRing<D>,
        polynomials: Vec<Polynomial<D>>,
        meta: D::BasisMeta,
    ) -> Self {
        GroebnerBasis {
            ring,
            polynomials,
            meta,
        }
    }

    /// The ring the basis lives in.
    pub fn ring(&self) -> &PolynomialRing<D> {
        &self.ring
    }

    /// Take the owned polynomials.
    pub fn into_polynomials(self) -> Vec<Polynomial<D>> {
        self.polynomials
    }

    /// Reduce `f` modulo the basis and return the remainder.
    ///
    /// The remainder holds no monomial divisible by the leading monomial
    /// of a basis element. The basis is a Gröbner basis, so the remainder
    /// is the unique normal form of `f` modulo the ideal, and two calls
    /// with the same arguments return the same value.
    ///
    /// The budget stops the division between two reduction steps, and
    /// [`NormalFormError`] reports what stopped it. A polynomial of
    /// another ring is [`NormalFormError::RingMismatch`]. Over `Q` the
    /// coefficients are exact and can grow, which is what the budget is
    /// for.
    ///
    /// ```
    /// use sylvester::{Budget, ComputeOptions, PolynomialRing};
    ///
    /// let ring = PolynomialRing::prime_field(32003, ["x", "y"])?;
    /// let f = ring.parse_polynomial("x^2 - 1")?;
    /// let ideal = ring.ideal([f.clone(), ring.parse_polynomial("x*y - 1")?])?;
    /// let basis = ideal.groebner_basis(ComputeOptions::new())?;
    /// assert!(basis.normal_form(&f, Budget::new())?.is_zero());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn normal_form(
        &self,
        f: &Polynomial<D>,
        budget: Budget,
    ) -> Result<Polynomial<D>, NormalFormError> {
        if f.ring() != &self.ring {
            return Err(NormalFormError::RingMismatch);
        }
        let limits = ComputeLimits::of_budget(&budget);
        normal_form::normal_form(&self.polynomials, f, &limits).map_err(NormalFormError::from)
    }

    /// Report whether `f` belongs to the ideal this basis generates.
    ///
    /// A polynomial belongs to that ideal exactly when its normal form is
    /// zero, so this is [`GroebnerBasis::normal_form`] with the remainder
    /// tested for zero. It carries the same budget and the same errors.
    ///
    /// The answer is about the ideal this basis generates, which is not
    /// always the ideal the caller started from.
    /// [`GroebnerBasis::from_polynomials`] takes the basis whole, so
    /// there the two ideals are the same by construction. Over a prime
    /// field the engine returns a basis of the input ideal, and
    /// [`Ideal::groebner_basis_certified`] is the path an isolated
    /// verifier stands behind. Over `Q` the driver is a heuristic, and
    /// [`GroebnerBasis::lift`] says what holds:
    ///
    /// - Under [`Established::Unchanged`] nothing relates the two ideals.
    ///   Neither answer is established for the input ideal.
    /// - Under [`Established::ContainsInput`] the ideal of this basis
    ///   contains the input ideal. `false` then holds for the input ideal
    ///   as well, because a polynomial outside the larger ideal is
    ///   outside the smaller one. `true` does not, because the polynomial
    ///   can lie in the larger ideal alone.
    ///
    /// [`Established::Unchanged`]: crate::Established::Unchanged
    /// [`Established::ContainsInput`]: crate::Established::ContainsInput
    pub fn contains(&self, f: &Polynomial<D>, budget: Budget) -> Result<bool, NormalFormError> {
        Ok(self.normal_form(f, budget)?.is_zero())
    }

    /// The Hilbert series of the quotient by the leading monomial ideal
    /// of this basis.
    ///
    /// The value is `N(t) / (1 - t)^n`, with `n` the number of variables
    /// of the ring. [`HilbertSeries`] states what it says about the ideal
    /// the basis generates: it is that ideal's own Hilbert series when the
    /// ideal is homogeneous, and the first difference of its affine
    /// Hilbert function otherwise.
    ///
    /// The computation reads leading monomials alone, so the same ideal
    /// over `F_p` and over `Q` gives the same series.
    ///
    /// Every statement above is about the ideal this basis generates,
    /// which is not always the ideal the caller started from. Over `Q`,
    /// [`GroebnerBasis::lift`] says what relates the two: under
    /// [`Established::Unchanged`] nothing does, and under
    /// [`Established::ContainsInput`] the series is that of an ideal
    /// containing the input ideal.
    ///
    /// [`Established::Unchanged`]: crate::Established::Unchanged
    /// [`Established::ContainsInput`]: crate::Established::ContainsInput
    ///
    /// The recursion behind the numerator branches, and this release
    /// offers no bound on its work, so the budget stops it at a recursion
    /// node and [`HilbertError`] reports what stopped it. The recursion
    /// runs on a heap allocated worklist, which the budget charges, and
    /// not on the calling stack: its depth reaches the sum of the degrees
    /// of the minimal generators of the leading monomial ideal, which the
    /// calling stack does not hold.
    ///
    /// ```
    /// use sylvester::{Budget, ComputeOptions, PolynomialRing};
    ///
    /// let ring = PolynomialRing::prime_field(32003, ["x", "y"])?;
    /// let ideal = ring.ideal([
    ///     ring.parse_polynomial("x^2")?,
    ///     ring.parse_polynomial("x*y")?,
    /// ])?;
    /// let series = ideal
    ///     .groebner_basis(ComputeOptions::new())?
    ///     .hilbert_series(Budget::new())?;
    /// assert_eq!(series.to_string(), "(1 - 2*t^2 + t^3)/(1 - t)^2");
    /// assert_eq!(series.dimension(), Some(1));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn hilbert_series(&self, budget: Budget) -> Result<HilbertSeries, HilbertError> {
        let limits = ComputeLimits::of_budget(&budget);
        hilbert::hilbert_series(&self.polynomials, self.ring.nvars(), &limits)
    }

    /// The Krull dimension of the quotient by the ideal this basis
    /// generates, or `None` for the unit ideal.
    ///
    /// The value is the one [`HilbertSeries::dimension`] reads off
    /// [`GroebnerBasis::hilbert_series`], and it carries the same budget
    /// and the same errors. Passing to the leading monomial ideal keeps
    /// the dimension, so the value is the dimension of the quotient by
    /// the ideal this basis generates.
    ///
    /// That ideal is not always the ideal the caller started from. Over
    /// `Q`, under [`Established::Unchanged`] nothing relates the two.
    /// Under [`Established::ContainsInput`] the ideal of this basis
    /// contains the input ideal, so the value is at most the Krull
    /// dimension of the quotient by the input ideal.
    ///
    /// `None` is the unit ideal, whose quotient ring is zero. Its
    /// dimension is -1 by the usual convention, which the return type
    /// cannot hold.
    ///
    /// [`Established::Unchanged`]: crate::Established::Unchanged
    /// [`Established::ContainsInput`]: crate::Established::ContainsInput
    pub fn krull_dimension(&self, budget: Budget) -> Result<Option<usize>, HilbertError> {
        Ok(self.hilbert_series(budget)?.dimension())
    }

    /// Report whether the ideal the basis generates is homogeneous.
    ///
    /// The test is every element of the basis, and it decides the ideal:
    /// under a graded order an ideal is homogeneous exactly when its
    /// reduced Gröbner basis is. Grevlex is graded.
    pub fn is_homogeneous(&self) -> bool {
        hilbert::all_homogeneous(&self.polynomials)
    }
}

impl GroebnerBasis<PrimeField> {
    /// Take a basis the caller supplies, after checking that it is one.
    ///
    /// The check runs under the budget: every polynomial belongs to
    /// `ring`, none is zero, each is monic, the leading monomials run
    /// strictly descending, no monomial is divisible by the leading
    /// monomial of another element, and every S-polynomial reduces to
    /// zero. [`BasisError`] names the first failure.
    ///
    /// A checked basis is a check and not a certificate. It establishes
    /// what the list is, and nothing about the ideal the caller meant.
    /// [`Ideal::groebner_basis_certified`] is the only path an
    /// independent verifier stands behind.
    ///
    /// ```
    /// use sylvester::{Budget, GroebnerBasis, PolynomialRing, PrimeField};
    ///
    /// let ring = PolynomialRing::prime_field(32003, ["x", "y"])?;
    /// let basis = GroebnerBasis::<PrimeField>::from_polynomials(
    ///     &ring,
    ///     vec![ring.parse_polynomial("y^2 - 1")?, ring.parse_polynomial("x - y")?],
    ///     Budget::new(),
    /// )?;
    /// assert_eq!(basis.len(), 2);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn from_polynomials(
        ring: &PolynomialRing<PrimeField>,
        polynomials: Vec<Polynomial<PrimeField>>,
        budget: Budget,
    ) -> Result<Self, BasisError> {
        let limits = ComputeLimits::of_budget(&budget);
        normal_form::check_basis(ring, &polynomials, &limits)?;
        Ok(GroebnerBasis::new(ring.clone(), polynomials, ()))
    }
}

impl GroebnerBasis<Rationals> {
    /// Take a basis the caller supplies, after checking that it is one.
    ///
    /// The checks are the ones
    /// [`GroebnerBasis::<PrimeField>::from_polynomials`] runs, over exact
    /// rational arithmetic. The basis records
    /// [`RationalMeta::Checked`], so [`GroebnerBasis::lift`] returns
    /// `None`: no modular run produced it.
    ///
    /// ```
    /// use sylvester::{Budget, GroebnerBasis, PolynomialRing, Rationals};
    ///
    /// let ring = PolynomialRing::rationals(["x", "y", "z"])?;
    /// let basis = GroebnerBasis::<Rationals>::from_polynomials(
    ///     &ring,
    ///     vec![
    ///         ring.parse_polynomial("z^3 - 1")?,
    ///         ring.parse_polynomial("y^2 + y*z + z^2")?,
    ///         ring.parse_polynomial("x + y + z")?,
    ///     ],
    ///     Budget::new(),
    /// )?;
    /// assert!(basis.lift().is_none());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// [`GroebnerBasis::<PrimeField>::from_polynomials`]: GroebnerBasis::from_polynomials
    pub fn from_polynomials(
        ring: &PolynomialRing<Rationals>,
        polynomials: Vec<Polynomial<Rationals>>,
        budget: Budget,
    ) -> Result<Self, BasisError> {
        let limits = ComputeLimits::of_budget(&budget);
        normal_form::check_basis(ring, &polynomials, &limits)?;
        Ok(GroebnerBasis::new(
            ring.clone(),
            polynomials,
            RationalMeta::Checked,
        ))
    }

    /// The record of the multimodular run that produced the basis.
    ///
    /// It is `None` for a basis the caller supplied and the checked
    /// constructor accepted, because no modular run produced one. The
    /// counters describe the run and establish nothing on their own.
    /// [`ModularLift::established`] is what the run establishes.
    pub fn lift(&self) -> Option<&ModularLift> {
        match &self.meta {
            RationalMeta::Lifted(lift) => Some(lift),
            RationalMeta::Checked => None,
        }
    }
}

impl<D: Domain> PartialEq for GroebnerBasis<D> {
    fn eq(&self, other: &Self) -> bool {
        self.ring == other.ring && self.polynomials == other.polynomials
    }
}

impl<D: Domain> Eq for GroebnerBasis<D> {}

impl<D: Domain> Deref for GroebnerBasis<D> {
    type Target = [Polynomial<D>];

    fn deref(&self) -> &[Polynomial<D>] {
        &self.polynomials
    }
}

impl<D: Domain> IntoIterator for GroebnerBasis<D> {
    type Item = Polynomial<D>;
    type IntoIter = std::vec::IntoIter<Polynomial<D>>;

    fn into_iter(self) -> Self::IntoIter {
        self.polynomials.into_iter()
    }
}

impl<'a, D: Domain> IntoIterator for &'a GroebnerBasis<D> {
    type Item = &'a Polynomial<D>;
    type IntoIter = std::slice::Iter<'a, Polynomial<D>>;

    fn into_iter(self) -> Self::IntoIter {
        self.polynomials.iter()
    }
}
