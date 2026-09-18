//! Exact rational checks for saved computation records.

use std::mem::size_of;

use crate::compute::{Budget, ComputeLimits};
use crate::ideal::GroebnerBasis;
use crate::rational_check::{self, RationalEqualityCheck};
use crate::ring::{PolynomialRing, RationalMeta, Rationals};

use super::{EnvelopeDomain, EnvelopeError, ResultEnvelope, check, retained_bytes_with_limits};

impl ResultEnvelope {
    /// Check the saved input and basis for exact equality over `Q`.
    ///
    /// The record's provenance claim has no effect. The input and basis are
    /// parsed into one rational ring, then checked under one absolute
    /// deadline and one cumulative memory cap. Rational records carry no
    /// certificate contract.
    pub fn check_rational(&self, budget: Budget) -> Result<RationalEqualityCheck, EnvelopeError> {
        let limits = ComputeLimits::of_budget(&budget);
        check(&limits, 0)?;
        self.check_rational_header()?;
        let record_bytes = self.heap_bytes_with_limits(&limits)?;
        let ring_bytes = ring_metadata_bytes(&self.variables, &limits)?;
        let external_bytes = record_bytes.saturating_add(ring_bytes);
        check(&limits, external_bytes)?;
        let ring = self.rational_ring()?;
        let mut held = external_bytes;
        let input = parse_polynomials(&ring, &self.input, &limits, &mut held)?;
        let basis = parse_polynomials(&ring, &self.basis, &limits, &mut held)?;
        let candidate = GroebnerBasis::new(ring.clone(), basis, RationalMeta::Checked);
        rational_check::check_basis_with_limits(&ring, &input, &candidate, &limits, external_bytes)
            .map_err(EnvelopeError::EqualityCheck)
    }

    fn check_rational_header(&self) -> Result<(), EnvelopeError> {
        if !matches!(self.domain, EnvelopeDomain::Rationals) {
            return Err(EnvelopeError::Format(
                "exact equality checks require a rational record".into(),
            ));
        }
        if self.schema != super::SCHEMA || self.order != super::ORDER {
            return Err(EnvelopeError::Format(
                "unsupported schema or monomial order".into(),
            ));
        }
        self.check_claim_shape()?;
        if self.certificate.is_some() {
            return Err(EnvelopeError::Format(
                "rational records have no certificate contract".into(),
            ));
        }
        Ok(())
    }

    fn rational_ring(&self) -> Result<PolynomialRing<Rationals>, EnvelopeError> {
        PolynomialRing::rationals(&self.variables).map_err(EnvelopeError::Ring)
    }
}

fn ring_metadata_bytes(
    variables: &[String],
    limits: &ComputeLimits,
) -> Result<usize, EnvelopeError> {
    check(limits, 0)?;
    let mut bytes = size_of::<PolynomialRing<Rationals>>()
        .saturating_add(size_of::<Vec<String>>())
        .saturating_add(size_of::<String>().saturating_mul(variables.len()));
    for (index, name) in variables.iter().enumerate() {
        if index.is_multiple_of(crate::compute::TICK) {
            check(limits, bytes)?;
        }
        bytes = bytes.saturating_add(name.capacity());
    }
    check(limits, bytes)?;
    Ok(bytes)
}

fn parse_polynomials(
    ring: &PolynomialRing<Rationals>,
    expressions: &[String],
    limits: &ComputeLimits,
    held: &mut usize,
) -> Result<Vec<crate::Polynomial<Rationals>>, EnvelopeError> {
    reserve_polynomial_vec(expressions.len(), limits, held)?;
    let mut polynomials = Vec::with_capacity(expressions.len());
    for (index, expression) in expressions.iter().enumerate() {
        if index.is_multiple_of(crate::compute::TICK) {
            check(limits, *held)?;
        }
        let parsing = residual_limits(limits, *held)?;
        let polynomial = ring
            .parse_polynomial_with_limits(expression, &parsing)
            .map_err(EnvelopeError::Expression)?;
        *held = held.saturating_add(retained_bytes_with_limits(&polynomial, limits)?);
        check(limits, *held)?;
        polynomials.push(polynomial);
    }
    Ok(polynomials)
}

fn reserve_polynomial_vec(
    count: usize,
    limits: &ComputeLimits,
    held: &mut usize,
) -> Result<(), EnvelopeError> {
    check(limits, *held)?;
    let allocation = size_of::<Vec<crate::Polynomial<Rationals>>>()
        .saturating_add(size_of::<crate::Polynomial<Rationals>>().saturating_mul(count));
    *held = held.saturating_add(allocation);
    check(limits, *held)
}

fn residual_limits(limits: &ComputeLimits, held: usize) -> Result<ComputeLimits, EnvelopeError> {
    let mut local = limits.clone();
    if let Some(cap) = limits.memory {
        local.memory = Some(
            cap.checked_sub(held)
                .ok_or(EnvelopeError::MemoryLimitExceeded)?,
        );
    }
    Ok(local)
}
