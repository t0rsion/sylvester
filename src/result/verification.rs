//! Match a saved record against independently accepted certificate data.

use std::sync::Arc;

use crate::compute::{Budget, ComputeLimits, verifier_limits};
use crate::ring::{PolynomialRing, PrimeField};
use crate::verify::{self, Limits, VerifiedGb};

use super::{EnvelopeDomain, EnvelopeError, ResultEnvelope, check};

impl ResultEnvelope {
    /// Verify the certificate and match its ring, input, and basis to this record.
    ///
    /// The saved provenance claim has no effect. Both verifier caps and the
    /// call's budget apply. If both limits supply cancellation flags, they
    /// must share the same flag. Variable names label the certificate's
    /// ordered variables; certificates store their count, not their names.
    pub fn verify_prime(&self, caps: &Limits, budget: Budget) -> Result<VerifiedGb, EnvelopeError> {
        let limits = ComputeLimits::of_budget(&budget);
        self.verify_prime_with_limits(caps, &limits)
    }

    pub(super) fn verify_prime_with_limits(
        &self,
        caps: &Limits,
        limits: &ComputeLimits,
    ) -> Result<VerifiedGb, EnvelopeError> {
        let held = self.heap_bytes_with_limits(limits)?.saturating_mul(2);
        check(limits, held)?;
        self.check_header()?;
        let ring = self.prime_ring()?;
        let bytes = self.certificate_bytes()?;
        let caps = combined_caps(caps, limits, bytes.len(), held, ring.nvars())?;
        let verified =
            verify::verify_with_limits(bytes, &caps).map_err(EnvelopeError::Verification)?;
        self.match_verified(&ring, &verified, limits, held)?;
        caps.check_deadline().map_err(EnvelopeError::Verification)?;
        Ok(verified)
    }

    fn certificate_bytes(&self) -> Result<&[u8], EnvelopeError> {
        self.certificate()
            .ok_or_else(|| EnvelopeError::Format("the record has no certificate".into()))
    }

    fn prime_ring(&self) -> Result<PolynomialRing, EnvelopeError> {
        match self.domain {
            EnvelopeDomain::PrimeField { modulus } => {
                PolynomialRing::prime_field(modulus, &self.variables).map_err(EnvelopeError::Ring)
            }
            EnvelopeDomain::Rationals => Err(EnvelopeError::Format(
                "rational records have no certificate contract".into(),
            )),
        }
    }

    fn match_verified(
        &self,
        ring: &PolynomialRing,
        verified: &VerifiedGb,
        limits: &ComputeLimits,
        held: usize,
    ) -> Result<(), EnvelopeError> {
        if ring.modulus() != verified.modulus() || ring.nvars() != verified.nvars() {
            return Err(EnvelopeError::CertificateMismatch);
        }
        let held = held.saturating_add(verified_bytes_with_limits(verified, limits)?);
        check(limits, held)?;
        let mut parsing = limits.clone();
        parsing.memory = parsing.memory.map(|cap| cap.saturating_sub(held));
        match_polynomials(ring, &self.input, verified.input(), &parsing)?;
        match_polynomials(ring, &self.basis, verified.basis(), &parsing)?;
        check(limits, held)
    }
}

fn match_polynomials(
    ring: &PolynomialRing,
    expressions: &[String],
    accepted: &[verify::Poly],
    limits: &ComputeLimits,
) -> Result<(), EnvelopeError> {
    check(limits, 0)?;
    if expressions.len() != accepted.len() {
        return Err(EnvelopeError::CertificateMismatch);
    }
    for (text, expected) in expressions.iter().zip(accepted) {
        let polynomial = ring
            .parse_polynomial_with_limits(text, limits)
            .map_err(EnvelopeError::Expression)?;
        if !same_polynomial(&polynomial, expected, limits)? {
            return Err(EnvelopeError::CertificateMismatch);
        }
    }
    Ok(())
}

fn same_polynomial(
    polynomial: &crate::Polynomial<PrimeField>,
    expected: &verify::Poly,
    limits: &ComputeLimits,
) -> Result<bool, EnvelopeError> {
    let terms = polynomial.terms();
    if terms.len() != expected.terms().len() {
        return Ok(false);
    }
    let bytes = polynomial.retained_bytes();
    for (index, ((coefficient, exponents), term)) in terms.zip(expected.terms()).enumerate() {
        if index.is_multiple_of(crate::compute::TICK) {
            check(limits, bytes)?;
        }
        if coefficient.value() != term.coeff()
            || !exponents.iter().map(|&value| u32::from(value)).eq(term
                .mono()
                .exps()
                .iter()
                .copied())
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn verified_bytes_with_limits(
    verified: &VerifiedGb,
    limits: &ComputeLimits,
) -> Result<usize, EnvelopeError> {
    check(limits, 0)?;
    let term_bytes = size_of::<verify::Term>()
        .saturating_add(verified.nvars().saturating_mul(size_of::<verify::Exp>()));
    let mut bytes = size_of::<VerifiedGb>();
    for (index, polynomial) in verified.input().iter().chain(verified.basis()).enumerate() {
        if index.is_multiple_of(crate::compute::TICK) {
            check(limits, bytes)?;
        }
        bytes = bytes
            .saturating_add(size_of::<verify::Poly>())
            .saturating_add(term_bytes.saturating_mul(polynomial.terms().len()));
    }
    check(limits, bytes)?;
    Ok(bytes)
}

fn combined_caps(
    caps: &Limits,
    limits: &ComputeLimits,
    certificate_len: usize,
    held: usize,
    nvars: usize,
) -> Result<Limits, EnvelopeError> {
    let mut result = caps.clone();
    let mapped = verifier_limits(certificate_len, held, nvars, limits);
    intersect_memory_caps(&mut result, &mapped);
    result.deadline = result.deadline.into_iter().chain(mapped.deadline).min();
    result.cancellation = combine_flags(&result.cancellation, &mapped.cancellation)?;
    Ok(result)
}

fn combine_flags(
    first: &Option<Arc<std::sync::atomic::AtomicBool>>,
    second: &Option<Arc<std::sync::atomic::AtomicBool>>,
) -> Result<Option<Arc<std::sync::atomic::AtomicBool>>, EnvelopeError> {
    if let (Some(first), Some(second)) = (first, second)
        && !Arc::ptr_eq(first, second)
    {
        return Err(EnvelopeError::Format(
            "the verifier and operation budgets use different cancellation flags".into(),
        ));
    }
    Ok(first.clone().or_else(|| second.clone()))
}

fn intersect_memory_caps(result: &mut Limits, other: &Limits) {
    macro_rules! cap {
        ($($field:ident),+ $(,)?) => { $(result.$field = result.$field.min(other.$field);)+ };
    }
    cap!(
        max_bytes,
        max_polys,
        max_terms_per_poly,
        max_total_terms,
        max_entries,
        max_intermediate_bytes,
        max_pool_monomials,
        max_pool_entries,
        max_input_polys,
        max_nodes,
        max_comb_steps,
        max_trace_steps,
        max_basis,
        max_pairs,
        max_division_steps,
        max_total_division_steps,
        max_live_bytes,
        max_work_units,
    );
}
