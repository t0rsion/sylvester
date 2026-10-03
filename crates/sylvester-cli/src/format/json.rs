//! The `json` format: a `sylv-result-v1` computation record.

use super::{Body, DomainClaim, Reading};
use sylvester::{Budget, EnvelopeDomain, EnvelopeError, ResultEnvelope};

/// Read a computation record and expose its basis as the polynomial body.
///
/// The record's provenance and certificate remain data. Commands that need
/// an input ideal read it from the record explicitly.
pub fn read(_origin: &str, text: &str, budget: Budget) -> Result<Reading, EnvelopeError> {
    let record = ResultEnvelope::from_json(text.as_bytes(), budget)?;
    let domain = match record.domain() {
        EnvelopeDomain::PrimeField { modulus } => DomainClaim::Prime(*modulus),
        EnvelopeDomain::Rationals => DomainClaim::Rationals,
    };
    let names = record.variables().to_vec();
    let nvars = Some(names.len());
    let body = Body::Expressions(record.basis().to_vec());
    Ok(Reading {
        names: Some(names),
        domain: Some(domain),
        nvars,
        body,
        record: Some(record),
    })
}
