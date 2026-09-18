//! Portable computation records with untrusted provenance claims.
//!
//! A record stores the input as well as the basis. Loading a record never
//! establishes ideal equality or verifies a certificate.

mod io;
mod rational;
mod verification;

use std::fmt;
use std::io::BufReader;

use serde::{Deserialize, Serialize};

use crate::compute::{Budget, ComputeError, ComputeLimits};
use crate::ideal::{GroebnerBasis, Ideal};
use crate::ring::{
    Domain, DomainOps, Established, PolynomialRing, PrimeField, Rationals, RingError,
};

use self::io::{BudgetReader, BudgetWriter, CHUNK};

const SCHEMA: &str = "sylv-result-v1";
const ORDER: &str = "grevlex-v1";

/// The coefficient domain named by a saved record.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum EnvelopeDomain {
    /// A prime field with the stated modulus.
    PrimeField {
        /// The prime modulus.
        modulus: u64,
    },
    /// The rational numbers.
    Rationals,
}

/// An untrusted statement recorded by the producer of a result.
///
/// Deserialization preserves this statement as data. It does not turn
/// the statement into an established fact.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ClaimedProvenance {
    /// The record makes no correctness claim.
    Unverified,
    /// The producer supplied and checked a basis of its own ideal.
    SuppliedBasis,
    /// A rational lift stopped changing. Nothing about the input follows.
    Unchanged,
    /// The producer checked containment of the input ideal in the basis ideal.
    ContainsInput,
    /// The producer checked equality with the input ideal.
    EqualsInput,
    /// The producer reports an accepted prime-field certificate.
    Certified,
}

/// A saved input, basis, and provenance claim under `grevlex-v1`.
///
/// The format is `sylv-result-v1`. Polynomial strings use the ring's text
/// syntax. A certificate is optional and stored as an array of bytes.
/// The record is not a certificate. Its claims remain untrusted after
/// loading, including when it holds certificate bytes.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResultEnvelope {
    schema: String,
    order: String,
    domain: EnvelopeDomain,
    variables: Vec<String>,
    input: Vec<String>,
    basis: Vec<String>,
    claimed_provenance: ClaimedProvenance,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    certificate: Option<Vec<u8>>,
}

/// Why a computation record cannot be constructed, read, or written.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EnvelopeError {
    /// The input and basis belong to different rings.
    RingMismatch,
    /// The record has an unsupported schema, order, or invalid JSON shape.
    Format(String),
    /// The record names an invalid ring.
    Ring(RingError),
    /// A saved polynomial expression is invalid or exceeds its budget.
    Expression(crate::ExpressionError),
    /// The attached certificate is invalid or exhausts its verifier caps.
    Verification(crate::verify::VerifyError),
    /// The exact rational equality check rejected the claim or exhausted its budget.
    EqualityCheck(crate::EqualityCheckError),
    /// Accepted certificate data differs from the saved input or basis.
    CertificateMismatch,
    /// The operation passed its deadline or was cancelled.
    Timeout,
    /// The operation passed its memory estimate limit.
    MemoryLimitExceeded,
}

impl fmt::Display for EnvelopeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RingMismatch => f.write_str("the input and basis belong to different rings"),
            Self::Format(message) => write!(f, "invalid computation record: {message}"),
            Self::Ring(error) => error.fmt(f),
            Self::Expression(error) => error.fmt(f),
            Self::Verification(error) => error.fmt(f),
            Self::EqualityCheck(error) => error.fmt(f),
            Self::CertificateMismatch => {
                f.write_str("the certificate does not match the saved result")
            }
            Self::Timeout => f.write_str("the computation record operation stopped"),
            Self::MemoryLimitExceeded => {
                f.write_str("the computation record passed its memory limit")
            }
        }
    }
}

impl std::error::Error for EnvelopeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Ring(error) => Some(error),
            Self::Expression(error) => Some(error),
            Self::Verification(error) => Some(error),
            Self::EqualityCheck(error) => Some(error),
            Self::RingMismatch
            | Self::Format(_)
            | Self::CertificateMismatch
            | Self::Timeout
            | Self::MemoryLimitExceeded => None,
        }
    }
}

impl ResultEnvelope {
    /// Record a prime-field computation and optional unverified certificate bytes.
    ///
    /// This does not verify the bytes or establish input ideal equality.
    pub fn from_prime(
        input: &Ideal<PrimeField>,
        basis: &GroebnerBasis<PrimeField>,
        certificate: Option<&[u8]>,
        budget: Budget,
    ) -> Result<Self, EnvelopeError> {
        let limits = ComputeLimits::of_budget(&budget);
        Self::build(
            input,
            basis,
            EnvelopeDomain::PrimeField {
                modulus: input.ring().modulus(),
            },
            ClaimedProvenance::Unverified,
            certificate,
            &limits,
        )
    }

    /// Record a certified basis after matching the certificate to this input.
    ///
    /// The certificate is verified again to bind its input to the saved
    /// statement. An exhausted budget does not produce a record.
    pub fn from_certified(
        input: &Ideal<PrimeField>,
        certified: &crate::CertifiedGroebnerBasis,
        budget: Budget,
    ) -> Result<Self, EnvelopeError> {
        let limits = ComputeLimits::of_budget(&budget);
        let mut record = Self::build(
            input,
            certified.basis(),
            EnvelopeDomain::PrimeField {
                modulus: input.ring().modulus(),
            },
            ClaimedProvenance::Unverified,
            Some(certified.certificate()),
            &limits,
        )?;
        let external = source_bytes_with_limits(input, certified.basis(), &limits)?
            .saturating_add(certified.certificate().len());
        check(
            &limits,
            external.saturating_add(record.heap_bytes_with_limits(&limits)?),
        )?;
        let caps = crate::verify::Limits {
            deadline: limits.deadline,
            cancellation: limits.cancel.clone(),
            ..Default::default()
        };
        let mut verification_limits = limits.clone();
        if let Some(cap) = limits.memory {
            verification_limits.memory = Some(cap.saturating_sub(external));
        }
        record.verify_prime_with_limits(&caps, &verification_limits)?;
        record.claimed_provenance = ClaimedProvenance::Certified;
        Ok(record)
    }

    /// Record a rational computation with the producer's stopping claim.
    ///
    /// Neither existing rational stopping rule establishes input equality.
    pub fn from_rational(
        input: &Ideal<Rationals>,
        basis: &GroebnerBasis<Rationals>,
        budget: Budget,
    ) -> Result<Self, EnvelopeError> {
        let limits = ComputeLimits::of_budget(&budget);
        let provenance = match basis.lift().map(|lift| lift.established) {
            None => ClaimedProvenance::SuppliedBasis,
            Some(Established::Unchanged) => ClaimedProvenance::Unchanged,
            Some(Established::ContainsInput) => ClaimedProvenance::ContainsInput,
        };
        Self::build(
            input,
            basis,
            EnvelopeDomain::Rationals,
            provenance,
            None,
            &limits,
        )
    }

    /// Record the input and basis associated with an exact rational check.
    ///
    /// The saved equality statement becomes an untrusted claim on loading.
    /// There is no independent rational certificate in this record.
    pub fn from_checked_rational(
        checked: &crate::RationalEqualityCheck,
        budget: Budget,
    ) -> Result<Self, EnvelopeError> {
        let limits = ComputeLimits::of_budget(&budget);
        Self::build(
            checked.input(),
            checked.basis(),
            EnvelopeDomain::Rationals,
            ClaimedProvenance::EqualsInput,
            None,
            &limits,
        )
    }

    fn build<D: Domain>(
        input: &Ideal<D>,
        basis: &GroebnerBasis<D>,
        domain: EnvelopeDomain,
        claimed_provenance: ClaimedProvenance,
        certificate: Option<&[u8]>,
        limits: &ComputeLimits,
    ) -> Result<Self, EnvelopeError> {
        let held = build_preflight(input, basis, certificate, limits)?;
        let result = Self {
            schema: SCHEMA.into(),
            order: ORDER.into(),
            domain,
            variables: input.ring().variables().to_vec(),
            input: render(input.generators(), limits, held)?,
            basis: render(basis, limits, held)?,
            claimed_provenance,
            certificate: certificate.map(<[u8]>::to_vec),
        };
        check(
            limits,
            held.saturating_add(result.heap_bytes_with_limits(limits)?),
        )?;
        Ok(result)
    }

    /// Read a record without trusting its provenance or certificate.
    ///
    /// The byte count bounds the JSON allocation estimate before decoding.
    /// Ring metadata is validated; polynomial expressions remain data.
    pub fn from_json(bytes: &[u8], budget: Budget) -> Result<Self, EnvelopeError> {
        let limits = ComputeLimits::of_budget(&budget);
        let read_buffer = bytes.len().clamp(1, CHUNK);
        let estimate = bytes
            .len()
            .saturating_mul(32)
            .saturating_add(size_of::<Self>())
            .saturating_add(read_buffer);
        check(&limits, estimate)?;
        let mut reader = BufReader::with_capacity(read_buffer, BudgetReader::new(bytes, &limits));
        let decoded: Result<Self, _> = serde_json::from_reader(&mut reader);
        let record = match decoded {
            Ok(record) => record,
            Err(error) => {
                if let Some(stop) = reader.get_mut().take_stop() {
                    return Err(stop);
                }
                check(&limits, estimate)?;
                return Err(EnvelopeError::Format(error.to_string()));
            }
        };
        check(&limits, estimate)?;
        if let Err(error) = record.check_header() {
            check(&limits, estimate)?;
            return Err(error);
        }
        check(
            &limits,
            estimate.saturating_add(record.heap_bytes_with_limits(&limits)?),
        )?;
        Ok(record)
    }

    /// Write a JSON record under the supplied budget.
    pub fn to_json(&self, budget: Budget) -> Result<Vec<u8>, EnvelopeError> {
        let limits = ComputeLimits::of_budget(&budget);
        let record_bytes = self.heap_bytes_with_limits(&limits)?;
        let estimate = record_bytes.saturating_mul(8).saturating_add(1024);
        check(&limits, estimate)?;
        let mut writer = BudgetWriter::new(&limits, estimate);
        let encoded = serde_json::to_writer_pretty(&mut writer, self);
        if let Some(stop) = writer.take_stop() {
            return Err(stop);
        }
        if let Err(error) = encoded {
            check(&limits, estimate.saturating_add(writer.capacity()))?;
            return Err(EnvelopeError::Format(error.to_string()));
        }
        check(&limits, estimate.saturating_add(writer.capacity()))?;
        Ok(writer.into_inner())
    }

    /// The coefficient domain named by the record.
    pub fn domain(&self) -> &EnvelopeDomain {
        &self.domain
    }

    /// The variable names in ring order.
    pub fn variables(&self) -> &[String] {
        &self.variables
    }

    /// The original input expressions claimed by the record.
    pub fn input(&self) -> &[String] {
        &self.input
    }

    /// The basis expressions claimed by the record.
    pub fn basis(&self) -> &[String] {
        &self.basis
    }

    /// The producer's untrusted claim.
    pub fn claimed_provenance(&self) -> ClaimedProvenance {
        self.claimed_provenance
    }

    /// Unverified certificate bytes, if present.
    pub fn certificate(&self) -> Option<&[u8]> {
        self.certificate.as_deref()
    }

    fn check_header(&self) -> Result<(), EnvelopeError> {
        if self.schema != SCHEMA || self.order != ORDER {
            return Err(EnvelopeError::Format(
                "unsupported schema or monomial order".into(),
            ));
        }
        self.check_claim_shape()?;
        match self.domain {
            EnvelopeDomain::PrimeField { modulus } => {
                PolynomialRing::prime_field(modulus, &self.variables)
                    .map_err(EnvelopeError::Ring)?;
            }
            EnvelopeDomain::Rationals => {
                PolynomialRing::rationals(&self.variables).map_err(EnvelopeError::Ring)?;
                if self.certificate.is_some() {
                    return Err(EnvelopeError::Format(
                        "rational records have no certificate contract".into(),
                    ));
                }
            }
        }
        Ok(())
    }

    fn check_claim_shape(&self) -> Result<(), EnvelopeError> {
        match (&self.domain, self.claimed_provenance) {
            (EnvelopeDomain::Rationals, ClaimedProvenance::Certified)
            | (
                EnvelopeDomain::PrimeField { .. },
                ClaimedProvenance::Unchanged
                | ClaimedProvenance::ContainsInput
                | ClaimedProvenance::EqualsInput,
            ) => Err(EnvelopeError::Format(
                "the provenance claim does not apply to this domain".into(),
            )),
            (_, ClaimedProvenance::Certified) if self.certificate.is_none() => Err(
                EnvelopeError::Format("a certified claim needs certificate bytes".into()),
            ),
            _ => Ok(()),
        }
    }

    pub(super) fn heap_bytes_with_limits(
        &self,
        limits: &ComputeLimits,
    ) -> Result<usize, EnvelopeError> {
        check(limits, 0)?;
        let mut bytes = size_of::<Self>()
            .saturating_add(
                self.variables
                    .capacity()
                    .saturating_mul(size_of::<String>()),
            )
            .saturating_add(self.input.capacity().saturating_mul(size_of::<String>()))
            .saturating_add(self.basis.capacity().saturating_mul(size_of::<String>()))
            .saturating_add(self.certificate.as_ref().map_or(0, Vec::capacity));
        for (index, text) in self
            .variables
            .iter()
            .chain(&self.input)
            .chain(&self.basis)
            .chain(std::iter::once(&self.schema))
            .chain(std::iter::once(&self.order))
            .enumerate()
        {
            if index.is_multiple_of(crate::compute::TICK) {
                check(limits, bytes)?;
            }
            bytes = bytes.saturating_add(text.capacity());
        }
        check(limits, bytes)?;
        Ok(bytes)
    }
}

fn record_metadata_bytes_with_limits<D: Domain>(
    ring: &PolynomialRing<D>,
    certificate: Option<&[u8]>,
    limits: &ComputeLimits,
) -> Result<usize, EnvelopeError> {
    check(limits, 0)?;
    let mut bytes = size_of::<ResultEnvelope>()
        .saturating_add(SCHEMA.len() + ORDER.len())
        .saturating_add(certificate.map_or(0, <[u8]>::len));
    for (index, name) in ring.variables().iter().enumerate() {
        if index.is_multiple_of(crate::compute::TICK) {
            check(limits, bytes)?;
        }
        bytes = bytes
            .saturating_add(size_of::<String>())
            .saturating_add(name.len());
    }
    check(limits, bytes)?;
    Ok(bytes)
}

fn build_preflight<D: Domain>(
    input: &Ideal<D>,
    basis: &GroebnerBasis<D>,
    certificate: Option<&[u8]>,
    limits: &ComputeLimits,
) -> Result<usize, EnvelopeError> {
    check(limits, 0)?;
    if input.ring() != basis.ring() {
        return Err(EnvelopeError::RingMismatch);
    }
    let held = source_bytes_with_limits(input, basis, limits)?
        .saturating_add(certificate.map_or(0, <[u8]>::len));
    check(limits, held)?;
    let rendered =
        rendered_bytes_with_limits(input.generators().iter().chain(basis.iter()), limits)?;
    let metadata = record_metadata_bytes_with_limits(input.ring(), certificate, limits)?;
    check(
        limits,
        held.saturating_add(rendered).saturating_add(metadata),
    )?;
    Ok(held)
}

fn source_bytes_with_limits<D: Domain>(
    input: &Ideal<D>,
    basis: &GroebnerBasis<D>,
    limits: &ComputeLimits,
) -> Result<usize, EnvelopeError> {
    check(limits, 0)?;
    let mut names = 0usize;
    for (index, name) in input.ring().variables().iter().enumerate() {
        if index.is_multiple_of(crate::compute::TICK) {
            check(limits, names)?;
        }
        names = names
            .saturating_add(name.capacity())
            .saturating_add(size_of::<String>());
    }
    let mut bytes = names.saturating_mul(2);
    for (index, polynomial) in input.generators().iter().chain(basis.iter()).enumerate() {
        if index.is_multiple_of(crate::compute::TICK) {
            check(limits, bytes)?;
        }
        bytes = bytes
            .saturating_add(retained_bytes_with_limits(polynomial, limits)?)
            .saturating_add(size_of::<crate::Polynomial<D>>());
    }
    check(limits, bytes)?;
    Ok(bytes)
}

fn rendered_bytes_with_limits<'a, D: Domain + 'a>(
    polynomials: impl Iterator<Item = &'a crate::Polynomial<D>>,
    limits: &ComputeLimits,
) -> Result<usize, EnvelopeError> {
    check(limits, 0)?;
    let mut bytes = 0usize;
    for (index, polynomial) in polynomials.enumerate() {
        if index.is_multiple_of(crate::compute::TICK) {
            check(limits, bytes)?;
        }
        bytes = bytes.saturating_add(render_bound_with_limits(polynomial, limits)?);
    }
    check(limits, bytes)?;
    Ok(bytes)
}

fn render_bound_with_limits<D: Domain>(
    polynomial: &crate::Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<usize, EnvelopeError> {
    let retained = retained_bytes_with_limits(polynomial, limits)?;
    let spare = polynomial
        .terms
        .capacity()
        .saturating_sub(polynomial.terms.len())
        .saturating_mul(size_of::<crate::poly::Term<D>>());
    let names = polynomial.ring().variables().iter().enumerate().try_fold(
        0usize,
        |bytes, (index, name)| {
            if index.is_multiple_of(crate::compute::TICK) {
                check(limits, bytes)?;
            }
            Ok::<usize, EnvelopeError>(bytes.saturating_add(name.len()))
        },
    )?;
    let heap = retained.saturating_sub(spare);
    let bound = heap
        .saturating_mul(8)
        .saturating_add(names.saturating_mul(polynomial.terms.len()))
        .saturating_add(size_of::<String>() + 1);
    check(limits, bound)?;
    Ok(bound)
}

fn render<D: Domain>(
    polynomials: &[crate::Polynomial<D>],
    limits: &ComputeLimits,
    held: usize,
) -> Result<Vec<String>, EnvelopeError> {
    let mut bytes = held.saturating_add(polynomials.len().saturating_mul(size_of::<String>()));
    check(limits, bytes)?;
    let mut strings = Vec::with_capacity(polynomials.len());
    for polynomial in polynomials {
        check(
            limits,
            bytes.saturating_add(render_bound_with_limits(polynomial, limits)?),
        )?;
        let text = polynomial.to_string();
        bytes = bytes.saturating_add(text.capacity());
        check(limits, bytes)?;
        strings.push(text);
    }
    Ok(strings)
}

pub(super) fn retained_bytes_with_limits<D: Domain>(
    polynomial: &crate::Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<usize, EnvelopeError> {
    check(limits, 0)?;
    let mut bytes = size_of::<crate::poly::Term<D>>()
        .saturating_mul(polynomial.terms.capacity())
        .saturating_add(
            crate::poly::heap_exps_bytes(polynomial.ring().nvars())
                .saturating_mul(polynomial.terms.len()),
        );
    for (index, (coefficient, _)) in polynomial.terms().enumerate() {
        if index.is_multiple_of(crate::compute::TICK) {
            check(limits, bytes)?;
        }
        bytes = bytes.saturating_add(polynomial.ring().ops().heap_bytes(coefficient));
    }
    check(limits, bytes)?;
    Ok(bytes)
}

fn check(limits: &ComputeLimits, bytes: usize) -> Result<(), EnvelopeError> {
    if let Some(stop) = limits.stop() {
        return Err(match stop.reported() {
            ComputeError::MemoryLimitExceeded => EnvelopeError::MemoryLimitExceeded,
            _ => EnvelopeError::Timeout,
        });
    }
    if limits.memory.is_some_and(|cap| bytes > cap) {
        return Err(EnvelopeError::MemoryLimitExceeded);
    }
    Ok(())
}
