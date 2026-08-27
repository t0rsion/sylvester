//! Independent verifiers for the v1 JSON and v2 binary certificate
//! contracts.
//!
//! The engines are untrusted candidate generators. This module is the trust
//! boundary. It reads certificate bytes and returns a [`VerifiedGb`] only
//! after every obligation of the selected contract holds.
//!
//! The module shares no code with the engines. It carries its own field
//! arithmetic, its own monomial and polynomial types, its own grevlex
//! comparison, its own decoder, and its own S-pair enumeration.
//!
//! Acceptance proves that the input and the basis generate the same ideal,
//! and that the basis is the reduced Gröbner basis of that ideal under
//! grevlex. Acceptance says nothing about the producer of the bytes.
//!
//! The verifier returns a value for every byte string. It rejects with a
//! typed error, and it does not panic. Caps bound the memory it allocates,
//! and a deadline bounds the time it runs. Once the deadline has passed,
//! the verifier reports [`VerifyError::DeadlineExceeded`], never acceptance
//! and never an invalidity it found on the way.
//!
//! ```rust
//! use sylvester::verify::verify;
//!
//! let bytes = br#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"input":[[[1,[1,0]]]],"basis":[[[1,[1,0]]]],"origin":[[[[1,[0,0]]]]],"membership":[[[[1,[0,0]]]]],"spairs":[]}"#;
//! let verified = verify(bytes).expect("the certificate holds");
//! assert_eq!(verified.modulus(), 7);
//! assert_eq!(verified.basis().len(), 1);
//! ```

mod algebra;
mod checks;
mod error;
mod json;
mod limits;
pub mod v2;

pub use algebra::{Exp, Mono, Poly, Term};
pub use error::{
    BasisFault, BinaryFault, Cap, DivisionFault, DivisionSite, Location, NodeFault, PolyFault,
    PoolFault, Syntax, VerifyError, WitnessFault,
};
pub use limits::Limits;

const V1_FIRST: u8 = b'{';

/// A basis that passed every obligation, with the data it was checked
/// against.
///
/// The values come from the certificate, not from an engine.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VerifiedGb {
    modulus: u64,
    nvars: usize,
    input: Vec<Poly>,
    basis: Vec<Poly>,
}

impl VerifiedGb {
    /// The prime the certificate names. The verifier proved it prime.
    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    /// The number of variables.
    pub fn nvars(&self) -> usize {
        self.nvars
    }

    /// The input polynomials, in certificate order.
    pub fn input(&self) -> &[Poly] {
        &self.input
    }

    /// The reduced basis, sorted strictly descending by leading monomial.
    pub fn basis(&self) -> &[Poly] {
        &self.basis
    }

    /// Take the owned input and basis.
    pub fn into_parts(self) -> (Vec<Poly>, Vec<Poly>) {
        (self.input, self.basis)
    }
}

/// Verify certificate bytes under the default caps.
pub fn verify(bytes: &[u8]) -> Result<VerifiedGb, VerifyError> {
    verify_with_limits(bytes, &Limits::default())
}

/// Verify certificate bytes under caller-supplied caps.
///
/// The caps stop the work before the verifier allocates past them. An
/// exhausted verifier reports [`VerifyError::CapExceeded`] or
/// [`VerifyError::DeadlineExceeded`]. Neither value says the certificate is
/// invalid.
pub fn verify_with_limits(bytes: &[u8], limits: &Limits) -> Result<VerifiedGb, VerifyError> {
    match bytes.first().copied() {
        Some(V1_FIRST) => verify_v1(bytes, limits),
        Some(v2::MAGIC_FIRST) => v2::verify_with_limits(bytes, limits),
        first => Err(VerifyError::Format { first }),
    }
}

fn verify_v1(bytes: &[u8], limits: &Limits) -> Result<VerifiedGb, VerifyError> {
    let outcome = json::decode(bytes, limits).and_then(|raw| checks::check(raw, limits));
    // A passed deadline outranks acceptance and invalidity. A cap already
    // reported remains a cap because the verifier stopped at that boundary.
    if !matches!(&outcome, Err(error) if error.is_exhaustion()) {
        limits.check_deadline()?;
    }
    let accepted = outcome?;
    Ok(VerifiedGb {
        modulus: accepted.modulus,
        nvars: accepted.nvars,
        input: accepted.input,
        basis: accepted.basis,
    })
}

/// Convert an accepted v2 value without sharing any verifier arithmetic.
fn from_v2(accepted: v2::Accepted) -> VerifiedGb {
    fn polys(raw: Vec<v2::DensePoly>) -> Vec<Poly> {
        raw.into_iter()
            .map(|terms| {
                Poly::new(
                    terms
                        .into_iter()
                        .map(|(coeff, exps)| Term::new(coeff, Mono::new(exps)))
                        .collect(),
                )
            })
            .collect()
    }
    VerifiedGb {
        modulus: accepted.modulus,
        nvars: accepted.nvars,
        input: polys(accepted.input),
        basis: polys(accepted.basis),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const UNIT_IDEAL: &str = concat!(
        r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
        r#""input":[[[1,[1,0]]],[[1,[1,0]],[1,[0,0]]]],"basis":[[[1,[0,0]]]],"#,
        r#""origin":[[[[6,[0,0]]],[[1,[0,0]]]]],"#,
        r#""membership":[[[[1,[1,0]]]],[[[1,[1,0]],[1,[0,0]]]]],"spairs":[]}"#
    );

    #[test]
    fn verify_accepts_the_unit_ideal_certificate() {
        let verified = verify(UNIT_IDEAL.as_bytes()).expect("the certificate holds");
        assert_eq!(verified.modulus(), 7);
        assert_eq!(verified.nvars(), 2);
        assert_eq!(verified.input().len(), 2);
        assert_eq!(verified.basis().len(), 1);
        assert_eq!(verified.basis()[0].terms()[0].coeff(), 1);
        assert_eq!(verified.basis()[0].terms()[0].mono().degree(), 0);
    }

    #[test]
    fn a_deadline_passed_during_the_last_check_still_stops_acceptance() {
        // Every array is empty, so no loop inside a check or the decoder
        // ever runs and calls its own deadline check. Only a check after
        // every obligation has run can catch a deadline exhausted here.
        const EMPTY: &str = concat!(
            r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":0,"#,
            r#""input":[],"basis":[],"origin":[],"membership":[],"spairs":[]}"#
        );
        let limits = Limits {
            deadline: Some(std::time::Instant::now() - std::time::Duration::from_secs(1)),
            ..Limits::default()
        };
        assert_eq!(
            verify_with_limits(EMPTY.as_bytes(), &limits),
            Err(VerifyError::DeadlineExceeded)
        );
    }

    #[test]
    fn verify_reports_exhaustion_apart_from_invalidity() {
        let limits = Limits {
            max_bytes: 16,
            ..Limits::default()
        };
        let error = verify_with_limits(UNIT_IDEAL.as_bytes(), &limits).expect_err("the cap holds");
        assert!(error.is_exhaustion());
        assert_eq!(
            error,
            VerifyError::CapExceeded {
                cap: Cap::Bytes,
                limit: 16
            }
        );
    }

    #[test]
    fn verify_names_the_basis_element_of_a_failed_origin_identity() {
        let bytes = UNIT_IDEAL.replace(r#""origin":[[[[6,[0,0]]]"#, r#""origin":[[[[5,[0,0]]]"#);
        assert_eq!(
            verify(bytes.as_bytes()),
            Err(VerifyError::OriginIdentity { basis: 0 })
        );
    }
}
