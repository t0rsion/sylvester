//! Certificates for a computed basis.
//!
//! A certificate carries the input, the basis, and the evidence needed to
//! check both ideal inclusions. Classic writes the JSON
//! `sylv-gb-cert-v1` cofactor contract; F4 writes the binary
//! `sylv-gb-cert-v2` trace contract. The writers are untrusted. The
//! independent verifier in [`crate::verify`] is the trust boundary.

use std::fmt;

use crate::compute::ComputeError;
use crate::ideal::GroebnerBasis;
use crate::verify::VerifyError;

/// A basis with the certificate the verifier accepted for it.
///
/// [`crate::Ideal::groebner_basis_certified`] returns one. The basis is
/// decoded from the accepted bytes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CertifiedGroebnerBasis {
    basis: GroebnerBasis,
    certificate: Vec<u8>,
}

impl CertifiedGroebnerBasis {
    pub(crate) fn new(basis: GroebnerBasis, certificate: Vec<u8>) -> Self {
        CertifiedGroebnerBasis { basis, certificate }
    }

    /// The reduced Gröbner basis the certificate carries.
    pub fn basis(&self) -> &GroebnerBasis {
        &self.basis
    }

    /// The certificate bytes, in the selected backend's contract.
    ///
    /// A certified F4 run records on one thread, so its thread count
    /// changes no byte.
    pub fn certificate(&self) -> &[u8] {
        &self.certificate
    }

    /// The basis and certificate bytes, consuming the value.
    pub fn into_parts(self) -> (GroebnerBasis, Vec<u8>) {
        (self.basis, self.certificate)
    }
}

/// The place a polynomial holds in the data the emitter reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Place {
    /// Input polynomial at this index.
    Input(usize),
    /// Basis polynomial at this index.
    Basis(usize),
    /// One cofactor of an origin identity.
    Origin {
        /// The basis element the origin belongs to.
        basis: usize,
        /// The input polynomial the cofactor multiplies.
        input: usize,
    },
}

impl fmt::Display for Place {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Place::Input(index) => write!(f, "input polynomial {index}"),
            Place::Basis(index) => write!(f, "basis polynomial {index}"),
            Place::Origin { basis, input } => {
                write!(f, "origin cofactor of basis {basis} for input {input}")
            }
        }
    }
}

/// Why certification stops.
///
/// [`CertifyError::Engine`] and [`CertifyError::VerifierExhausted`] report
/// an exhausted budget. Neither says anything about the basis.
/// [`CertifyError::Emitter`] reports a defect in the candidate or the
/// written certificate. [`CertifyError::InputMismatch`] reports an
/// accepted certificate that does not describe the caller's ideal.
/// [`CertifyError::Rejected`] reports a certificate the verifier refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CertifyError {
    /// The emitter could not write a certificate for the candidate it read.
    Emitter(EmitterFault),
    /// The accepted certificate does not describe the ideal the caller
    /// asked about. The emitter wrote the wrong prime, the wrong variable
    /// count, or the wrong input.
    InputMismatch,
    /// The engine stopped before it produced a basis.
    Engine(ComputeError),
    /// The verifier stopped before it reached a verdict.
    VerifierExhausted(VerifyError),
    /// The verifier rejected the certificate the emitter wrote.
    Rejected(VerifyError),
    /// The certificate would exceed a format cap, so the writer stopped.
    CapExceeded {
        /// The cap the certificate would exceed.
        cap: CertificateCap,
        /// The value the cap holds.
        limit: usize,
    },
}

/// The part of an F4 trace the v2 writer could not read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraceFault {
    /// A pivot row carries no summand.
    EmptyRow,
    /// A named row or pivot has no value.
    UnboundRow,
    /// A returned basis element has no value.
    UnboundBasis,
    /// The trace and returned basis have different lengths.
    BasisCount {
        /// The number of elements named by the trace.
        found: usize,
        /// The number of elements returned by the engine.
        expected: usize,
    },
}

/// A hard count cap in the v2 certificate contract.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CertificateCap {
    /// Monomials in the pool.
    PoolMonomials,
    /// Variable/exponent entries in the pool.
    PoolEntries,
    /// Input polynomials.
    InputPolys,
    /// Terms in one polynomial.
    TermsPerPoly,
    /// Terms in the certificate.
    TotalTerms,
    /// Trace nodes.
    Nodes,
    /// Steps in one combination node.
    CombSteps,
    /// Combination steps in the trace.
    TraceSteps,
    /// Basis elements.
    Basis,
    /// Basis pairs.
    Pairs,
    /// Steps in one division trace.
    DivisionSteps,
    /// Division steps in the certificate.
    TotalDivisionSteps,
}

impl fmt::Display for CertificateCap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            CertificateCap::PoolMonomials => "pool monomials",
            CertificateCap::PoolEntries => "pool entries",
            CertificateCap::InputPolys => "input polynomials",
            CertificateCap::TermsPerPoly => "terms in one polynomial",
            CertificateCap::TotalTerms => "terms in the certificate",
            CertificateCap::Nodes => "trace nodes",
            CertificateCap::CombSteps => "steps in one combination",
            CertificateCap::TraceSteps => "combination steps in the trace",
            CertificateCap::Basis => "basis elements",
            CertificateCap::Pairs => "basis pairs",
            CertificateCap::DivisionSteps => "steps in one division trace",
            CertificateCap::TotalDivisionSteps => "division steps in the certificate",
        };
        f.write_str(name)
    }
}

impl fmt::Display for TraceFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TraceFault::EmptyRow => f.write_str("a pivot row carries no summand"),
            TraceFault::UnboundRow => f.write_str("a named row carries no value"),
            TraceFault::UnboundBasis => f.write_str("a basis element carries no value"),
            TraceFault::BasisCount { found, expected } => write!(
                f,
                "the trace names {found} basis elements, the run returned {expected}"
            ),
        }
    }
}

/// Why the emitter could not write a certificate for a candidate.
///
/// Every value names a defect in the basis, the origins, or a polynomial
/// the emitter read, never an exhausted budget.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EmitterFault {
    /// A polynomial holds an exponent vector of the wrong width.
    ExponentCount {
        /// Where the polynomial sits.
        at: Place,
        /// The width the polynomial holds.
        found: usize,
        /// The number of variables in the ring.
        expected: usize,
    },
    /// A basis element is the zero polynomial.
    BasisElementZero {
        /// The index of the element.
        index: usize,
    },
    /// The origin list does not hold one entry per basis element.
    OriginCount {
        /// The number of entries.
        found: usize,
        /// The number of basis elements.
        expected: usize,
    },
    /// An origin entry does not hold one cofactor per input polynomial.
    OriginEntryCount {
        /// The basis element the entry belongs to.
        basis: usize,
        /// The number of cofactors in the entry.
        found: usize,
        /// The number of input polynomials.
        expected: usize,
    },
    /// The S-polynomial of this pair does not reduce to zero over the
    /// basis, so the basis is not a Gröbner basis.
    NotAGroebnerBasis {
        /// The first element of the pair.
        i: usize,
        /// The second element of the pair.
        j: usize,
    },
    /// This input polynomial has a nonzero remainder on division by the
    /// basis, so the basis does not generate an ideal that holds the
    /// input.
    InputHasRemainder {
        /// The index of the input polynomial.
        input: usize,
    },
    /// A basis element is not monic, so no v2 division trace holds.
    BasisElementNotMonic {
        /// The index of the element.
        index: usize,
    },
    /// The F4 trace does not describe the returned basis.
    Trace(TraceFault),
}

impl fmt::Display for CertifyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CertifyError::Emitter(fault) => write!(f, "{fault}"),
            CertifyError::InputMismatch => f.write_str(
                "the accepted certificate does not describe the ideal that was asked about",
            ),
            CertifyError::Engine(error) => write!(f, "the engine stopped: {error}"),
            CertifyError::VerifierExhausted(error) => {
                write!(f, "the verifier stopped: {error}")
            }
            CertifyError::Rejected(error) => {
                write!(f, "the verifier rejected the certificate: {error}")
            }
            CertifyError::CapExceeded { cap, limit } => {
                write!(f, "the certificate would hold more than {limit} {cap}")
            }
        }
    }
}

impl fmt::Display for EmitterFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EmitterFault::ExponentCount {
                at,
                found,
                expected,
            } => write!(
                f,
                "{at} holds {found} exponents, the certificate declares {expected}"
            ),
            EmitterFault::BasisElementZero { index } => {
                write!(f, "basis element {index} is zero")
            }
            EmitterFault::OriginCount { found, expected } => write!(
                f,
                "the origin list holds {found} entries, the basis holds {expected} elements"
            ),
            EmitterFault::OriginEntryCount {
                basis,
                found,
                expected,
            } => write!(
                f,
                "origin entry {basis} holds {found} cofactors, the input holds {expected} polynomials"
            ),
            EmitterFault::NotAGroebnerBasis { i, j } => write!(
                f,
                "the S-polynomial of pair ({i},{j}) does not reduce to zero over the basis"
            ),
            EmitterFault::InputHasRemainder { input } => write!(
                f,
                "input polynomial {input} has a nonzero remainder on division by the basis"
            ),
            EmitterFault::BasisElementNotMonic { index } => {
                write!(f, "basis element {index} is not monic")
            }
            EmitterFault::Trace(fault) => write!(f, "{fault}"),
        }
    }
}

impl std::error::Error for CertifyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            CertifyError::Engine(error) => Some(error),
            CertifyError::VerifierExhausted(error) | CertifyError::Rejected(error) => Some(error),
            CertifyError::Emitter(_)
            | CertifyError::InputMismatch
            | CertifyError::CapExceeded { .. } => None,
        }
    }
}

impl From<EmitterFault> for CertifyError {
    fn from(fault: EmitterFault) -> Self {
        CertifyError::Emitter(fault)
    }
}

impl From<TraceFault> for CertifyError {
    fn from(fault: TraceFault) -> Self {
        CertifyError::Emitter(EmitterFault::Trace(fault))
    }
}

impl From<crate::poly::ExponentOverflow> for CertifyError {
    fn from(overflow: crate::poly::ExponentOverflow) -> Self {
        CertifyError::Engine(ComputeError::from(overflow))
    }
}

impl From<ComputeError> for CertifyError {
    fn from(error: ComputeError) -> Self {
        CertifyError::Engine(error)
    }
}

impl From<VerifyError> for CertifyError {
    /// An exhausted verifier has reached no verdict, so it is not a
    /// rejection.
    fn from(error: VerifyError) -> Self {
        if error.is_exhaustion() {
            CertifyError::VerifierExhausted(error)
        } else {
            CertifyError::Rejected(error)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verify::Cap;

    #[test]
    fn an_exhausted_verifier_is_not_a_rejection() {
        let capped = VerifyError::CapExceeded {
            cap: Cap::Bytes,
            limit: 16,
        };
        assert_eq!(
            CertifyError::from(capped.clone()),
            CertifyError::VerifierExhausted(capped)
        );
        assert_eq!(
            CertifyError::from(VerifyError::DeadlineExceeded),
            CertifyError::VerifierExhausted(VerifyError::DeadlineExceeded)
        );

        let wrong = VerifyError::OriginIdentity { basis: 0 };
        assert_eq!(
            CertifyError::from(wrong.clone()),
            CertifyError::Rejected(wrong)
        );
    }

    #[test]
    fn an_exhausted_engine_keeps_its_value() {
        assert_eq!(
            CertifyError::from(ComputeError::Timeout),
            CertifyError::Engine(ComputeError::Timeout)
        );
    }
}
