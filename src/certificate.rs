//! Certificates for a computed basis.
//!
//! A certificate carries the input, the basis, and everything the
//! verifier needs to check both ideal inclusions. There are two
//! contracts, one per backend. `sylv-gb-cert-v1`, in
//! `docs/certificate-v1.md`, carries explicit cofactor polynomials and
//! belongs to the classic backend. `sylv-gb-cert-v2`, in
//! `docs/certificate-v2.md`, carries the operation trace of the run and
//! belongs to the F4 backend. The module that writes a certificate is
//! untrusted. The verifier in [`crate::verify`] is the trust boundary.
//! Acceptance proves two facts: the input and the basis generate the same
//! ideal, and the basis is the reduced Gröbner basis of that ideal under
//! grevlex.

use std::fmt;

use crate::compute::ComputeError;
use crate::ideal::GroebnerBasis;
use crate::verify::VerifyError;

/// A basis with the certificate the verifier accepted for it.
///
/// [`crate::Ideal::groebner_basis_certified`] returns one. The basis is
/// decoded from the accepted bytes, so the value the caller reads is the
/// value the verifier checked.
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

    /// The certificate bytes, in the contract of the backend that wrote
    /// them.
    ///
    /// The bytes are a function of the input, the backend, the compute
    /// options, and the target. Two runs of one build over one input, one
    /// backend, and one set of options write the same bytes. A memory
    /// limit is part of that: the F4 engine splits a batch it cannot hold,
    /// and a split batch gives a different trace, so a limit low enough to
    /// split one can change the bytes. The deadline only stops the run.
    /// `docs/certificate-v2.md` section 9 states the rule.
    pub fn certificate(&self) -> &[u8] {
        &self.certificate
    }

    /// Take the owned basis and bytes.
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
/// [`CertifyError::Engine`], [`CertifyError::WriterExhausted`], and
/// [`CertifyError::VerifierExhausted`] report an exhausted budget, one per
/// part of the run. None of them says anything about the basis.
/// [`CertifyError::Emitter`] reports a defect in the candidate the emitter
/// read or in the certificate it wrote. [`CertifyError::InputMismatch`]
/// reports an accepted certificate that does not describe the ideal the
/// caller asked about. [`CertifyError::Rejected`] reports a certificate the
/// verifier read and refused.
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
    /// The certificate writer stopped before it wrote the bytes.
    ///
    /// The deadline or the memory limit of the run ran out while the
    /// writer read the run, divided by the basis, or encoded a section.
    /// The writer wrote no bytes.
    WriterExhausted(ComputeError),
    /// The verifier stopped before it reached a verdict.
    VerifierExhausted(VerifyError),
    /// The verifier rejected the certificate the emitter wrote.
    Rejected(VerifyError),
    /// The certificate would exceed a cap the contract puts on it, so the
    /// writer wrote no bytes. This is exhaustion, not a verdict.
    CapExceeded {
        /// The cap the certificate would exceed.
        cap: CertificateCap,
        /// The value the cap holds.
        limit: usize,
    },
}

/// The part of a recorded run the writer could not read.
///
/// Every value names a defect in the trace the engine reported, never an
/// exhausted budget. The `sylv-gb-cert-v2` writer reports one of these
/// and writes no bytes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraceFault {
    /// A row the kernel installed as a pivot carries no summand.
    EmptyRow,
    /// A row or a pivot slot the kernel named carries no value.
    UnboundRow,
    /// A basis element the run returned carries no value.
    UnboundBasis,
    /// The run returned a basis of one length and the trace names another.
    BasisCount {
        /// The number of elements the trace names.
        found: usize,
        /// The number of elements the run returned.
        expected: usize,
    },
}

/// A cap the certificate contract puts on one certificate.
///
/// The writer stops rather than write bytes above a cap, because the
/// verifier applies the same cap before it allocates. The names are the
/// ones `docs/certificate-v2.md` section 8 uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CertificateCap {
    /// Monomials in the pool.
    PoolMonomials,
    /// Variable and exponent pairs in the pool.
    PoolEntries,
    /// Polynomials in the input section.
    InputPolys,
    /// Terms in one polynomial.
    TermsPerPoly,
    /// Terms in the certificate.
    TotalTerms,
    /// Nodes in the trace.
    Nodes,
    /// Steps in one `Comb` node.
    CombSteps,
    /// `Comb` steps in the whole trace.
    TraceSteps,
    /// Elements of the basis.
    Basis,
    /// Pairs of the basis.
    Pairs,
    /// Steps in one division trace.
    DivisionSteps,
    /// Division steps in the whole certificate.
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
            TraceFault::EmptyRow => {
                f.write_str("a row the kernel installed as a pivot carries no summand")
            }
            TraceFault::UnboundRow => f.write_str("a row the kernel named carries no value"),
            TraceFault::UnboundBasis => {
                f.write_str("a basis element the run returned carries no value")
            }
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
    /// A basis element is not monic, so no division trace over the basis
    /// holds.
    BasisElementNotMonic {
        /// The index of the element.
        index: usize,
    },
    /// The trace the engine recorded does not describe the basis it
    /// returned.
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
            CertifyError::WriterExhausted(error) => {
                write!(f, "the certificate writer stopped: {error}")
            }
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
            CertifyError::Engine(error) | CertifyError::WriterExhausted(error) => Some(error),
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
