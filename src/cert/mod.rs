//! Certificate emission for the `sylv-gb-cert-v1` contract.
//!
//! This module is untrusted. It reads engine values and writes the
//! canonical bytes of `docs/certificate-v1.md`. The verifier in
//! [`crate::verify`] is the trust boundary: it reads those bytes back and
//! accepts or rejects them. A wrong cofactor here becomes a rejection
//! there, never a silent repair. This module never imports `crate::verify`:
//! [`crate::compute::groebner_basis_certified`] runs the emitter here, then
//! hands the bytes to the verifier itself.
//!
//! The caller supplies the origin cofactors, which the engines track
//! through reduction. The emitter generates the membership and S-pair
//! cofactors by division with recorded quotients, after interreduction.
//!
//! [`assemble`] writes one certificate. It holds to the deadline and the
//! memory limit of the run, and it charges the certificate bytes to that
//! budget before it returns them.

mod divide;
mod json;
pub(crate) mod origin;
mod represent;
pub(crate) mod v2;

use crate::certificate::{CertifyError, EmitterFault, Place};
use crate::compute::{ComputeError, ComputeLimits};
use crate::poly::{Polynomial, Term};
use crate::ring::PolynomialRing;

use origin::Origin;
use represent::{membership_representations, spair_representations};

/// The limits one emission runs under.
///
/// The emitter charges the data it keeps before it builds more. The count
/// covers the input, the basis, the origins, the cofactors the emitter
/// generates, and the certificate bytes. The count estimates the live
/// emitter data. It is not a count of process memory. The limits are the
/// ones the engine ran under, so the deadline the writer holds to is what
/// the engine left of it.
///
/// An exhausted budget here is [`CertifyError::WriterExhausted`]. The
/// writer stopped; the engine did not.
///
/// One term costs its own size, one exponent per variable, and the heap
/// bytes of its coefficient. One polynomial costs the header as well,
/// because a basis of n elements carries n cofactors per S-pair and most
/// of them are zero.
pub(crate) struct WriterBudget {
    limits: ComputeLimits,
    poly_bytes: usize,
    term_bytes: usize,
    held: usize,
}

impl WriterBudget {
    pub(crate) fn new(limits: &ComputeLimits, nvars: usize) -> Self {
        WriterBudget {
            limits: limits.clone(),
            poly_bytes: size_of::<Polynomial>(),
            term_bytes: size_of::<Term>() + nvars * size_of::<u16>(),
            held: 0,
        }
    }

    /// A budget that stops nothing.
    #[cfg(test)]
    pub(crate) fn unlimited(nvars: usize) -> Self {
        WriterBudget::new(&ComputeLimits::default(), nvars)
    }

    /// Return an error once the deadline has passed or the run is
    /// cancelled.
    pub(crate) fn check_stop(&self) -> Result<(), CertifyError> {
        match self.limits.stop() {
            Some(stop) => Err(CertifyError::WriterExhausted(stop.reported())),
            None => Ok(()),
        }
    }

    /// The bytes of `polys` polynomials that hold `terms` terms together.
    fn bytes(&self, polys: usize, terms: usize) -> usize {
        polys
            .saturating_mul(self.poly_bytes)
            .saturating_add(terms.saturating_mul(self.term_bytes))
    }

    /// Return an error if `bytes` more would pass the limit.
    pub(crate) fn check_bytes(&self, bytes: usize) -> Result<(), CertifyError> {
        match self.limits.memory {
            Some(limit) if self.held.saturating_add(bytes) > limit => Err(
                CertifyError::WriterExhausted(ComputeError::MemoryLimitExceeded),
            ),
            _ => Ok(()),
        }
    }

    /// Return an error if `polys` polynomials of `terms` terms together
    /// would pass the limit.
    pub(crate) fn check(&self, polys: usize, terms: usize) -> Result<(), CertifyError> {
        self.check_bytes(self.bytes(polys, terms))
    }

    /// Add `polys` polynomials of `terms` terms to the data the emitter
    /// holds, and return an error if the limit is passed.
    pub(crate) fn hold(&mut self, polys: usize, terms: usize) -> Result<(), CertifyError> {
        self.held = self.held.saturating_add(self.bytes(polys, terms));
        self.check_bytes(0)
    }

    /// Charge one list of polynomials to the budget.
    pub(crate) fn hold_polys(&mut self, polys: &[Polynomial]) -> Result<(), CertifyError> {
        let terms = polys.iter().map(|poly| poly.terms.len()).sum();
        let coefficients: usize = polys.iter().map(|poly| poly.coefficient_bytes()).sum();
        self.hold_bytes(self.bytes(polys.len(), terms).saturating_add(coefficients))
    }

    /// Charge a raw byte count to the budget, and return an error if the
    /// limit is passed.
    ///
    /// [`assemble`] calls this with the length of the certificate bytes it
    /// wrote, so those bytes count against the limit.
    pub(crate) fn hold_bytes(&mut self, bytes: usize) -> Result<(), CertifyError> {
        self.held = self.held.saturating_add(bytes);
        self.check_bytes(0)
    }

    /// The bytes the emitter holds.
    pub(crate) fn held(&self) -> usize {
        self.held
    }

    /// Take `bytes` off the data the emitter holds.
    ///
    /// The caller releases what it dropped. A release of more than the
    /// budget holds leaves the budget at zero.
    pub(crate) fn release_bytes(&mut self, bytes: usize) {
        self.held = self.held.saturating_sub(bytes);
    }
}

/// Write the certificate for one computation.
///
/// `input` is the input F, `basis` the claimed reduced basis G, and
/// `origins` one cofactor list per basis element over the input, so that
/// g_j = sum_i c_ji * f_i. The emitter generates the membership and S-pair
/// cofactors itself.
///
/// The certificate carries the basis in the order the contract requires:
/// strictly descending by leading monomial. The origin entries move with
/// it. The bytes are a function of the arguments alone. `budget` stops the
/// work early; it never changes the bytes.
///
/// The emitter does not check the origin identities. It reports a defect
/// only where it cannot write a certificate at all: a division that leaves
/// a remainder, a zero basis element, a count that does not match, or an
/// exponent vector of the wrong width. A count error and a width error name
/// the index the caller passed; a division defect names the index the
/// certificate carries, after the sort.
///
/// The contract requires a monic basis, and an S-polynomial depends on the
/// leading coefficients. The emitter writes the basis it is given; the
/// verifier rejects one that is not monic.
pub(crate) fn assemble(
    ring: &PolynomialRing,
    input: &[Polynomial],
    basis: Vec<Polynomial>,
    origins: Vec<Origin>,
    budget: &mut WriterBudget,
) -> Result<Vec<u8>, CertifyError> {
    let modulus = ring.modulus();
    let nvars = ring.nvars();
    validate_parts(input, &basis, &origins, nvars)?;
    hold_parts(input, &basis, &origins, budget)?;
    let (basis, origins) = sort_basis(basis, origins);
    let membership = membership_representations(input, &basis, ring.ops(), budget)?;
    let spairs = spair_representations(&basis, ring.ops(), budget)?;
    budget.check_stop()?;
    let certificate = json::write(&json::Parts {
        modulus,
        nvars,
        input,
        basis: &basis,
        origin: &origins,
        membership: &membership,
        spairs: &spairs,
    });
    budget.hold_bytes(certificate.len())?;
    Ok(certificate)
}

fn validate_parts(
    input: &[Polynomial],
    basis: &[Polynomial],
    origins: &[Origin],
    nvars: usize,
) -> Result<(), CertifyError> {
    if origins.len() != basis.len() {
        return Err(CertifyError::Emitter(EmitterFault::OriginCount {
            found: origins.len(),
            expected: basis.len(),
        }));
    }
    validate_input(input, nvars)?;
    validate_basis(basis, nvars)?;
    validate_origins(origins, input.len(), nvars)
}

fn validate_input(input: &[Polynomial], nvars: usize) -> Result<(), CertifyError> {
    for (index, poly) in input.iter().enumerate() {
        check_width(poly, Place::Input(index), nvars)?;
    }
    Ok(())
}

fn validate_basis(basis: &[Polynomial], nvars: usize) -> Result<(), CertifyError> {
    for (index, poly) in basis.iter().enumerate() {
        if poly.is_zero() {
            return Err(CertifyError::Emitter(EmitterFault::BasisElementZero {
                index,
            }));
        }
        check_width(poly, Place::Basis(index), nvars)?;
    }
    Ok(())
}

fn validate_origins(
    origins: &[Origin],
    input_len: usize,
    nvars: usize,
) -> Result<(), CertifyError> {
    for (basis_index, entry) in origins.iter().enumerate() {
        if entry.len() != input_len {
            return Err(CertifyError::Emitter(EmitterFault::OriginEntryCount {
                basis: basis_index,
                found: entry.len(),
                expected: input_len,
            }));
        }
        for (input_index, poly) in entry.iter().enumerate() {
            let at = Place::Origin {
                basis: basis_index,
                input: input_index,
            };
            check_width(poly, at, nvars)?;
        }
    }
    Ok(())
}

fn hold_parts(
    input: &[Polynomial],
    basis: &[Polynomial],
    origins: &[Origin],
    budget: &mut WriterBudget,
) -> Result<(), CertifyError> {
    budget.check_stop()?;
    budget.hold_polys(input)?;
    budget.hold_polys(basis)?;
    for entry in origins {
        budget.hold_polys(entry)?;
    }
    Ok(())
}

fn sort_basis(basis: Vec<Polynomial>, origins: Vec<Origin>) -> (Vec<Polynomial>, Vec<Origin>) {
    let mut sorted: Vec<(Polynomial, Origin)> = basis.into_iter().zip(origins).collect();
    sorted.sort_by(|a, b| b.0.lm().cmp(&a.0.lm()));
    sorted.into_iter().unzip()
}

/// Report a term whose exponent vector is not `nvars` wide.
///
/// Every polynomial in a certificate carries one exponent per variable.
/// The engine arithmetic requires the same width, so the emitter checks it
/// before it divides.
fn check_width(poly: &Polynomial, at: Place, nvars: usize) -> Result<(), CertifyError> {
    for term in &poly.terms {
        if term.mono.nvars() != nvars {
            return Err(CertifyError::Emitter(EmitterFault::ExponentCount {
                at,
                found: term.mono.nvars(),
                expected: nvars,
            }));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
