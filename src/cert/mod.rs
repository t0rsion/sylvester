//! Certificate emission for the `sylv-gb-cert-v1` contract.
//!
//! This module is untrusted. It reads engine values and writes the
//! canonical bytes of the `sylv-gb-cert-v1` contract. The verifier in
//! [`crate::verify`] is the trust boundary: it reads those bytes back and
//! accepts or rejects them. A wrong cofactor here becomes a rejection
//! there, never a silent repair. This module never imports `crate::verify`,
//! in its code or in its tests:
//! [`crate::compute::groebner_basis_certified`] runs the emitter here, then
//! hands the bytes to the verifier itself. `tests/certified.rs` runs the
//! two sides against each other and enforces the import boundary.
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

use std::time::Instant;

use crate::certificate::{CertifyError, EmitterFault, Place};
use crate::compute::ComputeError;
use crate::poly::{Polynomial, Term};
use crate::ring::PolynomialRing;

use origin::Origin;
use represent::{membership_representations, spair_representations};

/// Return an error once the deadline has passed.
pub(crate) fn check_deadline(deadline: Option<Instant>) -> Result<(), ComputeError> {
    match deadline {
        Some(deadline) if Instant::now() >= deadline => Err(ComputeError::Timeout),
        _ => Ok(()),
    }
}

/// The deadline and the memory limit one emission runs under.
///
/// The emitter charges the data it keeps before it builds more. The count
/// covers the input, the basis, the origins, the cofactors the emitter
/// generates, and the certificate bytes. The count estimates the live
/// emitter data. It is not a count of process memory. The limit is the one
/// the engine ran under.
///
/// One term costs its own size plus one exponent per variable. One
/// polynomial costs the header as well, because a basis of n elements
/// carries n cofactors per S-pair and most of them are zero.
pub(crate) struct Budget {
    deadline: Option<Instant>,
    limit: Option<usize>,
    poly_bytes: usize,
    term_bytes: usize,
    held: usize,
}

impl Budget {
    pub(crate) fn new(deadline: Option<Instant>, limit: Option<usize>, nvars: usize) -> Self {
        Budget {
            deadline,
            limit,
            poly_bytes: size_of::<Polynomial>(),
            term_bytes: size_of::<Term>() + nvars * size_of::<u16>(),
            held: 0,
        }
    }

    /// A budget that stops nothing.
    #[cfg(test)]
    pub(crate) fn unlimited(nvars: usize) -> Self {
        Budget::new(None, None, nvars)
    }

    pub(crate) fn check_deadline(&self) -> Result<(), ComputeError> {
        check_deadline(self.deadline)
    }

    /// The bytes of `polys` polynomials that hold `terms` terms together.
    fn bytes(&self, polys: usize, terms: usize) -> usize {
        polys
            .saturating_mul(self.poly_bytes)
            .saturating_add(terms.saturating_mul(self.term_bytes))
    }

    /// Return an error if `bytes` more would pass the limit.
    pub(crate) fn check_bytes(&self, bytes: usize) -> Result<(), ComputeError> {
        match self.limit {
            Some(limit) if self.held.saturating_add(bytes) > limit => {
                Err(ComputeError::MemoryLimitExceeded)
            }
            _ => Ok(()),
        }
    }

    /// Return an error if `polys` polynomials of `terms` terms together
    /// would pass the limit.
    pub(crate) fn check(&self, polys: usize, terms: usize) -> Result<(), ComputeError> {
        self.check_bytes(self.bytes(polys, terms))
    }

    /// Add `polys` polynomials of `terms` terms to the data the emitter
    /// holds, and return an error if the limit is passed.
    fn hold(&mut self, polys: usize, terms: usize) -> Result<(), ComputeError> {
        self.held = self.held.saturating_add(self.bytes(polys, terms));
        self.check_bytes(0)
    }

    /// Charge one list of polynomials to the budget.
    pub(crate) fn hold_polys(&mut self, polys: &[Polynomial]) -> Result<(), ComputeError> {
        let terms = polys.iter().map(|poly| poly.terms.len()).sum();
        self.hold(polys.len(), terms)
    }

    /// Charge a raw byte count to the budget, and return an error if the
    /// limit is passed.
    ///
    /// [`assemble`] calls this with the length of the certificate bytes it
    /// wrote, so the certificate itself, not just the data it came from,
    /// counts against the limit.
    pub(crate) fn hold_bytes(&mut self, bytes: usize) -> Result<(), ComputeError> {
        self.held = self.held.saturating_add(bytes);
        self.check_bytes(0)
    }

    /// The bytes the emitter holds.
    pub(crate) fn held(&self) -> usize {
        self.held
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
    budget: &mut Budget,
) -> Result<Vec<u8>, CertifyError> {
    let modulus = ring.modulus();
    let nvars = ring.nvars();

    if origins.len() != basis.len() {
        return Err(CertifyError::Emitter(EmitterFault::OriginCount {
            found: origins.len(),
            expected: basis.len(),
        }));
    }

    for (index, poly) in input.iter().enumerate() {
        check_width(poly, Place::Input(index), nvars)?;
    }
    for (index, poly) in basis.iter().enumerate() {
        if poly.is_zero() {
            return Err(CertifyError::Emitter(EmitterFault::BasisElementZero {
                index,
            }));
        }
        check_width(poly, Place::Basis(index), nvars)?;
    }
    for (basis_index, entry) in origins.iter().enumerate() {
        if entry.len() != input.len() {
            return Err(CertifyError::Emitter(EmitterFault::OriginEntryCount {
                basis: basis_index,
                found: entry.len(),
                expected: input.len(),
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

    budget.check_deadline()?;
    budget.hold_polys(input)?;
    budget.hold_polys(&basis)?;
    for entry in &origins {
        budget.hold_polys(entry)?;
    }

    // The basis and the origins move into the sort together, so the sort
    // copies no cofactor list.
    let mut sorted: Vec<(Polynomial, Origin)> = basis.into_iter().zip(origins).collect();
    sorted.sort_by(|a, b| b.0.lm().cmp(&a.0.lm()));
    let (sorted_basis, sorted_origins): (Vec<Polynomial>, Vec<Origin>) = sorted.into_iter().unzip();

    let membership = membership_representations(input, &sorted_basis, modulus, budget)?;
    let spairs = spair_representations(&sorted_basis, modulus, budget)?;

    budget.check_deadline()?;
    // The writer charges the growing buffer per polynomial, so a
    // certificate the limit forbids stops the writer instead of being
    // built in full and rejected after the fact.
    let certificate = json::write(
        &json::Parts {
            modulus,
            nvars,
            input,
            basis: &sorted_basis,
            origin: &sorted_origins,
            membership: &membership,
            spairs: &spairs,
        },
        budget,
    )?;
    // Charge the certificate its own byte length, not the estimate already
    // held for the same data: that would double count. Charging nothing
    // would let a limit too small for even an empty certificate pass.
    budget.hold_bytes(certificate.len())?;
    Ok(certificate)
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
