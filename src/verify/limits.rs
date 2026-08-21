//! Resource caps for one verification.

use std::time::Instant;

use super::algebra::{Exp, Term};
use super::error::{Cap, VerifyError};

/// The number of steps between two deadline polls inside a loop.
///
/// Reading the clock costs more than one step of decoding or one term of
/// arithmetic, so a loop polls once every stride. A deadline can overrun by
/// one stride of work.
pub(crate) const WORK_STRIDE: usize = 1024;

/// Caps the verifier enforces while it reads and while it checks.
///
/// The verifier applies a cap before it allocates the item the cap covers,
/// so a certificate that breaches a cap never allocates past it. Lower the
/// caps for bytes from an untrusted peer.
///
/// The caps bound memory. A deadline bounds time, and every loop over the
/// certificate polls it on a fixed stride of work, so a run can overrun its
/// deadline by one stride.
#[derive(Clone, Debug)]
pub struct Limits {
    /// The largest certificate the verifier reads, in bytes.
    pub max_bytes: usize,
    /// The largest number of polynomials in one certificate. The count
    /// covers the input, the basis, and every cofactor.
    pub max_polys: usize,
    /// The largest number of terms in one polynomial.
    pub max_terms_per_poly: usize,
    /// The largest number of terms in one certificate.
    pub max_total_terms: usize,
    /// The largest number of entries in the `origin`, `membership`, and
    /// `spairs` arrays together. One entry holds one list of cofactors.
    /// An entry with no cofactor still costs memory, so this cap covers it.
    pub max_entries: usize,
    /// The largest number of bytes the verifier arithmetic holds at once.
    ///
    /// The verifier multiplies and adds polynomials to check the
    /// identities. Those buffers are not part of the certificate, so no
    /// other cap covers them. One prospective term costs
    /// `size_of::<Term>() + nvars * size_of::<Exp>()`, 48 bytes over two
    /// variables and 1064 bytes over the contract maximum of 256, so a term
    /// count bounds memory only once the variable count is fixed.
    ///
    /// The cap covers the buffers live at once, not one result alone. It
    /// counts term and exponent payloads, not allocator metadata, so it is
    /// not a cap on resident memory.
    pub max_intermediate_bytes: usize,
    /// The instant the work must stop at.
    pub deadline: Option<Instant>,
}

impl Default for Limits {
    fn default() -> Self {
        Limits {
            max_bytes: 64 << 20,
            max_polys: 4_000_000,
            max_terms_per_poly: 1 << 20,
            max_total_terms: 8_000_000,
            max_entries: 4_000_000,
            max_intermediate_bytes: 512 << 20,
            deadline: None,
        }
    }
}

/// The memory one term over `nvars` variables takes.
///
/// A term owns its exponent vector, so the cost grows with the variable
/// count. The value leaves out allocator metadata, which no cap models.
pub(crate) fn term_bytes(nvars: usize) -> usize {
    size_of::<Term>() + nvars * size_of::<Exp>()
}

impl Limits {
    /// Return an error if the deadline has passed.
    pub fn check_deadline(&self) -> Result<(), VerifyError> {
        self.ticker().now()
    }

    /// The budget the arithmetic runs under, over `nvars` variables.
    ///
    /// Every polynomial the arithmetic touches carries one exponent per
    /// variable of the certificate, checked before any identity check runs,
    /// so one `nvars` prices every term.
    pub(crate) fn budget(&self, nvars: usize) -> Budget {
        Budget {
            max_bytes: self.max_intermediate_bytes,
            term_bytes: term_bytes(nvars),
            ticker: self.ticker(),
        }
    }

    /// The ticker the decoder runs under.
    pub(crate) fn ticker(&self) -> Ticker {
        Ticker {
            deadline: self.deadline,
            left: WORK_STRIDE,
        }
    }
}

/// A deadline poll that reads the clock once every [`WORK_STRIDE`] steps.
#[derive(Clone, Debug)]
pub(crate) struct Ticker {
    deadline: Option<Instant>,
    left: usize,
}

impl Ticker {
    /// Count one step of work, and poll the deadline every stride.
    pub(crate) fn step(&mut self) -> Result<(), VerifyError> {
        self.left -= 1;
        if self.left > 0 {
            return Ok(());
        }
        self.left = WORK_STRIDE;
        self.now()
    }

    /// Poll the deadline now.
    pub(crate) fn now(&self) -> Result<(), VerifyError> {
        match self.deadline {
            Some(deadline) if Instant::now() >= deadline => Err(VerifyError::DeadlineExceeded),
            _ => Ok(()),
        }
    }
}

/// The byte budget and the deadline one obligation runs under.
///
/// The budget covers every buffer the arithmetic reserves and every buffer
/// the caller holds live beside it, charged for its capacity until it is
/// dropped. It excludes the decoded certificate, which the term caps cover.
/// The verifier charges before it reserves, so a check that needs more than
/// the cap fails before it allocates.
#[derive(Clone, Debug)]
pub(crate) struct Budget {
    max_bytes: usize,
    term_bytes: usize,
    ticker: Ticker,
}

impl Budget {
    fn exceeded<T>(&self) -> Result<T, VerifyError> {
        Err(VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit: self.max_bytes,
        })
    }

    /// Count one step of work, and poll the deadline every stride.
    pub(crate) fn step(&mut self) -> Result<(), VerifyError> {
        self.ticker.step()
    }

    /// Poll the deadline now.
    pub(crate) fn check_deadline(&self) -> Result<(), VerifyError> {
        self.ticker.now()
    }

    /// Return an error if buffers of these term counts, held at once, cost
    /// more than the budget.
    ///
    /// A size that overflows a `usize` is above every budget.
    pub(crate) fn check_live(&self, counts: &[usize]) -> Result<(), VerifyError> {
        let mut terms = 0usize;
        for &count in counts {
            terms = match terms.checked_add(count) {
                Some(sum) => sum,
                None => return self.exceeded(),
            };
        }
        match terms.checked_mul(self.term_bytes) {
            Some(bytes) if bytes <= self.max_bytes => Ok(()),
            _ => self.exceeded(),
        }
    }

    /// Return the term count of the raw product of two term lists.
    ///
    /// The count is the capacity the product buffer needs before the merge
    /// drops equal monomials. `live` is the term count the caller holds
    /// beside that buffer. Return an error if the two together are above
    /// the budget.
    pub(crate) fn check_product(
        &self,
        left: usize,
        right: usize,
        live: usize,
    ) -> Result<usize, VerifyError> {
        let Some(count) = left.checked_mul(right) else {
            return self.exceeded();
        };
        self.check_live(&[count, live])?;
        Ok(count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    /// A budget that admits `terms` terms over two variables, and no more.
    fn budget_for(terms: usize) -> Budget {
        Limits {
            max_intermediate_bytes: terms * term_bytes(2),
            ..Limits::default()
        }
        .budget(2)
    }

    #[test]
    fn a_passed_deadline_stops_the_work() {
        let limits = Limits {
            deadline: Some(Instant::now() - Duration::from_secs(1)),
            ..Limits::default()
        };
        assert_eq!(limits.check_deadline(), Err(VerifyError::DeadlineExceeded));
    }

    #[test]
    fn no_deadline_never_stops_the_work() {
        assert_eq!(Limits::default().check_deadline(), Ok(()));
    }

    #[test]
    fn a_ticker_polls_the_deadline_once_every_stride() {
        let mut ticker = Limits {
            deadline: Some(Instant::now() - Duration::from_secs(1)),
            ..Limits::default()
        }
        .ticker();
        for _ in 1..WORK_STRIDE {
            assert_eq!(ticker.step(), Ok(()));
        }
        assert_eq!(ticker.step(), Err(VerifyError::DeadlineExceeded));
    }

    #[test]
    fn the_budget_admits_buffers_up_to_its_cap() {
        let budget = budget_for(4);
        assert_eq!(budget.check_live(&[4]), Ok(()));
        assert_eq!(budget.check_live(&[2, 2]), Ok(()));
        assert_eq!(budget.check_product(2, 2, 0), Ok(4));
    }

    #[test]
    fn the_budget_rejects_buffers_above_its_cap() {
        let budget = budget_for(4);
        let exceeded = VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit: 4 * term_bytes(2),
        };
        assert_eq!(budget.check_live(&[5]), Err(exceeded.clone()));
        assert_eq!(budget.check_live(&[2, 2, 1]), Err(exceeded.clone()));
        assert_eq!(budget.check_product(2, 2, 1), Err(exceeded));
    }

    #[test]
    fn a_wider_ring_prices_the_same_term_count_higher() {
        let limits = Limits {
            max_intermediate_bytes: 100 * term_bytes(2),
            ..Limits::default()
        };
        assert_eq!(limits.budget(2).check_live(&[100]), Ok(()));
        assert_eq!(
            limits.budget(256).check_live(&[100]),
            Err(VerifyError::CapExceeded {
                cap: Cap::IntermediateBytes,
                limit: 100 * term_bytes(2),
            })
        );
    }

    #[test]
    fn a_size_that_overflows_is_above_every_budget() {
        let budget = Limits::default().budget(2);
        let exceeded = VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit: Limits::default().max_intermediate_bytes,
        };
        assert_eq!(
            budget.check_product(usize::MAX, 2, 0),
            Err(exceeded.clone())
        );
        assert_eq!(budget.check_live(&[usize::MAX, 2]), Err(exceeded.clone()));
        assert_eq!(budget.check_live(&[usize::MAX]), Err(exceeded));
    }
}
