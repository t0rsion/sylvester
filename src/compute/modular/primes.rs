//! The prime sequence of the multimodular driver.

use crate::compute::{ComputeLimits, RunError};
use crate::ring::field::is_prime;

/// The first prime of the sequence, which is the largest modulus a ring
/// accepts.
pub(crate) const FIRST: u64 = (1 << 31) - 1;

/// One prime of the sequence, with what the driver skipped before it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Prime {
    /// The prime itself.
    pub(crate) value: u64,
    /// The primes the driver skipped before this one.
    pub(crate) skipped_before: usize,
}

/// The primes the driver runs, largest first.
///
/// The sequence walks down the odd numbers from 2^31 - 1 and keeps the
/// primes a predicate accepts. It is a function of its position and of
/// that predicate, so two runs over one input take the same primes in the
/// same order at any thread count. The prime 2 is not in it: each prime of
/// the range carries about 31 bits of the modulus and 2 carries one, so
/// including it would cost a whole run for one bit.
///
/// The sequence is finite. It ends below 3, which is
/// [`ComputeError::PrimesExhausted`] at the driver.
///
/// [`ComputeError::PrimesExhausted`]: crate::ComputeError::PrimesExhausted
#[derive(Clone, Debug)]
pub(crate) struct Sequence {
    /// The primes the predicate accepted, in sequence order. The driver
    /// consumes them by index and can leave an index unconsumed, so the
    /// sequence holds what it has already found.
    accepted: Vec<Prime>,
    /// The next odd number the search tests, or 1 once the search is over.
    next: u64,
    /// The primes the predicate rejected so far.
    skipped: usize,
}

impl Sequence {
    /// The sequence from its start.
    pub(crate) fn new() -> Self {
        Sequence {
            accepted: Vec::new(),
            next: FIRST,
            skipped: 0,
        }
    }

    /// The prime at `index`, or `None` once the sequence ends.
    ///
    /// `accept` reports whether the driver runs a prime. It is called once
    /// per prime of the raw sequence, in order, and the result is held, so
    /// a caller that asks for one index twice gets the same value and
    /// costs one search.
    ///
    /// The search walks about 2^30 odd numbers before the sequence ends,
    /// so `limits` stops it.
    pub(crate) fn at(
        &mut self,
        index: usize,
        mut accept: impl FnMut(u64) -> bool,
        limits: &ComputeLimits,
    ) -> Result<Option<Prime>, RunError> {
        while self.accepted.len() <= index {
            let Some(candidate) = self.search(limits)? else {
                return Ok(None);
            };
            if accept(candidate) {
                self.accepted.push(Prime {
                    value: candidate,
                    skipped_before: self.skipped,
                });
            } else {
                self.skipped += 1;
            }
        }
        Ok(Some(self.accepted[index]))
    }

    /// The bytes the held primes cost.
    pub(crate) fn heap_bytes(&self) -> usize {
        self.accepted
            .len()
            .saturating_mul(std::mem::size_of::<Prime>())
    }

    /// The next prime below the last one the search returned.
    fn search(&mut self, limits: &ComputeLimits) -> Result<Option<u64>, RunError> {
        let mut step = 0usize;
        while self.next >= 3 {
            if let Some(stop) = limits.stop_every(step) {
                return Err(stop);
            }
            step += 1;
            let candidate = self.next;
            self.next -= 2;
            if is_prime(candidate) {
                return Ok(Some(candidate));
            }
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Budget;
    use std::time::Duration;

    fn all(_: u64) -> bool {
        true
    }

    fn open() -> ComputeLimits {
        ComputeLimits::default()
    }

    fn prime_at(sequence: &mut Sequence, index: usize) -> Option<Prime> {
        sequence
            .at(index, all, &open())
            .expect("no limit stops the search")
    }

    #[test]
    fn the_sequence_starts_at_the_largest_modulus_and_descends() {
        let mut sequence = Sequence::new();
        let first: Vec<u64> = (0..4)
            .map(|index| {
                prime_at(&mut sequence, index)
                    .expect("the sequence holds four primes")
                    .value
            })
            .collect();
        assert_eq!(first[0], FIRST);
        assert!(first.windows(2).all(|pair| pair[0] > pair[1]));
        assert!(first.iter().all(|&p| is_prime(p)));
        assert!(sequence.heap_bytes() > 0);
    }

    #[test]
    fn one_index_always_reads_one_prime() {
        let mut sequence = Sequence::new();
        let third = prime_at(&mut sequence, 2).expect("the sequence holds three primes");
        assert_eq!(prime_at(&mut sequence, 2), Some(third));
        assert_eq!(
            prime_at(&mut sequence, 0).map(|prime| prime.value),
            Some(FIRST)
        );
    }

    #[test]
    fn a_rejected_prime_counts_and_leaves_the_sequence() {
        let mut sequence = Sequence::new();
        let first = sequence
            .at(0, |value| value != FIRST, &open())
            .expect("no limit stops the search")
            .expect("a prime is left");
        assert_ne!(first.value, FIRST);
        assert_eq!(first.skipped_before, 1);
        assert_eq!(
            sequence
                .at(1, |value| value != FIRST, &open())
                .expect("no limit stops the search")
                .map(|prime| prime.skipped_before),
            Some(1)
        );
    }

    /// The search reads the deadline, so an exhausted one stops it.
    #[test]
    fn the_search_stops_on_an_exhausted_deadline() {
        let mut sequence = Sequence::new();
        let limits = ComputeLimits::of_budget(&Budget::new().timeout(Duration::ZERO));
        assert_eq!(
            sequence.at(0, all, &limits).map(|prime| prime.is_some()),
            Err(RunError::Compute(crate::ComputeError::Timeout))
        );
    }
}
