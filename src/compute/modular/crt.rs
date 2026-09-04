//! The Chinese remainder accumulator of the multimodular driver.

use std::collections::BTreeMap;
use std::mem::size_of;

use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};

use crate::compute::{ComputeLimits, RunError};
use crate::poly::{Monomial, Polynomial};
use crate::ring::field::Felt;
use crate::ring::rational::limb_bytes;

/// The residues of one basis element modulo the product of the folded
/// primes.
#[derive(Clone, Debug)]
struct Entry {
    /// The leading monomial, which the category fixes.
    lead: Monomial,
    /// One residue per monomial the folded runs carried, ascending under
    /// grevlex. A residue is in `[1, M)`; a monomial no run carried is
    /// absent and reads as zero.
    residues: BTreeMap<Monomial, BigInt>,
}

/// The Chinese remainder accumulator over one category.
///
/// The accumulator holds one entry per basis element of the category, in
/// basis order, keyed by the element's leading monomial. A monomial a run
/// does not carry contributes the residue 0, and a monomial the
/// accumulator has not seen starts at 0 for the whole previous modulus,
/// which is what the earlier runs said about it.
///
/// Runs fold in the order the prime sequence produced them, never in
/// completion order.
#[derive(Clone, Debug)]
pub(crate) struct Accumulator {
    /// The product of the folded primes.
    modulus: BigInt,
    entries: Vec<Entry>,
    folded: usize,
}

impl Accumulator {
    /// The empty accumulator over the category `leads`, which is the
    /// leading monomial of each basis element in basis order.
    pub(crate) fn new(leads: &[Monomial]) -> Self {
        Accumulator {
            modulus: BigInt::one(),
            entries: leads
                .iter()
                .map(|lead| Entry {
                    lead: lead.clone(),
                    residues: BTreeMap::new(),
                })
                .collect(),
            folded: 0,
        }
    }

    /// The product of the folded primes.
    pub(crate) fn modulus(&self) -> &BigInt {
        &self.modulus
    }

    /// The number of folded runs.
    pub(crate) fn folded(&self) -> usize {
        self.folded
    }

    /// The residues of the basis element at `index`, ascending under
    /// grevlex.
    pub(crate) fn residues(&self, index: usize) -> impl Iterator<Item = (&Monomial, &BigInt)> {
        self.entries[index].residues.iter()
    }

    /// The number of basis elements.
    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }

    /// Fold the run `basis` modulo `p` into the accumulator.
    ///
    /// `basis` belongs to the category the accumulator was built over, so
    /// its leading monomials are the ones the entries hold. The limits stop
    /// the fold at a term: one basis element carries as many terms as the
    /// whole basis has, and each one costs a big-integer step.
    pub(crate) fn fold(
        &mut self,
        basis: &[Polynomial],
        p: u64,
        limits: &ComputeLimits,
    ) -> Result<(), RunError> {
        debug_assert_eq!(
            basis.len(),
            self.entries.len(),
            "the category fixes the shape"
        );
        let modulus_residue = (&self.modulus % p)
            .to_u64()
            .expect("a residue below the modulus fits a u64");
        // The primes are distinct, so the product of the earlier ones is
        // invertible modulo this one.
        let inverse = Felt::from_residue(modulus_residue).inv(p);

        let mut term = 0usize;
        for (entry, element) in self.entries.iter_mut().zip(basis) {
            debug_assert_eq!(
                element.lm(),
                Some(&entry.lead),
                "the category fixes the leading monomial"
            );
            let mut fresh: BTreeMap<&Monomial, u64> = BTreeMap::new();
            for held in &element.terms {
                if let Some(stop) = limits.stop_every(term) {
                    return Err(stop);
                }
                term += 1;
                fresh.insert(&held.mono, held.coeff.value());
            }
            for (mono, residue) in &mut entry.residues {
                if let Some(stop) = limits.stop_every(term) {
                    return Err(stop);
                }
                term += 1;
                let value = fresh.remove(mono).unwrap_or(0);
                *residue = step(residue, value, p, &self.modulus, inverse);
            }
            let mut added: Vec<(Monomial, BigInt)> = Vec::with_capacity(fresh.len());
            for (mono, value) in fresh {
                if let Some(stop) = limits.stop_every(term) {
                    return Err(stop);
                }
                term += 1;
                added.push((
                    mono.clone(),
                    step(&BigInt::zero(), value, p, &self.modulus, inverse),
                ));
            }
            entry.residues.extend(added);
        }

        self.modulus *= p;
        self.folded += 1;
        Ok(())
    }

    /// The heap bytes the accumulator holds.
    pub(crate) fn heap_bytes(&self, nvars: usize) -> usize {
        let per_residue = size_of::<(Monomial, BigInt)>() + crate::poly::heap_exps_bytes(nvars);
        self.entries
            .iter()
            .fold(limb_bytes(&self.modulus), |bytes, entry| {
                entry.residues.values().fold(bytes, |bytes, residue| {
                    bytes
                        .saturating_add(per_residue)
                        .saturating_add(limb_bytes(residue))
                })
            })
    }
}

/// The residue modulo `modulus * p` of the value that is `residue` modulo
/// `modulus` and `value` modulo `p`.
///
/// `inverse` is the inverse of `modulus` modulo `p`, which the caller
/// computes once per fold.
fn step(residue: &BigInt, value: u64, p: u64, modulus: &BigInt, inverse: Felt) -> BigInt {
    let held = (residue % p)
        .to_u64()
        .expect("a residue below the modulus fits a u64");
    let difference = Felt::from_residue(value).sub(Felt::from_residue(held), p);
    let step = difference.mul(inverse, p).value();
    residue + modulus * step
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::PolynomialRing;
    use num_traits::Signed;

    /// The residue in `[0, M)` that is `residues[i]` modulo `primes[i]`,
    /// solved directly.
    fn direct(residues: &[u64], primes: &[u64]) -> BigInt {
        let modulus: BigInt = primes.iter().fold(BigInt::one(), |held, p| held * p);
        let mut value = BigInt::ZERO;
        for (&residue, &p) in residues.iter().zip(primes) {
            let rest = &modulus / p;
            let inverse = Felt::from_residue(
                (&rest % p)
                    .to_u64()
                    .expect("a residue below the modulus fits a u64"),
            )
            .inv(p);
            value += &rest * inverse.mul(Felt::from_residue(residue), p).value();
        }
        value % modulus
    }

    fn ring(p: u64) -> PolynomialRing {
        PolynomialRing::prime_field(p, ["x", "y"]).expect("the modulus is prime")
    }

    #[test]
    fn folding_agrees_with_a_direct_solve() {
        let primes = [10007u64, 10009, 10037];
        let coefficients = [7u64, 5000, 9001];
        let bases: Vec<Vec<Polynomial>> = primes
            .iter()
            .zip(coefficients)
            .map(|(&p, coeff)| {
                let ring = ring(p);
                vec![
                    ring.polynomial([(1i64, [2, 0]), (coeff as i64, [0, 1])])
                        .expect("the exponent vectors match the ring"),
                ]
            })
            .collect();
        let leads = vec![bases[0][0].lm().expect("the element is not zero").clone()];
        let mut accumulator = Accumulator::new(&leads);
        let limits = ComputeLimits::default();
        for (basis, &p) in bases.iter().zip(&primes) {
            accumulator
                .fold(basis, p, &limits)
                .expect("no limit stops it");
        }
        assert_eq!(accumulator.folded(), 3);
        let held: Vec<BigInt> = accumulator
            .residues(0)
            .map(|(_, residue)| residue.clone())
            .collect();
        // The monomials run ascending: y first, then x^2.
        assert_eq!(held[0], direct(&coefficients, &primes));
        assert_eq!(held[1], direct(&[1, 1, 1], &primes));
    }

    #[test]
    fn a_fold_checks_limits_at_its_first_term() {
        let p = 10007u64;
        let ring = ring(p);
        let element = ring
            .polynomial([(1, [1, 0])])
            .expect("the exponent vectors match the ring");
        let leads = vec![element.lm().expect("the element is not zero").clone()];
        let mut accumulator = Accumulator::new(&leads);
        let limits = ComputeLimits {
            deadline: Some(std::time::Instant::now()),
            ..ComputeLimits::default()
        };

        assert_eq!(
            accumulator.fold(&[element], p, &limits),
            Err(RunError::Compute(crate::ComputeError::Timeout))
        );
    }

    #[test]
    fn a_monomial_of_a_later_run_starts_at_zero() {
        let primes = [10007u64, 10009];
        let first = ring(primes[0]);
        let second = ring(primes[1]);
        let bases = [
            vec![
                first
                    .polynomial([(1i64, [2, 0])])
                    .expect("the exponent vectors match the ring"),
            ],
            vec![
                second
                    .polynomial([(1i64, [2, 0]), (3i64, [0, 1])])
                    .expect("the exponent vectors match the ring"),
            ],
        ];
        let leads = vec![bases[0][0].lm().expect("the element is not zero").clone()];
        let mut accumulator = Accumulator::new(&leads);
        let limits = ComputeLimits::default();
        for (basis, &p) in bases.iter().zip(&primes) {
            accumulator
                .fold(basis, p, &limits)
                .expect("no limit stops it");
        }
        let held: Vec<BigInt> = accumulator
            .residues(0)
            .map(|(_, residue)| residue.clone())
            .collect();
        // The first run said the coefficient of y is zero, so the value is
        // 0 modulo the first prime and 3 modulo the second.
        assert_eq!(held[0], direct(&[0, 3], &primes));
        assert!(held.iter().all(|residue| residue.is_positive()));
        assert!(accumulator.heap_bytes(2) > 0);
    }
}
