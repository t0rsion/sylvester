//! Critical pairs and the Gebauer-Moller update (design section 3.4).
//!
//! [`PairSet`] holds every queued pair in one flat vector.
//! [`PairSet::update`] runs when a basis element arrives, and
//! [`PairSet::take_lowest_degree`] hands the run loop one batch.
//!
//! This module names basis elements by index and reads their leading
//! monomials from a slice, so it does not depend on the basis
//! representation.

use super::monomial::{Lanes, MonomialId, MonomialTable, intern_from, lcm_into};
use super::{Deadline, F4Error};

/// The slot of a basis element that is not a candidate of this update.
const NO_SLOT: u32 = u32::MAX;

/// One critical pair.
///
/// `lcm` is `lcm(lm_i, lm_j)` in the basis table, and `deg` is its total
/// degree. Indices are ordered `i < j`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Pair {
    /// The lcm of the two leading monomials, in the basis table.
    pub(crate) lcm: MonomialId,
    /// The total degree of `lcm`.
    pub(crate) deg: u32,
    /// The smaller basis index.
    pub(crate) i: u32,
    /// The larger basis index.
    pub(crate) j: u32,
}

/// How many pairs each step of the update made or dropped.
///
/// `generated` counts every candidate the update forms, one per live
/// element. The four discard counters count candidates the product
/// criterion, criterion M, and criterion F drop, and queued pairs
/// criterion B deletes. So the number of pairs the update queues is
/// `generated - product - criterion_m - criterion_f`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct PairCounters {
    /// Candidate pairs formed.
    pub(crate) generated: u64,
    /// Candidates dropped because the two leading monomials are coprime.
    pub(crate) product: u64,
    /// Queued pairs deleted by criterion B.
    pub(crate) criterion_b: u64,
    /// Candidates dropped by criterion M.
    pub(crate) criterion_m: u64,
    /// Candidates dropped by criterion F.
    pub(crate) criterion_f: u64,
}

/// The limits [`PairSet::take_lowest_degree`] reads.
///
/// `max_batch_pairs` caps one batch, so one degree class cannot build a
/// matrix past the memory budget. The default is no cap. A cap of 0 acts
/// as 1, because an empty batch would stall the run loop.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct SelectOptions {
    /// The largest number of pairs one batch takes, or `None` for no cap.
    pub(crate) max_batch_pairs: Option<usize>,
}

/// One candidate pair `(i, t)` of one update.
#[derive(Clone, Copy, Debug)]
struct Candidate {
    /// `lcm(lm_i, lm_t)`, in the update workspace table.
    lcm: MonomialId,
    /// The other basis index.
    i: u32,
    /// Whether the two leading monomials share no variable.
    coprime: bool,
    /// Whether a criterion has dropped this candidate.
    dropped: bool,
}

/// The queued critical pairs of one run.
pub(crate) struct PairSet {
    pairs: Vec<Pair>,
    counters: PairCounters,
    candidates: Vec<Candidate>,
    /// The candidate slot of each basis index, or [`NO_SLOT`]. Reused
    /// across updates, cleared at the end of each one.
    slot_of: Vec<u32>,
}

impl Default for PairSet {
    fn default() -> Self {
        Self::new()
    }
}

impl PairSet {
    /// An empty pair set.
    pub(crate) fn new() -> Self {
        PairSet {
            pairs: Vec::new(),
            counters: PairCounters::default(),
            candidates: Vec::new(),
            slot_of: Vec::new(),
        }
    }

    /// The number of queued pairs.
    pub(crate) fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Report whether the queue is empty.
    pub(crate) fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// The queued pairs, in queue order.
    pub(crate) fn pairs(&self) -> &[Pair] {
        &self.pairs
    }

    /// The counters of every update so far.
    pub(crate) fn counters(&self) -> PairCounters {
        self.counters
    }

    /// Drop every queued pair, keeping the counters.
    ///
    /// The run loop calls this when a basis element is a nonzero constant,
    /// because the basis is then the unit ideal and no pair is left to
    /// process.
    pub(crate) fn clear(&mut self) {
        self.pairs.clear();
    }

    /// Put a batch's deferred pairs back in the queue.
    ///
    /// The run loop calls this when the projected batch passes the memory
    /// budget and it retries with half the pairs (design section 3.9). The
    /// returned pairs keep their lcm ids, which belong to the basis table
    /// and stay valid.
    pub(crate) fn requeue(&mut self, deferred: Vec<Pair>) -> Result<(), F4Error> {
        self.pairs
            .try_reserve(deferred.len())
            .map_err(|_| F4Error::MemoryLimitExceeded)?;
        self.pairs.extend(deferred);
        Ok(())
    }

    /// The bytes the pair set holds, counting capacities.
    pub(crate) fn memory_bytes(&self) -> usize {
        self.pairs
            .capacity()
            .saturating_mul(size_of::<Pair>())
            .saturating_add(
                self.candidates
                    .capacity()
                    .saturating_mul(size_of::<Candidate>()),
            )
            .saturating_add(self.slot_of.capacity().saturating_mul(size_of::<u32>()))
    }

    /// Run the Gebauer-Moller update for the new element `t`.
    ///
    /// `leads` holds the leading monomial of every element up to and
    /// including `t`, in `table`. `live` holds the live element indices as
    /// they were before `t` arrived, ascending, and none of them is `t`.
    /// Retirement runs after this call, so an element this update names is
    /// still live here (design section 3.3).
    ///
    /// Candidate lcms go into `ws`, which this call clears when it
    /// finishes. The lcm of a surviving pair moves into `table`.
    pub(crate) fn update<L: Lanes>(
        &mut self,
        table: &mut MonomialTable<L>,
        ws: &mut MonomialTable<L>,
        leads: &[MonomialId],
        live: &[u32],
        t: u32,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        self.generate(table, ws, leads, live, t, clock)?;
        self.criterion_b(table, ws, leads, t, clock)?;
        self.sort_candidates(ws);
        self.criterion_m(ws, clock)?;
        self.criterion_f();
        self.queue_survivors(table, ws, t)?;

        // Step 3 dropped the coprime candidates, so the live set, not the
        // candidate array, names every slot step 1 took.
        for &i in live {
            self.slot_of[i as usize] = NO_SLOT;
        }
        self.candidates.clear();
        ws.clear();
        Ok(())
    }

    /// Step 1, the candidate array.
    fn generate<L: Lanes>(
        &mut self,
        table: &MonomialTable<L>,
        ws: &mut MonomialTable<L>,
        leads: &[MonomialId],
        live: &[u32],
        t: u32,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        debug_assert!((t as usize) < leads.len(), "element t must have a lead");
        if self.slot_of.len() < leads.len() {
            self.slot_of.resize(leads.len(), NO_SLOT);
        }
        let lm_t = leads[t as usize];
        for &i in live {
            clock.tick()?;
            debug_assert_ne!(i, t, "the live set excludes the new element");
            let lm_i = leads[i as usize];
            let lcm = lcm_into(ws, table, lm_i, lm_t)?;
            let coprime = table.is_coprime(lm_i, lm_t);
            self.slot_of[i as usize] = self.candidates.len() as u32;
            self.candidates.push(Candidate {
                lcm,
                i,
                coprime,
                dropped: coprime,
            });
            self.counters.generated += 1;
            self.counters.product += u64::from(coprime);
        }
        Ok(())
    }

    /// Step 2, criterion B over the old queue.
    ///
    /// The third condition is the coverage guard. B deletes `(i, j)`
    /// because `(i, t)` and `(j, t)` cover it, so both must be candidates
    /// of this update. An endpoint retired before `t` arrived produces no
    /// candidate, and without the guard B would name a pair that was never
    /// formed.
    fn criterion_b<L: Lanes>(
        &mut self,
        table: &mut MonomialTable<L>,
        ws: &mut MonomialTable<L>,
        leads: &[MonomialId],
        t: u32,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        let lm_t = leads[t as usize];
        let mut kept = 0;
        for index in 0..self.pairs.len() {
            clock.tick()?;
            let pair = self.pairs[index];
            if !self.deleted_by_b(table, ws, lm_t, pair)? {
                self.pairs[kept] = pair;
                kept += 1;
            }
        }
        self.counters.criterion_b += (self.pairs.len() - kept) as u64;
        self.pairs.truncate(kept);
        Ok(())
    }

    /// Report whether criterion B deletes the queued pair `pair`.
    fn deleted_by_b<L: Lanes>(
        &self,
        table: &mut MonomialTable<L>,
        ws: &mut MonomialTable<L>,
        lm_t: MonomialId,
        pair: Pair,
    ) -> Result<bool, F4Error> {
        if !table.divides(lm_t, pair.lcm) {
            return Ok(false);
        }
        let (si, sj) = (self.slot_of[pair.i as usize], self.slot_of[pair.j as usize]);
        if si == NO_SLOT || sj == NO_SLOT {
            return Ok(false);
        }
        // Both candidate lcms and the queued lcm must be comparable, and a
        // candidate lcm lives in ws. Moving the queued lcm into ws keeps the
        // three in one table, so the divisibility test uses that table's
        // masks.
        let target = intern_from(ws, table, pair.lcm)?;
        let (li, lj) = (
            self.candidates[si as usize].lcm,
            self.candidates[sj as usize].lcm,
        );
        Ok(proper_divisor(ws, li, target) && proper_divisor(ws, lj, target))
    }

    /// Step 3, drop the coprime candidates and sort what remains.
    fn sort_candidates<L: Lanes>(&mut self, ws: &MonomialTable<L>) {
        self.candidates.retain(|candidate| !candidate.coprime);
        self.candidates
            .sort_unstable_by(|a, b| ws.cmp(a.lcm, b.lcm).then(a.i.cmp(&b.i)));
    }

    /// Step 4, criterion M.
    ///
    /// The witness set is the candidate array as step 3 left it, not the
    /// array as M shrinks it. Strict divisibility is a strict partial
    /// order, so a dropped candidate never justifies dropping its own
    /// witness.
    ///
    /// The loop is quadratic in the candidate count, so it counts every
    /// divisibility test against the clock, not every candidate.
    fn criterion_m<L: Lanes>(
        &mut self,
        ws: &MonomialTable<L>,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        for index in 0..self.candidates.len() {
            let lcm = self.candidates[index].lcm;
            let mut covered = false;
            for other in 0..self.candidates.len() {
                clock.tick()?;
                if other != index && proper_divisor(ws, self.candidates[other].lcm, lcm) {
                    covered = true;
                    break;
                }
            }
            if covered {
                self.candidates[index].dropped = true;
                self.counters.criterion_m += 1;
            }
        }
        Ok(())
    }

    /// Step 5, criterion F.
    ///
    /// The candidates are sorted by lcm and then by index, so candidates
    /// with equal lcm are one run and the first survivor of the run has the
    /// smallest index.
    fn criterion_f(&mut self) {
        let mut kept: Option<MonomialId> = None;
        for candidate in &mut self.candidates {
            if candidate.dropped {
                continue;
            }
            if kept == Some(candidate.lcm) {
                candidate.dropped = true;
                self.counters.criterion_f += 1;
            } else {
                kept = Some(candidate.lcm);
            }
        }
    }

    /// Step 6, queue the survivors.
    fn queue_survivors<L: Lanes>(
        &mut self,
        table: &mut MonomialTable<L>,
        ws: &MonomialTable<L>,
        t: u32,
    ) -> Result<(), F4Error> {
        for index in 0..self.candidates.len() {
            let candidate = self.candidates[index];
            if candidate.dropped {
                continue;
            }
            debug_assert!(candidate.i < t, "t is the newest element");
            let lcm = intern_from(table, ws, candidate.lcm)?;
            self.pairs.push(Pair {
                lcm,
                deg: table.degree(lcm),
                i: candidate.i,
                j: t,
            });
        }
        Ok(())
    }

    /// Take every pair of the lowest degree present, up to the cap.
    ///
    /// The batch is sorted by lcm and then by the two indices, so the
    /// selection is a deterministic function of the queue. When the cap
    /// splits the degree class, the pairs past the cap stay in the queue
    /// and the next call takes them at the same degree.
    pub(crate) fn take_lowest_degree<L: Lanes>(
        &mut self,
        table: &MonomialTable<L>,
        options: &SelectOptions,
    ) -> Option<Vec<Pair>> {
        let degree = self.pairs.iter().map(|pair| pair.deg).min()?;
        // Partitioning to the end and splitting there moves only the batch.
        // Partitioning to the front would move the whole remainder.
        let mut split = self.pairs.len();
        for index in (0..self.pairs.len()).rev() {
            if self.pairs[index].deg == degree {
                split -= 1;
                self.pairs.swap(index, split);
            }
        }
        let mut batch = self.pairs.split_off(split);
        batch.sort_unstable_by(|a, b| {
            table
                .cmp(a.lcm, b.lcm)
                .then(a.i.cmp(&b.i))
                .then(a.j.cmp(&b.j))
        });
        if let Some(cap) = options.max_batch_pairs
            && batch.len() > cap.max(1)
        {
            self.pairs.extend(batch.split_off(cap.max(1)));
        }
        Some(batch)
    }
}

/// Report whether `a` divides `b` and differs from it.
fn proper_divisor<L: Lanes>(table: &MonomialTable<L>, a: MonomialId, b: MonomialId) -> bool {
    a != b && table.divides(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compute::f4::monomial::Lanes8;

    /// A run's two tables and the leading monomials of its elements.
    struct Fixture {
        table: MonomialTable<Lanes8>,
        ws: MonomialTable<Lanes8>,
        leads: Vec<MonomialId>,
    }

    impl Fixture {
        fn new(nvars: usize) -> Self {
            Fixture {
                table: MonomialTable::new(nvars),
                ws: MonomialTable::new(nvars),
                leads: Vec::new(),
            }
        }

        /// Add an element with the leading monomial `exps` and return its
        /// index.
        fn push(&mut self, exps: &[u32]) -> u32 {
            let lead = self.table.intern(exps).unwrap();
            self.leads.push(lead);
            (self.leads.len() - 1) as u32
        }

        fn update(&mut self, pairs: &mut PairSet, live: &[u32], t: u32) {
            pairs
                .update(
                    &mut self.table,
                    &mut self.ws,
                    &self.leads,
                    live,
                    t,
                    &mut Deadline::none(),
                )
                .unwrap();
        }

        /// The queued pairs as index pairs, sorted.
        fn queued(&self, pairs: &PairSet) -> Vec<(u32, u32)> {
            let mut out: Vec<(u32, u32)> = pairs.pairs().iter().map(|p| (p.i, p.j)).collect();
            out.sort_unstable();
            out
        }
    }

    #[test]
    fn the_product_criterion_drops_coprime_candidates() {
        let mut fixture = Fixture::new(2);
        fixture.push(&[2, 0]);
        fixture.push(&[0, 3]);
        let mut pairs = PairSet::new();
        fixture.update(&mut pairs, &[0], 1);

        assert!(pairs.is_empty());
        let counters = pairs.counters();
        assert_eq!(counters.generated, 1);
        assert_eq!(counters.product, 1);
    }

    #[test]
    fn criterion_m_keeps_the_dividing_lcm() {
        // Leads x^2, x*y, y^2 against the new lead x*y*z. The candidate lcms
        // are x^2*y*z, x*y*z, and x*y^2*z, so the middle one divides both
        // others and M drops them.
        let mut fixture = Fixture::new(3);
        fixture.push(&[2, 0, 0]);
        fixture.push(&[1, 1, 0]);
        fixture.push(&[0, 2, 0]);
        fixture.push(&[1, 1, 1]);
        let mut pairs = PairSet::new();
        fixture.update(&mut pairs, &[0, 1, 2], 3);

        assert_eq!(fixture.queued(&pairs), vec![(1, 3)]);
        let counters = pairs.counters();
        assert_eq!(counters.generated, 3);
        assert_eq!(counters.criterion_m, 2);
        assert_eq!(counters.criterion_f, 0);
    }

    #[test]
    fn criterion_f_keeps_the_smallest_index() {
        // Two live elements share a leading monomial, which only a retired
        // element can do, so both candidates have the same lcm.
        let mut fixture = Fixture::new(2);
        fixture.push(&[1, 0]);
        fixture.push(&[1, 0]);
        fixture.push(&[1, 1]);
        let mut pairs = PairSet::new();
        fixture.update(&mut pairs, &[0, 1], 2);

        assert_eq!(fixture.queued(&pairs), vec![(0, 2)]);
        let counters = pairs.counters();
        assert_eq!(counters.generated, 2);
        assert_eq!(counters.criterion_f, 1);
    }

    #[test]
    fn criterion_b_deletes_a_covered_pair() {
        // Leads x^2*y and x*y^2, then x*y. The queued pair (0, 1) has lcm
        // x^2*y^2, which x*y divides, and the two new lcms x^2*y and x*y^2
        // are proper divisors of it.
        let mut fixture = Fixture::new(2);
        fixture.push(&[2, 1]);
        fixture.push(&[1, 2]);
        let mut pairs = PairSet::new();
        fixture.update(&mut pairs, &[0], 1);
        assert_eq!(fixture.queued(&pairs), vec![(0, 1)]);

        fixture.push(&[1, 1]);
        fixture.update(&mut pairs, &[0, 1], 2);

        assert_eq!(fixture.queued(&pairs), vec![(0, 2), (1, 2)]);
        assert_eq!(pairs.counters().criterion_b, 1);
    }

    #[test]
    fn the_coverage_guard_keeps_a_pair_with_a_retired_endpoint() {
        // The same leads as above, except that element 1 is retired before
        // element 2 arrives. It produces no candidate, so the pair (0, 1)
        // has no cover and B must keep it.
        let mut fixture = Fixture::new(2);
        fixture.push(&[2, 1]);
        fixture.push(&[1, 2]);
        let mut pairs = PairSet::new();
        fixture.update(&mut pairs, &[0], 1);

        fixture.push(&[1, 1]);
        fixture.update(&mut pairs, &[0], 2);

        assert_eq!(fixture.queued(&pairs), vec![(0, 1), (0, 2)]);
        assert_eq!(pairs.counters().criterion_b, 0);
    }

    #[test]
    fn take_lowest_degree_returns_one_degree_class_in_order() {
        let mut fixture = Fixture::new(2);
        fixture.push(&[3, 0]);
        fixture.push(&[0, 3]);
        fixture.push(&[2, 1]);
        fixture.push(&[1, 2]);
        let mut pairs = PairSet::new();
        fixture.update(&mut pairs, &[0], 1);
        fixture.update(&mut pairs, &[0, 1], 2);
        fixture.update(&mut pairs, &[0, 1, 2], 3);

        let mut degrees: Vec<u32> = pairs.pairs().iter().map(|p| p.deg).collect();
        degrees.sort_unstable();
        let lowest = degrees[0];
        let expected = degrees.iter().filter(|&&d| d == lowest).count();

        let batch = pairs
            .take_lowest_degree(&fixture.table, &SelectOptions::default())
            .unwrap();
        assert_eq!(batch.len(), expected);
        assert!(batch.iter().all(|p| p.deg == lowest));
        for window in batch.windows(2) {
            let (a, b) = (window[0], window[1]);
            let key = fixture.table.cmp(a.lcm, b.lcm).then(a.i.cmp(&b.i));
            assert_ne!(key, std::cmp::Ordering::Greater, "the batch is sorted");
        }
        assert!(pairs.pairs().iter().all(|p| p.deg > lowest));
    }

    #[test]
    fn the_batch_cap_defers_the_rest_of_the_degree_class() {
        let mut fixture = Fixture::new(2);
        fixture.push(&[3, 0]);
        fixture.push(&[0, 3]);
        fixture.push(&[2, 1]);
        fixture.push(&[1, 2]);
        let mut pairs = PairSet::new();
        fixture.update(&mut pairs, &[0], 1);
        fixture.update(&mut pairs, &[0, 1], 2);
        fixture.update(&mut pairs, &[0, 1, 2], 3);
        let total = pairs.len();

        let options = SelectOptions {
            max_batch_pairs: Some(1),
        };
        let full = {
            let mut copy = PairSet::new();
            copy.pairs = pairs.pairs().to_vec();
            copy.take_lowest_degree(&fixture.table, &SelectOptions::default())
                .unwrap()
        };
        let first = pairs.take_lowest_degree(&fixture.table, &options).unwrap();

        assert_eq!(first.len(), 1);
        assert_eq!(first[0], full[0]);
        assert_eq!(pairs.len(), total - 1);
        if full.len() > 1 {
            let second = pairs.take_lowest_degree(&fixture.table, &options).unwrap();
            assert_eq!(second[0], full[1]);
        }
    }

    #[test]
    fn an_empty_queue_yields_no_batch() {
        let mut pairs = PairSet::new();
        let table = MonomialTable::<Lanes8>::new(2);
        assert!(
            pairs
                .take_lowest_degree(&table, &SelectOptions::default())
                .is_none()
        );
    }

    #[test]
    fn the_update_clears_the_workspace() {
        let mut fixture = Fixture::new(2);
        fixture.push(&[2, 0]);
        fixture.push(&[1, 1]);
        let mut pairs = PairSet::new();
        fixture.update(&mut pairs, &[0], 1);
        assert_eq!(fixture.ws.len(), 1, "only the identity is left");
    }
}
