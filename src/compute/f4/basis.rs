//! The working basis, redundancy, and lead-divisor index (design sections 4
//! and 5).
//!
//! [`Basis`] owns the basis table of design section 3.4, so every
//! [`MonomialId`] a basis element holds belongs to that table.

use rustc_hash::FxHashSet;

use super::field::FieldOps;
use super::monomial::{Lanes, MonomialId, MonomialTable, intern_from};
use super::pairs::PairSet;
use super::symbolic::Reducer;
use super::{Deadline, F4Error};

/// One generator, interned into the basis table.
///
/// The index of the generator in the caller's list is its public input
/// index, which [`Basis::seed`] records in the source map. `monos` is
/// strictly descending under the table order and is empty for the zero
/// polynomial. `coeffs` is parallel to `monos`.
pub(crate) struct InputPoly {
    /// The monomials, strictly descending.
    pub(crate) monos: Vec<MonomialId>,
    /// The coefficients, parallel to `monos`, each below `p`.
    pub(crate) coeffs: Vec<u32>,
}

/// One polynomial of the basis.
///
/// The polynomial is monic and its monomials are strictly descending, so
/// `monos[0]` is the leading monomial.
pub(crate) struct BasisPoly {
    /// The monomials, strictly descending, in the basis table.
    pub(crate) monos: Vec<MonomialId>,
    /// The coefficients, parallel to `monos`.
    pub(crate) coeffs: Vec<u32>,
    /// The multiply precomputation of `coeffs`, or empty under a kernel
    /// that needs none. Rebuild it after any change to `coeffs`.
    pub(crate) shoup: Vec<u64>,
}

impl BasisPoly {
    /// The number of terms.
    pub(crate) fn len(&self) -> usize {
        self.monos.len()
    }

    /// The leading monomial.
    pub(crate) fn lead(&self) -> MonomialId {
        self.monos[0]
    }

    /// The bytes the polynomial holds, counting capacities.
    fn memory_bytes(&self) -> usize {
        self.monos
            .capacity()
            .saturating_mul(size_of::<MonomialId>())
            .saturating_add(self.coeffs.capacity().saturating_mul(size_of::<u32>()))
            .saturating_add(self.shoup.capacity().saturating_mul(size_of::<u64>()))
    }
}

/// One element of the basis.
struct BasisElem {
    poly: BasisPoly,
    redundant: bool,
}

/// One entry of the lead-divisor index.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LeadEntry {
    lm: MonomialId,
    /// The divisor mask of `lm`, in the table `lm` belongs to.
    mask: u32,
    deg: u32,
    #[cfg(test)]
    nterms: u32,
    basis: u32,
}

/// The lead-divisor index over the live basis elements.
///
/// Buckets are sparse and sorted by degree. A dense vector indexed by
/// degree would need one slot for every degree up to the largest lead, so
/// one high-degree lead would force a vector of billions of empty slots.
///
/// The index holds ids of the basis table. Symbolic preprocessing asks
/// about a monomial of the batch-local table, and an id and a mask are
/// both table-local, so the scan runs on a [`LeadView`] that
/// [`LeadIndex::project`] builds inside the asking table.
#[derive(Default)]
pub(crate) struct LeadIndex {
    buckets: Vec<(u32, Vec<LeadEntry>)>,
}

impl LeadIndex {
    /// An empty index.
    pub(crate) fn new() -> Self {
        LeadIndex {
            buckets: Vec::new(),
        }
    }

    /// The number of entries.
    pub(crate) fn len(&self) -> usize {
        self.buckets.iter().map(|(_, entries)| entries.len()).sum()
    }

    /// Drop every entry.
    pub(crate) fn clear(&mut self) {
        self.buckets.clear();
    }

    /// Add one entry, keeping the buckets sorted by degree.
    fn insert(&mut self, entry: LeadEntry) -> Result<(), F4Error> {
        match self
            .buckets
            .binary_search_by_key(&entry.deg, |&(deg, _)| deg)
        {
            Ok(at) => {
                reserve(&mut self.buckets[at].1, 1)?;
                self.buckets[at].1.push(entry);
            }
            Err(at) => {
                reserve(&mut self.buckets, 1)?;
                self.buckets.insert(at, (entry.deg, vec![entry]));
            }
        }
        Ok(())
    }

    /// Remove the entry of basis element `basis`, whose lead has degree
    /// `deg`.
    fn remove(&mut self, basis: u32, deg: u32) {
        let Ok(at) = self.buckets.binary_search_by_key(&deg, |&(d, _)| d) else {
            return;
        };
        self.buckets[at].1.retain(|entry| entry.basis != basis);
        if self.buckets[at].1.is_empty() {
            self.buckets.remove(at);
        }
    }

    /// Copy the index into `query`, so a scan runs inside one table.
    ///
    /// The copy costs one intern per live element per batch. In exchange
    /// the scan runs inside one table, on that table's masks and packed
    /// words, instead of comparing exponents variable by variable across
    /// two tables.
    pub(crate) fn project<L: Lanes>(
        &self,
        basis_table: &MonomialTable<L>,
        query: &mut MonomialTable<L>,
    ) -> Result<LeadView, F4Error> {
        let mut entries = Vec::new();
        reserve(&mut entries, self.len())?;
        for (deg, bucket) in &self.buckets {
            for entry in bucket {
                let lm = intern_from(query, basis_table, entry.lm)?;
                entries.push(LeadEntry {
                    lm,
                    mask: query.mask(lm),
                    deg: *deg,
                    #[cfg(test)]
                    nterms: entry.nterms,
                    basis: entry.basis,
                });
            }
        }
        Ok(LeadView { entries })
    }

    /// The bytes the index holds, counting capacities.
    pub(crate) fn memory_bytes(&self) -> usize {
        self.buckets
            .capacity()
            .saturating_mul(size_of::<(u32, Vec<LeadEntry>)>())
            .saturating_add(
                self.buckets
                    .iter()
                    .map(|(_, entries)| entries.capacity().saturating_mul(size_of::<LeadEntry>()))
                    .fold(0, usize::saturating_add),
            )
    }
}

/// A reducer found for one monomial.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Divisor {
    pub(crate) basis: u32,
    /// The multiplier `m / lm(basis)`, in the asking table.
    pub(crate) mult: MonomialId,
}

/// The lead-divisor index projected into one batch's table.
///
/// Entries are in scan order: degree ascending, and basis index ascending
/// inside one degree. [`Reducer::First`] takes the first divisor in that
/// order.
pub(crate) struct LeadView {
    entries: Vec<LeadEntry>,
}

impl LeadView {
    /// The number of entries.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }

    /// Find a divisor of `m` among the live basis elements.
    ///
    /// `m` and every entry belong to `query`. The result's multiplier is
    /// interned into `query`. A degree above `deg(m)` ends the scan,
    /// because the entries are sorted by degree.
    pub(crate) fn find_divisor<L: Lanes>(
        &self,
        query: &mut MonomialTable<L>,
        m: MonomialId,
        reducer: Reducer,
    ) -> Result<Option<Divisor>, F4Error> {
        let deg = query.degree(m);
        let mask = query.mask(m);
        let mut best: Option<LeadEntry> = None;
        for entry in &self.entries {
            if entry.deg > deg {
                break;
            }
            // Every threshold the divisor passes, the multiple must pass
            // too, so this test rejects without touching the words.
            if entry.mask & !mask != 0 || !query.divides(entry.lm, m) {
                continue;
            }
            match reducer {
                Reducer::First => {
                    best = Some(*entry);
                    break;
                }
                #[cfg(test)]
                Reducer::Shortest => {
                    let better = best
                        .is_none_or(|old| (entry.nterms, entry.basis) < (old.nterms, old.basis));
                    if better {
                        best = Some(*entry);
                    }
                }
            }
        }
        match best {
            None => Ok(None),
            Some(entry) => {
                let mult = query
                    .quotient(m, entry.lm)?
                    .expect("the scan accepted this entry as a divisor");
                Ok(Some(Divisor {
                    basis: entry.basis,
                    mult,
                }))
            }
        }
    }
}

/// The working basis of one run.
pub(crate) struct Basis<L: Lanes> {
    table: MonomialTable<L>,
    elems: Vec<BasisElem>,
    leads: Vec<MonomialId>,
    source: Vec<Option<u32>>,
    live: Vec<u32>,
    lead: LeadIndex,
    unit: bool,
}

impl<L: Lanes> Basis<L> {
    /// An empty basis over `table`.
    ///
    /// The caller interns its generators into `table` first, because a
    /// basis takes them as [`InputPoly`] values over that table.
    pub(crate) fn new(table: MonomialTable<L>) -> Self {
        Basis {
            table,
            elems: Vec::new(),
            leads: Vec::new(),
            source: Vec::new(),
            live: Vec::new(),
            lead: LeadIndex::new(),
            unit: false,
        }
    }

    /// The basis table.
    pub(crate) fn table(&self) -> &MonomialTable<L> {
        &self.table
    }

    /// The number of elements, live and retired.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.elems.len()
    }

    /// The polynomial of element `index`.
    pub(crate) fn poly(&self, index: u32) -> &BasisPoly {
        &self.elems[index as usize].poly
    }

    /// The leading monomial of element `index`.
    pub(crate) fn lead(&self, index: u32) -> MonomialId {
        self.leads[index as usize]
    }

    /// Report whether element `index` is retired.
    ///
    /// A retired element is never a reducer and never takes part in a new
    /// pair. Its queued pairs stay in the queue.
    #[cfg(test)]
    pub(crate) fn is_redundant(&self, index: u32) -> bool {
        self.elems[index as usize].redundant
    }

    /// The live element indices, ascending.
    pub(crate) fn live(&self) -> &[u32] {
        &self.live
    }

    /// The public input index of element `index`, or `None` when the
    /// element comes from a batch.
    pub(crate) fn source(&self, index: u32) -> Option<u32> {
        self.source[index as usize]
    }

    /// The lead-divisor index over the live elements.
    pub(crate) fn lead_index(&self) -> &LeadIndex {
        &self.lead
    }

    /// Report whether the basis is the unit ideal.
    ///
    /// The run stops as soon as this holds. The live set is then the one
    /// constant element (design section 4).
    pub(crate) fn is_unit(&self) -> bool {
        self.unit
    }

    /// Seed the basis from the generators (design section 10).
    ///
    /// The index of a generator in `inputs` is its public input index.
    /// Zero generators and exact duplicates are dropped, each survivor is
    /// made monic, and the survivors are inserted in increasing
    /// leading-monomial order. A nonzero constant generator gives the unit
    /// basis and stops the seeding.
    pub(crate) fn seed<F: FieldOps<Coeff = u32>>(
        &mut self,
        field: &F,
        inputs: Vec<InputPoly>,
        pairs: &mut PairSet,
        ws: &mut MonomialTable<L>,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        debug_assert!(
            self.elems.is_empty(),
            "seeding runs once, on an empty basis"
        );
        let mut kept: Vec<(u32, BasisPoly)> = Vec::new();
        let mut seen: FxHashSet<(Vec<MonomialId>, Vec<u32>)> = FxHashSet::default();
        for (index, input) in inputs.into_iter().enumerate() {
            if input.monos.is_empty() {
                continue;
            }
            let poly = monic(field, input);
            if !seen.insert((poly.monos.clone(), poly.coeffs.clone())) {
                continue;
            }
            kept.push((index as u32, poly));
        }
        kept.sort_by(|a, b| self.table.cmp(a.1.lead(), b.1.lead()).then(a.0.cmp(&b.0)));
        for (public, poly) in kept {
            let index = self.elems.len() as u32;
            self.insert(poly, pairs, ws, clock)?;
            self.source[index as usize] = Some(public);
            // The identity is the smallest monomial, so a constant
            // generator sorts first and stops the seeding here.
            if self.unit {
                break;
            }
        }
        Ok(())
    }

    /// Insert one element and generate its pairs (design sections 4 and 5).
    ///
    /// `poly` is monic and strictly descending. Pair generation runs
    /// before any retirement, so the pair that justifies retiring an
    /// element exists before the element retires. With `g = x^2 + y` and
    /// `h = x`, the pair `(g, h)` is what produces `y`, and generating
    /// after retirement would lose it.
    pub(crate) fn insert(
        &mut self,
        poly: BasisPoly,
        pairs: &mut PairSet,
        ws: &mut MonomialTable<L>,
        clock: &mut Deadline,
    ) -> Result<(), F4Error> {
        debug_assert!(!poly.monos.is_empty(), "a basis element is nonzero");
        debug_assert_eq!(poly.monos.len(), poly.coeffs.len(), "parallel vectors");
        debug_assert_eq!(poly.coeffs[0], 1, "a basis element is monic");
        debug_assert!(
            poly.monos
                .windows(2)
                .all(|w| self.table.cmp(w[0], w[1]).is_gt()),
            "the monomials are strictly descending"
        );

        let index = self.elems.len() as u32;
        let lm = poly.lead();
        reserve(&mut self.elems, 1)?;
        reserve(&mut self.leads, 1)?;
        reserve(&mut self.source, 1)?;
        self.leads.push(lm);
        self.source.push(None);
        self.elems.push(BasisElem {
            poly,
            redundant: false,
        });

        if lm == MonomialId::ONE {
            return self.take_unit(index, pairs);
        }

        pairs.update(&mut self.table, ws, &self.leads, &self.live, index, clock)?;
        self.retire(index)
    }

    /// Insert several elements of one batch, smallest lead first.
    ///
    /// Two new elements of one batch can have leads that divide each
    /// other, because two non-pivot columns can divide one another. Each
    /// element therefore goes through the full insertion, one at a time.
    ///
    /// `assigned`, when present, receives the global basis index each
    /// element of `polys` got, in the order of `polys`, and `None` for an
    /// element the unit path stopped before (contract section 9.4, report
    /// 3).
    pub(crate) fn insert_batch(
        &mut self,
        polys: Vec<BasisPoly>,
        pairs: &mut PairSet,
        ws: &mut MonomialTable<L>,
        clock: &mut Deadline,
        assigned: Option<&mut Vec<Option<u32>>>,
    ) -> Result<(), F4Error> {
        let Some(assigned) = assigned else {
            let mut polys = polys;
            polys.sort_by(|a, b| self.table.cmp(a.lead(), b.lead()));
            for poly in polys {
                self.insert(poly, pairs, ws, clock)?;
                if self.unit {
                    break;
                }
            }
            return Ok(());
        };
        // The recorder names an element by its place in `polys`, so the
        // insertion order lives in a permutation and `polys` keeps its own
        // order.
        let mut order: Vec<usize> = (0..polys.len()).collect();
        order.sort_by(|&a, &b| self.table.cmp(polys[a].lead(), polys[b].lead()));
        let mut polys: Vec<Option<BasisPoly>> = polys.into_iter().map(Some).collect();
        assigned.clear();
        assigned.resize(polys.len(), None);
        for &at in &order {
            let poly = polys[at].take().expect("a position is taken once");
            let index = self.elems.len() as u32;
            self.insert(poly, pairs, ws, clock)?;
            assigned[at] = Some(index);
            if self.unit {
                break;
            }
        }
        Ok(())
    }

    /// Intern one reduced row into the basis table.
    ///
    /// `monos` holds the row's monomials in `from`, strictly descending,
    /// and `coeffs` holds its coefficients, monic and parallel to `monos`.
    /// `from` is the batch-local table, whose ids mean nothing here, so
    /// every monomial is interned again. The multiply precomputation is
    /// built here, before the element reduces anything.
    pub(crate) fn adopt_row<F: FieldOps<Coeff = u32>>(
        &mut self,
        from: &MonomialTable<L>,
        monos: &[MonomialId],
        coeffs: &[u32],
        field: &F,
    ) -> Result<BasisPoly, F4Error> {
        debug_assert_eq!(monos.len(), coeffs.len(), "parallel vectors");
        debug_assert_eq!(coeffs[0], 1, "a new element is monic");
        let mut interned = Vec::new();
        reserve(&mut interned, monos.len())?;
        for &mono in monos {
            interned.push(intern_from(&mut self.table, from, mono)?);
        }
        let mut owned = Vec::new();
        reserve(&mut owned, coeffs.len())?;
        owned.extend_from_slice(coeffs);
        let mut shoup = Vec::new();
        field.precompute(&owned, &mut shoup);
        Ok(BasisPoly {
            monos: interned,
            coeffs: owned,
            shoup,
        })
    }

    /// The bytes the basis holds, counting capacities.
    pub(crate) fn memory_bytes(&self) -> usize {
        let elems = self
            .elems
            .iter()
            .map(|elem| elem.poly.memory_bytes())
            .fold(0, usize::saturating_add);
        self.table
            .memory_bytes()
            .saturating_add(elems)
            .saturating_add(self.elems.capacity().saturating_mul(size_of::<BasisElem>()))
            .saturating_add(
                self.leads
                    .capacity()
                    .saturating_mul(size_of::<MonomialId>()),
            )
            .saturating_add(
                self.source
                    .capacity()
                    .saturating_mul(size_of::<Option<u32>>()),
            )
            .saturating_add(self.live.capacity().saturating_mul(size_of::<u32>()))
            .saturating_add(self.lead.memory_bytes())
    }

    /// Retire elements, using leading monomials only (step 3).
    fn retire(&mut self, index: u32) -> Result<(), F4Error> {
        let lm = self.leads[index as usize];
        // Equality counts, so two live elements never share a leading
        // monomial.
        let covered = self
            .live
            .iter()
            .any(|&g| self.table.divides(self.leads[g as usize], lm));
        if covered {
            self.elems[index as usize].redundant = true;
            return Ok(());
        }
        for &g in &self.live {
            let lm_g = self.leads[g as usize];
            if lm_g != lm && self.table.divides(lm, lm_g) {
                self.elems[g as usize].redundant = true;
                self.lead.remove(g, self.table.degree(lm_g));
            }
        }
        let elems = &self.elems;
        self.live.retain(|&g| !elems[g as usize].redundant);
        reserve(&mut self.live, 1)?;
        self.live.push(index);
        self.lead.insert(LeadEntry {
            lm,
            mask: self.table.mask(lm),
            deg: self.table.degree(lm),
            #[cfg(test)]
            nterms: self.elems[index as usize].poly.len() as u32,
            basis: index,
        })
    }

    /// Take the unit basis, because element `index` is a nonzero constant.
    fn take_unit(&mut self, index: u32, pairs: &mut PairSet) -> Result<(), F4Error> {
        self.unit = true;
        for elem in &mut self.elems {
            elem.redundant = true;
        }
        self.elems[index as usize].redundant = false;
        self.live.clear();
        reserve(&mut self.live, 1)?;
        self.live.push(index);
        self.lead.clear();
        // The identity divides every leading monomial, so every other
        // element is retired and no pair can add anything.
        pairs.clear();
        self.lead.insert(LeadEntry {
            lm: MonomialId::ONE,
            mask: self.table.mask(MonomialId::ONE),
            deg: 0,
            #[cfg(test)]
            nterms: 1,
            basis: index,
        })
    }
}

/// Scale a generator to make it monic.
fn monic<F: FieldOps<Coeff = u32>>(field: &F, input: InputPoly) -> BasisPoly {
    let InputPoly { monos, mut coeffs } = input;
    if coeffs[0] != 1 {
        let inv = field.inv(coeffs[0]);
        for c in &mut coeffs {
            // Both factors are below 2^31, so the product fits the
            // accumulator lane `reduce_acc` takes.
            *c = field.reduce_acc(u64::from(*c) * u64::from(inv));
        }
    }
    let mut shoup = Vec::new();
    field.precompute(&coeffs, &mut shoup);
    BasisPoly {
        monos,
        coeffs,
        shoup,
    }
}

/// Grow `vec` by `extra` slots, or report an exhausted memory budget.
///
/// Every growth of an engine structure goes through this, so a failed
/// allocation is a typed budget report and not a panic (design section
/// 11).
pub(crate) fn reserve<T>(vec: &mut Vec<T>, extra: usize) -> Result<(), F4Error> {
    vec.try_reserve(extra)
        .map_err(|_| F4Error::MemoryLimitExceeded)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compute::f4::field::Small31;
    use crate::compute::f4::monomial::Lanes8;
    use crate::compute::f4::pairs::PairSet;

    const P: u32 = 32003;

    /// A basis under construction, with the workspace table beside it.
    struct Fixture {
        basis: Basis<Lanes8>,
        pairs: PairSet,
        ws: MonomialTable<Lanes8>,
        field: Small31,
    }

    impl Fixture {
        fn new(nvars: usize) -> Self {
            Fixture {
                basis: Basis::new(MonomialTable::new(nvars)),
                pairs: PairSet::new(),
                ws: MonomialTable::new(nvars),
                field: Small31::new(P),
            }
        }

        /// Intern one term list into the basis table.
        fn input(&mut self, terms: &[(u32, &[u32])]) -> InputPoly {
            let mut monos = Vec::new();
            let mut coeffs = Vec::new();
            for (coeff, exps) in terms {
                monos.push(self.basis.table.intern(exps).unwrap());
                coeffs.push(*coeff);
            }
            InputPoly { monos, coeffs }
        }

        /// Insert a monic polynomial with the given monomials.
        fn insert(&mut self, exps: &[&[u32]]) {
            let mut monos = Vec::new();
            for e in exps {
                monos.push(self.basis.table.intern(e).unwrap());
            }
            let coeffs = vec![1; monos.len()];
            let mut shoup = Vec::new();
            self.field.precompute(&coeffs, &mut shoup);
            self.basis
                .insert(
                    BasisPoly {
                        monos,
                        coeffs,
                        shoup,
                    },
                    &mut self.pairs,
                    &mut self.ws,
                    &mut Deadline::none(),
                )
                .unwrap();
        }

        fn queued(&self) -> Vec<(u32, u32)> {
            let mut out: Vec<(u32, u32)> = self.pairs.pairs().iter().map(|p| (p.i, p.j)).collect();
            out.sort_unstable();
            out
        }
    }

    #[test]
    fn pairs_are_generated_before_the_retirement_they_justify() {
        // g = x^2 + y and h = x. Retiring g on the grounds that lm(h)
        // divides lm(g) loses y, which only the pair (g, h) produces.
        let mut fixture = Fixture::new(2);
        fixture.insert(&[&[2, 0], &[0, 1]]);
        fixture.insert(&[&[1, 0]]);

        assert_eq!(fixture.queued(), vec![(0, 1)]);
        assert!(fixture.basis.is_redundant(0));
        assert!(!fixture.basis.is_redundant(1));
        assert_eq!(fixture.basis.live(), &[1]);
        assert_eq!(fixture.basis.lead_index().len(), 1);
    }

    #[test]
    fn a_new_element_under_a_live_lead_retires_itself() {
        let mut fixture = Fixture::new(2);
        fixture.insert(&[&[1, 0]]);
        fixture.insert(&[&[2, 0], &[0, 1]]);

        assert!(fixture.basis.is_redundant(1));
        assert_eq!(fixture.basis.live(), &[0]);
    }

    #[test]
    fn two_live_elements_never_share_a_lead() {
        let mut fixture = Fixture::new(2);
        fixture.insert(&[&[1, 1], &[1, 0]]);
        fixture.insert(&[&[1, 1], &[0, 1]]);

        assert_eq!(fixture.basis.live(), &[0]);
        assert!(fixture.basis.is_redundant(1));
        // The pair of the two equal leads is generated before the
        // retirement, so the element that retires is still reachable.
        assert_eq!(fixture.queued(), vec![(0, 1)]);
    }

    #[test]
    fn a_constant_generator_gives_the_unit_basis() {
        let mut fixture = Fixture::new(2);
        let inputs = vec![
            fixture.input(&[(3, &[2, 0]), (1, &[0, 1])]),
            fixture.input(&[(7, &[0, 0])]),
        ];
        let mut pairs = PairSet::new();
        let mut ws = MonomialTable::new(2);
        fixture
            .basis
            .seed(
                &fixture.field,
                inputs,
                &mut pairs,
                &mut ws,
                &mut Deadline::none(),
            )
            .unwrap();

        assert!(fixture.basis.is_unit());
        assert_eq!(fixture.basis.live().len(), 1);
        let only = fixture.basis.live()[0];
        assert_eq!(fixture.basis.lead(only), MonomialId::ONE);
        assert_eq!(fixture.basis.poly(only).coeffs, vec![1]);
        assert_eq!(fixture.basis.source(only), Some(1));
        assert!(pairs.is_empty());
    }

    #[test]
    fn seeding_drops_zero_and_duplicate_generators() {
        let mut fixture = Fixture::new(2);
        let inputs = vec![
            // 2*x^2 + 2*y, which is the same polynomial as input 2 once
            // both are monic.
            fixture.input(&[(2, &[2, 0]), (2, &[0, 1])]),
            fixture.input(&[]),
            fixture.input(&[(5, &[2, 0]), (5, &[0, 1])]),
            fixture.input(&[(1, &[0, 3])]),
        ];
        let mut pairs = PairSet::new();
        let mut ws = MonomialTable::new(2);
        fixture
            .basis
            .seed(
                &fixture.field,
                inputs,
                &mut pairs,
                &mut ws,
                &mut Deadline::none(),
            )
            .unwrap();

        assert_eq!(fixture.basis.len(), 2);
        // Grevlex compares the degree first, so x^2 comes before y^3.
        assert_eq!(fixture.basis.source(0), Some(0));
        assert_eq!(fixture.basis.source(1), Some(3));
        assert_eq!(fixture.basis.poly(0).coeffs, vec![1, 1]);
    }

    #[test]
    fn seeding_makes_every_generator_monic() {
        let mut fixture = Fixture::new(1);
        let inputs = vec![fixture.input(&[(5, &[2]), (3, &[0])])];
        let mut pairs = PairSet::new();
        let mut ws = MonomialTable::new(1);
        fixture
            .basis
            .seed(
                &fixture.field,
                inputs,
                &mut pairs,
                &mut ws,
                &mut Deadline::none(),
            )
            .unwrap();

        let poly = fixture.basis.poly(0);
        assert_eq!(poly.coeffs[0], 1);
        // 3 / 5 mod 32003.
        let expected = fixture
            .field
            .reduce_acc(u64::from(3u32) * u64::from(fixture.field.inv(5)));
        assert_eq!(poly.coeffs[1], expected);
    }

    /// A small deterministic generator, so the tests carry no dependency.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
            self.0 >> 33
        }

        fn below(&mut self, bound: u64) -> u64 {
            self.next() % bound
        }
    }

    /// Random exponent vectors, each of total degree at most `max_deg`.
    fn random_monomials(rng: &mut Rng, count: usize, nvars: usize, max_deg: u32) -> Vec<Vec<u32>> {
        (0..count)
            .map(|_| {
                (0..nvars)
                    .map(|_| rng.below(u64::from(max_deg) + 1) as u32)
                    .collect()
            })
            .collect()
    }

    /// Random exponent vectors of total degree at least 1.
    ///
    /// A degree of 0 is the constant 1, which takes the unit basis and
    /// empties the queue, so a pair test must not draw it.
    fn random_leads(rng: &mut Rng, count: usize, nvars: usize, max_deg: u32) -> Vec<Vec<u32>> {
        let mut out = Vec::new();
        while out.len() < count {
            let candidate = random_monomials(rng, 1, nvars, max_deg).pop().unwrap();
            if candidate.iter().any(|&e| e > 0) {
                out.push(candidate);
            }
        }
        out
    }

    #[test]
    fn the_lead_index_finds_the_same_divisor_as_a_linear_scan() {
        let mut rng = Rng(0x5eed);
        for _ in 0..64 {
            let nvars = 3;
            let mut fixture = Fixture::new(nvars);
            let leads = random_leads(&mut rng, 8, nvars, 3);
            for lead in &leads {
                let tail: Vec<&[u32]> = vec![lead.as_slice()];
                fixture.insert(&tail);
            }

            let mut query = MonomialTable::<Lanes8>::new(nvars);
            let view = fixture
                .basis
                .lead_index()
                .project(fixture.basis.table(), &mut query)
                .unwrap();
            assert_eq!(view.len(), fixture.basis.live().len());

            for target in random_monomials(&mut rng, 16, nvars, 5) {
                let m = query.intern(&target).unwrap();
                // The scan order is degree ascending, then basis index
                // ascending, so the reference sorts the live divisors the
                // same way.
                let mut divisors: Vec<(u32, u32)> = fixture
                    .basis
                    .live()
                    .iter()
                    .filter(|&&g| {
                        let lm = &leads[g as usize];
                        lm.iter().zip(&target).all(|(&a, &b)| a <= b)
                    })
                    .map(|&g| (leads[g as usize].iter().sum::<u32>(), g))
                    .collect();
                divisors.sort_unstable();

                let found = view.find_divisor(&mut query, m, Reducer::First).unwrap();
                match divisors.first() {
                    None => assert!(found.is_none()),
                    Some(&(_, expected)) => {
                        let found = found.expect("a divisor exists");
                        assert_eq!(found.basis, expected);
                        let lm = fixture.basis.lead(expected);
                        let lm_query = intern_from(&mut query, fixture.basis.table(), lm).unwrap();
                        assert_eq!(query.mul(found.mult, lm_query).unwrap(), m);
                    }
                }
            }
        }
    }

    #[test]
    fn the_shortest_reducer_has_the_fewest_terms() {
        let mut fixture = Fixture::new(2);
        // Two live elements with leads x and y, of three and two terms.
        fixture.insert(&[&[1, 0], &[0, 1], &[0, 0]]);
        fixture.insert(&[&[0, 1], &[0, 0]]);

        let mut query = MonomialTable::<Lanes8>::new(2);
        let view = fixture
            .basis
            .lead_index()
            .project(fixture.basis.table(), &mut query)
            .unwrap();
        let m = query.intern(&[1, 1]).unwrap();

        let first = view
            .find_divisor(&mut query, m, Reducer::First)
            .unwrap()
            .unwrap();
        let shortest = view
            .find_divisor(&mut query, m, Reducer::Shortest)
            .unwrap()
            .unwrap();
        assert_eq!(first.basis, 0);
        assert_eq!(shortest.basis, 1);
    }

    /// The naive divisibility test over exponent vectors.
    fn divides(a: &[u32], b: &[u32]) -> bool {
        a.iter().zip(b).all(|(&x, &y)| x <= y)
    }

    /// The naive lcm over exponent vectors.
    fn lcm(a: &[u32], b: &[u32]) -> Vec<u32> {
        a.iter().zip(b).map(|(&x, &y)| x.max(y)).collect()
    }

    /// The naive coprimality test over exponent vectors.
    fn coprime(a: &[u32], b: &[u32]) -> bool {
        a.iter().zip(b).all(|(&x, &y)| x == 0 || y == 0)
    }

    /// Report whether every pair the update dropped is covered.
    ///
    /// A pair is covered when it is queued, when the product criterion
    /// applies, or when some element `k` has a leading monomial dividing
    /// the pair's lcm and both pairs through `k` are already covered. The
    /// closure grows from queued and coprime pairs only, so no pair ever
    /// justifies itself.
    fn every_pair_is_covered(leads: &[Vec<u32>], queued: &[(u32, u32)]) -> bool {
        let mut covered = coprime_pairs(leads);
        for &(i, j) in queued {
            covered[i as usize][j as usize] = true;
            covered[j as usize][i as usize] = true;
        }
        while grow_coverage(leads, &mut covered) {}
        all_pairs_covered(&covered)
    }

    fn coprime_pairs(leads: &[Vec<u32>]) -> Vec<Vec<bool>> {
        let count = leads.len();
        let mut covered = vec![vec![false; count]; count];
        for i in 0..count {
            for j in 0..count {
                covered[i][j] = i != j && coprime(&leads[i], &leads[j]);
            }
        }
        covered
    }

    fn grow_coverage(leads: &[Vec<u32>], covered: &mut [Vec<bool>]) -> bool {
        let mut grew = false;
        for i in 0..leads.len() {
            for j in 0..leads.len() {
                if uncovered_chain_exists(leads, covered, i, j) {
                    covered[i][j] = true;
                    covered[j][i] = true;
                    grew = true;
                }
            }
        }
        grew
    }

    fn uncovered_chain_exists(
        leads: &[Vec<u32>],
        covered: &[Vec<bool>],
        i: usize,
        j: usize,
    ) -> bool {
        if i == j || covered[i][j] {
            return false;
        }
        let multiple = lcm(&leads[i], &leads[j]);
        (0..leads.len()).any(|k| {
            k != i && k != j && divides(&leads[k], &multiple) && covered[i][k] && covered[k][j]
        })
    }

    fn all_pairs_covered(covered: &[Vec<bool>]) -> bool {
        (0..covered.len()).all(|i| (0..covered.len()).all(|j| i == j || covered[i][j]))
    }

    #[test]
    fn every_pair_the_update_drops_stays_covered() {
        let mut rng = Rng(0xf5);
        for round in 0..200 {
            let nvars = 2 + (round % 2);
            let mut fixture = Fixture::new(nvars);
            let leads = random_leads(&mut rng, 6, nvars, 3);
            for lead in &leads {
                fixture.insert(&[lead.as_slice()]);
            }

            let queued = fixture.queued();
            assert!(
                every_pair_is_covered(&leads, &queued),
                "uncovered pair for leads {leads:?} with queue {queued:?}"
            );

            let counters = fixture.pairs.counters();
            let queued_now = counters.generated
                - counters.product
                - counters.criterion_m
                - counters.criterion_f
                - counters.criterion_b;
            assert_eq!(queued_now, queued.len() as u64, "the counters balance");
        }
    }
}
