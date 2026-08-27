//! Interned monomials for the F4 engine (design section 3).
//!
//! A [`MonomialTable`] packs an exponent vector into `u64` words, hash-conses
//! it, and returns one [`MonomialId`] per distinct vector. The table is
//! append-only, so an id stays valid until [`MonomialTable::clear`] resets the
//! table. An id is table-local: it means nothing in another table.
//!
//! [`MonomialId`] has no [`Ord`], because numeric id order has nothing to do
//! with the monomial order. Comparison goes through [`MonomialTable::cmp`].

use std::cmp::Ordering;
use std::marker::PhantomData;

use hashbrown::HashTable;

use super::F4Error;

/// The largest exponent a monomial may hold.
///
/// This is the public contract of [`crate::ring::PolynomialRing`], and the F4
/// engine holds to it at every lane width.
pub(crate) const MAX_EXPONENT: u32 = u16::MAX as u32;

/// The seed of the linear hash of design section 3.1.
///
/// It is a constant, not a random value, so two runs of the same input agree
/// on every hash and on every id.
const HASH_SEED: u64 = 0x9e37_79b9_7f4a_7c15;

/// One packing of exponents into `u64` words.
///
/// A lane holds a value in `0 ..= BOUND`, and `BOUND` is `min(H - 1,
/// MAX_EXPONENT)` for `H` half the lane range. The high bit of every lane
/// stays clear. That reserved bit is what makes the word tricks in this
/// module exact, because it absorbs the borrow of a lanewise subtraction.
pub(crate) trait Lanes: Copy + Clone + Send + Sync + 'static {
    /// Bits per lane.
    const BITS: u32;
    /// Lanes per `u64` word.
    const PER_WORD: usize;
    /// The largest exponent one lane holds.
    const BOUND: u32;
    /// The high bit of every lane.
    const GUARD: u64;
    /// The low bit of every lane.
    const ONES: u64;
    /// Every bit of a lane above `BOUND`.
    const OVER: u64;
}

/// Eight-bit lanes, eight per word. Exponents up to 127.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Lanes8;

/// Sixteen-bit lanes, four per word. Exponents up to 32767.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Lanes16;

/// Thirty-two-bit lanes, two per word. Exponents up to 65535.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Lanes32;

impl Lanes for Lanes8 {
    const BITS: u32 = 8;
    const PER_WORD: usize = 8;
    const BOUND: u32 = 127;
    const GUARD: u64 = 0x8080_8080_8080_8080;
    const ONES: u64 = 0x0101_0101_0101_0101;
    const OVER: u64 = 0x8080_8080_8080_8080;
}

impl Lanes for Lanes16 {
    const BITS: u32 = 16;
    const PER_WORD: usize = 4;
    const BOUND: u32 = 32767;
    const GUARD: u64 = 0x8000_8000_8000_8000;
    const ONES: u64 = 0x0001_0001_0001_0001;
    const OVER: u64 = 0x8000_8000_8000_8000;
}

impl Lanes for Lanes32 {
    const BITS: u32 = 32;
    const PER_WORD: usize = 2;
    const BOUND: u32 = 65535;
    const GUARD: u64 = 0x8000_0000_8000_0000;
    const ONES: u64 = 0x0000_0001_0000_0001;
    const OVER: u64 = 0xffff_0000_ffff_0000;
}

/// One monomial of one [`MonomialTable`].
///
/// `MonomialId(0)` is the identity monomial of every table.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub(crate) struct MonomialId(u32);

impl MonomialId {
    /// The identity monomial.
    pub(crate) const ONE: MonomialId = MonomialId(0);

    /// The dense index of this id, for a vector indexed by monomial.
    pub(crate) fn index(self) -> usize {
        self.0 as usize
    }
}

/// The cached data of one interned monomial.
#[derive(Clone, Copy, Debug)]
struct MonoMeta {
    deg: u32,
    mask: u32,
    hash: u64,
}

/// A hash-consed store of monomials over `nvars` variables.
///
/// Construct with [`MonomialTable::new`] or
/// [`MonomialTable::with_max_exponent`]. Both start with the identity
/// monomial at [`MonomialId::ONE`]. Interning is deterministic: the same
/// insertion sequence gives the same ids, because the table appends and never
/// iterates its hash index.
pub(crate) struct MonomialTable<L: Lanes> {
    nvars: usize,
    words: usize,
    mask_bits: usize,
    thresholds: [u32; 32],
    coeffs: Vec<u64>,
    exps: Vec<u64>,
    meta: Vec<MonoMeta>,
    index: HashTable<u32>,
    scratch: Vec<u64>,
    limit: Option<usize>,
    lanes: PhantomData<L>,
}

impl<L: Lanes> MonomialTable<L> {
    /// A table over `nvars` variables, with the mask ladder `t_k = k`.
    ///
    /// This ladder suits an input whose exponents stay small. When the input
    /// carries larger exponents, use [`MonomialTable::with_max_exponent`].
    #[cfg(test)]
    pub(crate) fn new(nvars: usize) -> Self {
        Self::with_max_exponent(nvars, 0)
    }

    /// A table whose mask ladder spreads over `1 ..= max_exponent`.
    ///
    /// `max_exponent` is the largest exponent the input holds. The ladder is
    /// fixed here and never changes, so the mask stays a pure function of the
    /// exponent vector and the table.
    pub(crate) fn with_max_exponent(nvars: usize, max_exponent: u32) -> Self {
        let words = nvars.div_ceil(L::PER_WORD);
        let mask_bits = if nvars == 0 { 1 } else { (32 / nvars).max(1) };
        let mut thresholds = [0u32; 32];
        let step = (max_exponent / mask_bits as u32).max(1);
        for (k, t) in thresholds.iter_mut().enumerate() {
            *t = k as u32 * step;
        }
        let mut table = MonomialTable {
            nvars,
            words,
            mask_bits,
            thresholds,
            coeffs: hash_coefficients(nvars),
            exps: Vec::new(),
            meta: Vec::new(),
            index: HashTable::new(),
            scratch: vec![0; words],
            limit: None,
            lanes: PhantomData,
        };
        table.push_identity();
        table
    }

    /// Stop the table from growing past `limit` bytes.
    ///
    /// `limit` is the memory limit of the whole run. One table alone never
    /// holds more than the run does, so a projected table above the limit
    /// is a run above the limit. The check runs before the table grows its
    /// exponent slab, its metadata, or its hash index, so the table passes
    /// the limit by at most one growth step of one of the three.
    pub(crate) fn with_limit(mut self, limit: Option<usize>) -> Self {
        self.limit = limit;
        self
    }

    /// The number of distinct monomials the table holds.
    pub(crate) fn len(&self) -> usize {
        self.meta.len()
    }

    /// Drop every monomial except the identity, keeping the allocations.
    ///
    /// The update workspace and symbolic tables of design section 3.4 are
    /// reset this way. Every id the caller still holds becomes meaningless.
    pub(crate) fn clear(&mut self) {
        self.exps.clear();
        self.meta.clear();
        self.index.clear();
        self.push_identity();
    }

    /// Release every byte the table's allocations hold above their
    /// contents.
    ///
    /// The memory retry of design section 11 calls this after
    /// [`MonomialTable::clear`], because the estimate counts capacity and
    /// a retry that holds the failed attempt's capacity cannot fit a
    /// smaller batch.
    pub(crate) fn shrink(&mut self) {
        self.exps.shrink_to_fit();
        self.meta.shrink_to_fit();
        let meta = std::mem::take(&mut self.meta);
        self.index.shrink_to_fit(|&i| meta[i as usize].hash);
        self.meta = meta;
    }

    /// The bytes the table's allocations hold, counting capacities.
    pub(crate) fn memory_bytes(&self) -> usize {
        let word = size_of::<u64>();
        self.exps
            .capacity()
            .saturating_mul(word)
            .saturating_add(self.meta.capacity().saturating_mul(size_of::<MonoMeta>()))
            .saturating_add(self.scratch.capacity().saturating_mul(word))
            .saturating_add(self.coeffs.capacity().saturating_mul(word))
            .saturating_add(self.index.capacity().saturating_mul(size_of::<u32>() + 1))
    }

    /// Make room for one more monomial, inside the table's limit.
    ///
    /// The three containers grow together, so one check covers all three.
    fn reserve_one(&mut self) -> Result<(), F4Error> {
        let words = self.words;
        let grows = self.exps.len() + words > self.exps.capacity()
            || self.meta.len() == self.meta.capacity()
            || self.index.len() == self.index.capacity();
        if !grows {
            return Ok(());
        }
        if let Some(limit) = self.limit
            && self.projected_bytes() > limit
        {
            return Err(F4Error::MemoryLimitExceeded);
        }
        self.exps
            .try_reserve(words)
            .map_err(|_| F4Error::MemoryLimitExceeded)?;
        self.meta
            .try_reserve(1)
            .map_err(|_| F4Error::MemoryLimitExceeded)?;
        let meta = &self.meta;
        self.index
            .try_reserve(1, |&i| meta[i as usize].hash)
            .map_err(|_| F4Error::MemoryLimitExceeded)?;
        Ok(())
    }

    /// The bytes the table holds after one more growth step.
    fn projected_bytes(&self) -> usize {
        let word = size_of::<u64>();
        let exps = grown(self.exps.capacity(), self.exps.len() + self.words);
        let meta = grown(self.meta.capacity(), self.meta.len() + 1);
        let index = grown(self.index.capacity(), self.index.len() + 1);
        exps.saturating_mul(word)
            .saturating_add(meta.saturating_mul(size_of::<MonoMeta>()))
            .saturating_add(self.scratch.capacity().saturating_mul(word))
            .saturating_add(self.coeffs.capacity().saturating_mul(word))
            .saturating_add(index.saturating_mul(size_of::<u32>() + 1))
    }

    /// The total degree of `a`.
    pub(crate) fn degree(&self, a: MonomialId) -> u32 {
        self.meta[a.index()].deg
    }

    /// The divisor mask of `a`, for a caller that caches it next to an id.
    ///
    /// The mask belongs to this table's ladder. Never compare it against a
    /// mask from another table.
    pub(crate) fn mask(&self, a: MonomialId) -> u32 {
        self.meta[a.index()].mask
    }

    /// The linear hash of `a`.
    ///
    /// The hash is `sum_i r_i * e_i` with wrapping arithmetic, so it is
    /// additive: `hash(a * b) == hash(a) + hash(b)`. Two tables over the same
    /// variable count agree on it.
    #[cfg(test)]
    pub(crate) fn hash(&self, a: MonomialId) -> u64 {
        self.meta[a.index()].hash
    }

    /// Write the exponent vector of `a` into `out`, replacing its contents.
    pub(crate) fn unpack(&self, a: MonomialId, out: &mut Vec<u32>) {
        out.clear();
        let words = self.words_of(a);
        out.extend((0..self.nvars).map(|i| lane::<L>(words, i)));
    }

    /// Intern an exponent vector.
    ///
    /// `exps` holds one exponent per variable. An exponent above
    /// [`MAX_EXPONENT`] gives [`F4Error::ExponentOverflow`]. An exponent above
    /// the lane bound gives [`F4Error::LaneOverflow`].
    pub(crate) fn intern(&mut self, exps: &[u32]) -> Result<MonomialId, F4Error> {
        debug_assert_eq!(exps.len(), self.nvars, "exponent vectors must share nvars");
        let mut deg = 0u32;
        for &e in exps {
            if e > L::BOUND {
                return Err(width_error(e));
            }
            // A ring holds at most 256 variables and an exponent at most
            // 65535, so the sum stays below u32::MAX (design section 3.2).
            deg += e;
        }
        let mut hash = 0u64;
        for (r, &e) in self.coeffs.iter().zip(exps) {
            hash = hash.wrapping_add(r.wrapping_mul(e as u64));
        }
        self.scratch.fill(0);
        for (i, &e) in exps.iter().enumerate() {
            set_lane::<L>(&mut self.scratch, i, e);
        }
        self.intern_scratch(deg, hash)
    }

    /// Compare two monomials under grevlex, the larger monomial greater.
    pub(crate) fn cmp(&self, a: MonomialId, b: MonomialId) -> Ordering {
        // Hash-consing makes ids unique, so equal ids are equal monomials.
        if a == b {
            return Ordering::Equal;
        }
        let (ia, ib) = (a.index(), b.index());
        match self.meta[ia].deg.cmp(&self.meta[ib].deg) {
            Ordering::Equal => self.tail_cmp(ia, ib),
            ord => ord,
        }
    }

    /// Compare the grevlex tails of two monomials of equal degree.
    ///
    /// The table stores no separate key slab. The reversed scan below gives
    /// the grevlex tail order from the exponent words alone (design section
    /// 3.1).
    fn tail_cmp(&self, ia: usize, ib: usize) -> Ordering {
        let wa = self.words_at(ia);
        let wb = self.words_at(ib);
        // Variable i sits in lane i % PER_WORD of word i / PER_WORD, so a
        // higher variable index sits in a higher bit. The first difference
        // found from the last word down is therefore the highest differing
        // variable, which is what grevlex compares, and the larger exponent
        // there makes the smaller monomial. Variable 0 is the last lane
        // scanned. Equal degree and equal higher lanes force it equal, so it
        // never decides.
        for k in (0..self.words).rev() {
            match wa[k].cmp(&wb[k]) {
                Ordering::Equal => continue,
                ord => return ord.reverse(),
            }
        }
        Ordering::Equal
    }

    /// Report whether `a` divides `b`.
    pub(crate) fn divides(&self, a: MonomialId, b: MonomialId) -> bool {
        let (ma, mb) = (&self.meta[a.index()], &self.meta[b.index()]);
        if ma.deg > mb.deg || ma.mask & !mb.mask != 0 {
            return false;
        }
        divides_words::<L>(self.words_at(a.index()), self.words_at(b.index()))
    }

    /// Report whether no variable appears in both monomials.
    pub(crate) fn is_coprime(&self, a: MonomialId, b: MonomialId) -> bool {
        let wa = self.words_at(a.index());
        let wb = self.words_at(b.index());
        wa.iter()
            .zip(wb)
            .all(|(&x, &y)| nonzero_lanes::<L>(x) & nonzero_lanes::<L>(y) == 0)
    }

    /// Intern the product `a * b`.
    #[cfg(test)]
    pub(crate) fn mul(&mut self, a: MonomialId, b: MonomialId) -> Result<MonomialId, F4Error> {
        let (ia, ib) = (a.index() * self.words, b.index() * self.words);
        let mut over = 0u64;
        for k in 0..self.words {
            // Both operands hold every lane below half the lane range, so a
            // lanewise sum never carries into the next lane.
            let w = self.exps[ia + k] + self.exps[ib + k];
            over |= w & L::OVER;
            self.scratch[k] = w;
        }
        if over != 0 {
            return Err(self.overflow_error());
        }
        let deg = self.meta[a.index()].deg + self.meta[b.index()].deg;
        let hash = self.meta[a.index()]
            .hash
            .wrapping_add(self.meta[b.index()].hash);
        self.intern_scratch(deg, hash)
    }

    /// Intern the quotient `a / b`, or report that `b` does not divide `a`.
    pub(crate) fn quotient(
        &mut self,
        a: MonomialId,
        b: MonomialId,
    ) -> Result<Option<MonomialId>, F4Error> {
        if !self.divides(b, a) {
            return Ok(None);
        }
        let (ia, ib) = (a.index() * self.words, b.index() * self.words);
        for k in 0..self.words {
            // Every lane of b is at most the matching lane of a, so no borrow
            // crosses a lane.
            let w = self.exps[ia + k] - self.exps[ib + k];
            self.scratch[k] = w;
        }
        let deg = self.meta[a.index()].deg - self.meta[b.index()].deg;
        let hash = self.meta[a.index()]
            .hash
            .wrapping_sub(self.meta[b.index()].hash);
        self.intern_scratch(deg, hash).map(Some)
    }

    /// Intern the least common multiple of `a` and `b`.
    #[cfg(test)]
    pub(crate) fn lcm(&mut self, a: MonomialId, b: MonomialId) -> Result<MonomialId, F4Error> {
        let (ia, ib) = (a.index() * self.words, b.index() * self.words);
        for k in 0..self.words {
            let w = lane_max::<L>(self.exps[ia + k], self.exps[ib + k]);
            self.scratch[k] = w;
        }
        let (deg, hash) = self.degree_and_hash_of_scratch();
        self.intern_scratch(deg, hash)
    }

    /// Intern the product of `lhs` in this table and `rhs_id` in `rhs`.
    ///
    /// Matrix construction runs this once per term of every reducer, with the
    /// multiplier in this table and the basis monomial in the basis table.
    /// Both tables must hold the same variable count. Their hash coefficients
    /// then agree, so this adds the two cached hashes instead of reading
    /// exponents.
    pub(crate) fn mul_external(
        &mut self,
        lhs: MonomialId,
        rhs: &Self,
        rhs_id: MonomialId,
    ) -> Result<MonomialId, F4Error> {
        debug_assert_eq!(self.nvars, rhs.nvars, "tables must share nvars");
        let base = lhs.index() * self.words;
        let other = rhs.words_at(rhs_id.index());
        let mut over = 0u64;
        for (k, &rhs_word) in other.iter().enumerate() {
            let w = self.exps[base + k] + rhs_word;
            over |= w & L::OVER;
            self.scratch[k] = w;
        }
        if over != 0 {
            return Err(self.overflow_error());
        }
        let deg = self.meta[lhs.index()].deg + rhs.meta[rhs_id.index()].deg;
        let hash = self.meta[lhs.index()]
            .hash
            .wrapping_add(rhs.meta[rhs_id.index()].hash);
        self.intern_scratch(deg, hash)
    }

    fn push_identity(&mut self) {
        debug_assert!(self.exps.is_empty());
        self.exps.resize(self.words, 0);
        self.meta.push(MonoMeta {
            deg: 0,
            mask: 0,
            hash: 0,
        });
        let meta = &self.meta;
        self.index.insert_unique(0, 0, |&i| meta[i as usize].hash);
    }

    fn words_of(&self, a: MonomialId) -> &[u64] {
        self.words_at(a.index())
    }

    fn words_at(&self, index: usize) -> &[u64] {
        let start = index * self.words;
        &self.exps[start..start + self.words]
    }

    /// Name the overflow the `OVER` mask caught, from the first lane past the
    /// bound.
    fn overflow_error(&self) -> F4Error {
        for i in 0..self.nvars {
            let e = lane::<L>(&self.scratch, i);
            if e > L::BOUND {
                return width_error(e);
            }
        }
        debug_assert!(false, "overflow mask fired with every lane in bounds");
        F4Error::LaneOverflow
    }

    fn degree_and_hash_of_scratch(&self) -> (u32, u64) {
        let mut deg = 0u32;
        let mut hash = 0u64;
        for (i, r) in self.coeffs.iter().enumerate() {
            let e = lane::<L>(&self.scratch, i);
            deg += e;
            hash = hash.wrapping_add(r.wrapping_mul(e as u64));
        }
        (deg, hash)
    }

    fn mask_of_words(&self, words: &[u64]) -> u32 {
        let mut mask = 0u32;
        for i in 0..self.nvars {
            let e = lane::<L>(words, i);
            if e == 0 {
                continue;
            }
            let base = (i % 32) * self.mask_bits;
            for k in 0..self.mask_bits {
                // The ladder is nondecreasing, so the first threshold the
                // exponent fails ends the group.
                if e <= self.thresholds[k] {
                    break;
                }
                mask |= 1 << (base + k);
            }
        }
        mask
    }

    fn intern_scratch(&mut self, deg: u32, hash: u64) -> Result<MonomialId, F4Error> {
        let scratch = std::mem::take(&mut self.scratch);
        let interned = self.intern_words(&scratch, deg, hash);
        self.scratch = scratch;
        interned
    }

    /// Intern packed words that do not alias the exponent slab.
    fn intern_words(&mut self, words: &[u64], deg: u32, hash: u64) -> Result<MonomialId, F4Error> {
        debug_assert_eq!(words.len(), self.words);
        debug_assert_eq!(
            words.iter().fold(0, |acc, w| acc | (w & L::GUARD)),
            0,
            "a stored lane never sets its guard bit"
        );
        let width = self.words;
        let slab = &self.exps;
        if let Some(&found) = self
            .index
            .find(hash, |&i| &slab[i as usize * width..][..width] == words)
        {
            return Ok(MonomialId(found));
        }
        if self.meta.len() >= (u32::MAX - 1) as usize {
            return Err(F4Error::TableFull);
        }
        self.reserve_one()?;
        let id = self.meta.len() as u32;
        let mask = self.mask_of_words(words);
        self.exps.extend_from_slice(words);
        self.meta.push(MonoMeta { deg, mask, hash });
        let meta = &self.meta;
        self.index
            .insert_unique(hash, id, |&i| meta[i as usize].hash);
        Ok(MonomialId(id))
    }
}

/// The capacity one growth step gives a container that holds `capacity`
/// items and must hold `needed`.
///
/// A `Vec` and a `HashTable` both double, so the projection doubles too.
fn grown(capacity: usize, needed: usize) -> usize {
    if needed <= capacity {
        return capacity;
    }
    capacity.saturating_mul(2).max(needed)
}

/// Copy the monomial `id` of `src` into `dst`.
///
/// The mask is recomputed from the exponents, because the two tables may hold
/// different ladders.
pub(crate) fn intern_from<L: Lanes>(
    dst: &mut MonomialTable<L>,
    src: &MonomialTable<L>,
    id: MonomialId,
) -> Result<MonomialId, F4Error> {
    debug_assert_eq!(dst.nvars, src.nvars, "tables must share nvars");
    dst.scratch.copy_from_slice(src.words_at(id.index()));
    let meta = src.meta[id.index()];
    dst.intern_scratch(meta.deg, meta.hash)
}

/// Intern the least common multiple of two monomials of `src` into `dst`.
pub(crate) fn lcm_into<L: Lanes>(
    dst: &mut MonomialTable<L>,
    src: &MonomialTable<L>,
    a: MonomialId,
    b: MonomialId,
) -> Result<MonomialId, F4Error> {
    debug_assert_eq!(dst.nvars, src.nvars, "tables must share nvars");
    let wa = src.words_at(a.index());
    let wb = src.words_at(b.index());
    for k in 0..dst.words {
        dst.scratch[k] = lane_max::<L>(wa[k], wb[k]);
    }
    let (deg, hash) = dst.degree_and_hash_of_scratch();
    dst.intern_scratch(deg, hash)
}

/// Intern `a / b` into `dst`, with `a` in `num` and `b` in `den`.
///
/// The two source tables may hold different ladders, so the divisibility test
/// here uses the degree and the packed words only, never a mask.
pub(crate) fn quotient_into<L: Lanes>(
    dst: &mut MonomialTable<L>,
    num: &MonomialTable<L>,
    a: MonomialId,
    den: &MonomialTable<L>,
    b: MonomialId,
) -> Result<Option<MonomialId>, F4Error> {
    debug_assert_eq!(dst.nvars, num.nvars, "tables must share nvars");
    debug_assert_eq!(dst.nvars, den.nvars, "tables must share nvars");
    let wa = num.words_at(a.index());
    let wb = den.words_at(b.index());
    let (da, db) = (num.meta[a.index()].deg, den.meta[b.index()].deg);
    if db > da || !divides_words::<L>(wb, wa) {
        return Ok(None);
    }
    for k in 0..dst.words {
        dst.scratch[k] = wa[k] - wb[k];
    }
    let hash = num.meta[a.index()]
        .hash
        .wrapping_sub(den.meta[b.index()].hash);
    dst.intern_scratch(da - db, hash).map(Some)
}

/// Report whether the packed monomial `a` divides the packed monomial `b`.
fn divides_words<L: Lanes>(a: &[u64], b: &[u64]) -> bool {
    // The guard bit of every lane of b is clear, so setting it gives each lane
    // a borrow to spend. The bit survives the subtraction exactly when the
    // lane of b is at least the lane of a. Padding lanes are zero in both, so
    // they always survive.
    a.iter()
        .zip(b)
        .all(|(&x, &y)| ((y | L::GUARD) - x) & L::GUARD == L::GUARD)
}

/// The guard bit of every lane that holds a nonzero exponent.
///
/// The guard bit of a stored lane is clear, so setting it gives the lane a
/// borrow to spend. The bit survives exactly when the lane is at least 1.
fn nonzero_lanes<L: Lanes>(w: u64) -> u64 {
    ((w | L::GUARD) - L::ONES) & L::GUARD
}

/// The lanewise maximum of two packed words.
fn lane_max<L: Lanes>(a: u64, b: u64) -> u64 {
    let mask = lane_mask::<L>();
    let mut out = 0u64;
    for slot in 0..L::PER_WORD {
        let shift = slot as u32 * L::BITS;
        let x = (a >> shift) & mask;
        let y = (b >> shift) & mask;
        out |= x.max(y) << shift;
    }
    out
}

fn lane_mask<L: Lanes>() -> u64 {
    (1u64 << L::BITS) - 1
}

fn lane<L: Lanes>(words: &[u64], i: usize) -> u32 {
    let shift = (i % L::PER_WORD) as u32 * L::BITS;
    ((words[i / L::PER_WORD] >> shift) & lane_mask::<L>()) as u32
}

/// Write one lane. The lane must be zero, and `e` must be within the bound.
fn set_lane<L: Lanes>(words: &mut [u64], i: usize, e: u32) {
    let shift = (i % L::PER_WORD) as u32 * L::BITS;
    words[i / L::PER_WORD] |= (e as u64) << shift;
}

fn width_error(e: u32) -> F4Error {
    if e > MAX_EXPONENT {
        F4Error::ExponentOverflow
    } else {
        F4Error::LaneOverflow
    }
}

/// The hash coefficients `r_i`, from the fixed seed.
fn hash_coefficients(nvars: usize) -> Vec<u64> {
    let mut state = HASH_SEED;
    (0..nvars)
        .map(|_| {
            state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            (z ^ (z >> 31)) | 1
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::Monomial;
    use proptest::prelude::*;

    /// The v0.1 monomial of `src/poly.rs`, the reference for every operation.
    fn reference(exps: &[u32]) -> Monomial {
        Monomial::from_exps(exps.iter().map(|&e| e as u16).collect())
    }

    fn clamp<L: Lanes>(exps: &[u32]) -> Vec<u32> {
        exps.iter().map(|&e| e.min(L::BOUND)).collect()
    }

    fn exps_of<L: Lanes>(table: &MonomialTable<L>, id: MonomialId) -> Vec<u32> {
        let mut out = Vec::new();
        table.unpack(id, &mut out);
        out
    }

    /// Compare every table operation on one pair against the reference.
    fn check_pair<L: Lanes>(left: &[u32], right: &[u32]) {
        let a = clamp::<L>(left);
        let b = clamp::<L>(right);
        let (ra, rb) = (reference(&a), reference(&b));
        let mut table = MonomialTable::<L>::new(a.len());
        let ia = table.intern(&a).unwrap();
        let ib = table.intern(&b).unwrap();
        check_basic_operations(&table, ia, ib, &a, &ra, &rb);
        check_mask(&table, ia, ib, &ra, &rb);
        check_quotient(&mut table, ia, ib, &ra, &rb);
        check_lcm(&mut table, ia, ib, &ra, &rb);
        check_product(&mut table, ia, ib, &a, &b, &ra, &rb);
    }

    fn check_basic_operations<L: Lanes>(
        table: &MonomialTable<L>,
        left: MonomialId,
        right: MonomialId,
        exps: &[u32],
        left_ref: &Monomial,
        right_ref: &Monomial,
    ) {
        assert_eq!(table.degree(left), left_ref.deg);
        assert_eq!(table.degree(right), right_ref.deg);
        assert_eq!(exps_of(table, left), exps);
        assert_eq!(table.cmp(left, right), left_ref.cmp(right_ref));
        assert_eq!(table.cmp(right, left), right_ref.cmp(left_ref));
        assert_eq!(table.divides(left, right), left_ref.divides(right_ref));
        assert_eq!(
            table.is_coprime(left, right),
            left_ref.is_coprime(right_ref)
        );
    }

    fn check_mask<L: Lanes>(
        table: &MonomialTable<L>,
        left: MonomialId,
        right: MonomialId,
        left_ref: &Monomial,
        right_ref: &Monomial,
    ) {
        if left_ref.divides(right_ref) {
            assert_eq!(table.mask(left) & !table.mask(right), 0);
        }
    }

    fn check_quotient<L: Lanes>(
        table: &mut MonomialTable<L>,
        left: MonomialId,
        right: MonomialId,
        left_ref: &Monomial,
        right_ref: &Monomial,
    ) {
        match table.quotient(right, left).unwrap() {
            Some(q) => {
                let expected = right_ref.quotient(left_ref).unwrap();
                assert_eq!(exps_of(table, q), clamp::<L>(&widen(&expected)));
                assert_eq!(table.degree(q), expected.deg);
            }
            None => assert!(!left_ref.divides(right_ref)),
        }
    }

    fn check_lcm<L: Lanes>(
        table: &mut MonomialTable<L>,
        left: MonomialId,
        right: MonomialId,
        left_ref: &Monomial,
        right_ref: &Monomial,
    ) {
        let expected = left_ref.lcm(right_ref);
        let actual = table.lcm(left, right).unwrap();
        assert_eq!(exps_of(table, actual), widen(&expected));
        assert_eq!(table.degree(actual), expected.deg);
    }

    #[allow(clippy::too_many_arguments)]
    fn check_product<L: Lanes>(
        table: &mut MonomialTable<L>,
        left: MonomialId,
        right: MonomialId,
        left_exps: &[u32],
        right_exps: &[u32],
        left_ref: &Monomial,
        right_ref: &Monomial,
    ) {
        let fits = left_exps
            .iter()
            .zip(right_exps)
            .all(|(&x, &y)| x + y <= L::BOUND);
        match table.mul(left, right) {
            Ok(product) => {
                assert!(fits);
                assert_eq!(
                    exps_of(table, product),
                    widen(
                        &left_ref
                            .checked_mul(right_ref)
                            .expect("the reference product fits")
                    )
                );
                assert_eq!(table.degree(product), left_ref.deg + right_ref.deg);
            }
            Err(error) => {
                assert!(!fits);
                let worst = left_exps
                    .iter()
                    .zip(right_exps)
                    .map(|(&x, &y)| x + y)
                    .max()
                    .unwrap();
                assert_eq!(error, width_error(worst.min(MAX_EXPONENT + 1)));
            }
        }
    }

    fn widen(monomial: &Monomial) -> Vec<u32> {
        monomial.exps.iter().map(|&e| e as u32).collect()
    }

    fn exponent_vectors() -> impl Strategy<Value = (Vec<u32>, Vec<u32>)> {
        let value = prop_oneof![
            4 => 0u32..=3,
            2 => 0u32..=40,
            1 => 0u32..=65535,
            1 => Just(65535u32),
        ];
        (0usize..=40).prop_flat_map(move |nvars| {
            (
                prop::collection::vec(value.clone(), nvars),
                prop::collection::vec(value.clone(), nvars),
            )
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1024))]

        #[test]
        fn matches_the_reference_at_every_width((a, b) in exponent_vectors()) {
            check_pair::<Lanes8>(&a, &b);
            check_pair::<Lanes16>(&a, &b);
            check_pair::<Lanes32>(&a, &b);
        }
    }

    #[test]
    fn no_variables_gives_one_monomial() {
        let mut table = MonomialTable::<Lanes8>::new(0);
        let one = table.intern(&[]).unwrap();
        assert_eq!(one, MonomialId::ONE);
        assert_eq!(table.len(), 1);
        assert_eq!(table.cmp(one, one), Ordering::Equal);
        assert!(table.divides(one, one));
        assert!(table.is_coprime(one, one));
        assert_eq!(table.mul(one, one).unwrap(), one);
        assert_eq!(table.lcm(one, one).unwrap(), one);
        assert_eq!(table.quotient(one, one).unwrap(), Some(one));
        assert_eq!(table.degree(one), 0);
    }

    #[test]
    fn identity_is_the_first_id() {
        let mut table = MonomialTable::<Lanes8>::new(3);
        assert_eq!(table.intern(&[0, 0, 0]).unwrap(), MonomialId::ONE);
        assert_eq!(table.degree(MonomialId::ONE), 0);
        let x = table.intern(&[1, 0, 0]).unwrap();
        assert_eq!(table.cmp(MonomialId::ONE, x), Ordering::Less);
        assert!(table.divides(MonomialId::ONE, x));
    }

    #[test]
    fn padding_lanes_stay_zero() {
        for nvars in [1usize, 7, 8, 9, 15, 16, 17] {
            let mut table = MonomialTable::<Lanes8>::new(nvars);
            let exps: Vec<u32> = (0..nvars).map(|i| (i % 100) as u32).collect();
            let id = table.intern(&exps).unwrap();
            let words = table.words_of(id);
            for slot in nvars..words.len() * Lanes8::PER_WORD {
                assert_eq!(lane::<Lanes8>(words, slot), 0, "padding lane {slot} set");
            }
        }
    }

    #[test]
    fn word_boundaries_compare_like_the_reference() {
        for nvars in [1usize, 2, 3, 4, 5, 7, 8, 9, 16, 17] {
            let mut table = MonomialTable::<Lanes8>::new(nvars);
            let base: Vec<u32> = (0..nvars).map(|i| (i % 5) as u32 + 1).collect();
            let id = table.intern(&base).unwrap();
            for var in 0..nvars {
                let mut other = base.clone();
                other[var] += 1;
                let up = table.intern(&other).unwrap();
                assert_eq!(table.cmp(id, up), reference(&base).cmp(&reference(&other)));
                assert!(table.divides(id, up));
                assert!(!table.divides(up, id));
                let mut down = base.clone();
                down[var] -= 1;
                let low = table.intern(&down).unwrap();
                assert_eq!(table.cmp(id, low), reference(&base).cmp(&reference(&down)));
                assert!(table.divides(low, id));
                assert!(!table.divides(id, low));
            }
        }
    }

    #[test]
    fn hash_consing_returns_one_id() {
        let mut table = MonomialTable::<Lanes16>::new(4);
        let first = table.intern(&[3, 0, 7, 1]).unwrap();
        let again = table.intern(&[3, 0, 7, 1]).unwrap();
        assert_eq!(first, again);
        assert_eq!(table.len(), 2);
        let other = table.intern(&[3, 0, 7, 2]).unwrap();
        assert_ne!(first, other);
        assert_eq!(table.len(), 3);
        let product = table.mul(first, MonomialId::ONE).unwrap();
        assert_eq!(product, first);
        assert_eq!(table.len(), 3);
    }

    #[test]
    fn ids_repeat_for_the_same_insertion_sequence() {
        let sequence: Vec<Vec<u32>> = (0..200)
            .map(|k: u32| vec![k % 7, (k * 3) % 5, k % 2, (k * 11) % 13])
            .collect();
        let build = || {
            let mut table = MonomialTable::<Lanes8>::new(4);
            let ids: Vec<MonomialId> = sequence
                .iter()
                .map(|exps| table.intern(exps).unwrap())
                .collect();
            (table.len(), ids)
        };
        let (len_a, ids_a) = build();
        let (len_b, ids_b) = build();
        assert_eq!(len_a, len_b);
        assert_eq!(ids_a, ids_b);
    }

    #[test]
    fn clear_keeps_the_identity_only() {
        let mut table = MonomialTable::<Lanes8>::new(3);
        table.intern(&[1, 2, 3]).unwrap();
        table.intern(&[4, 0, 0]).unwrap();
        assert_eq!(table.len(), 3);
        table.clear();
        assert_eq!(table.len(), 1);
        assert_eq!(exps_of(&table, MonomialId::ONE), vec![0, 0, 0]);
        assert_eq!(table.intern(&[4, 0, 0]).unwrap().index(), 1);
    }

    #[test]
    fn an_exponent_past_the_lane_bound_restarts_the_run() {
        let mut table = MonomialTable::<Lanes8>::new(2);
        assert_eq!(table.intern(&[128, 0]), Err(F4Error::LaneOverflow));
        let high = table.intern(&[127, 1]).unwrap();
        let one_more = table.intern(&[1, 0]).unwrap();
        assert_eq!(table.mul(high, one_more), Err(F4Error::LaneOverflow));

        let mut wider = MonomialTable::<Lanes16>::new(2);
        let same = wider.intern(&[128, 0]).unwrap();
        assert_eq!(wider.degree(same), 128);
        assert_eq!(wider.intern(&[32768, 0]), Err(F4Error::LaneOverflow));
    }

    #[test]
    fn an_exponent_past_the_contract_has_no_width() {
        let mut table = MonomialTable::<Lanes32>::new(2);
        assert_eq!(table.intern(&[65536, 0]), Err(F4Error::ExponentOverflow));
        let high = table.intern(&[65535, 0]).unwrap();
        let one = table.intern(&[1, 0]).unwrap();
        assert_eq!(table.mul(high, one), Err(F4Error::ExponentOverflow));
        assert_eq!(
            table.lcm(high, one).map(|id| exps_of(&table, id)),
            Ok(vec![65535, 0])
        );
    }

    #[test]
    fn cross_store_operations_match_one_table() {
        let mut basis = MonomialTable::<Lanes8>::with_max_exponent(4, 12);
        let mut symbolic = MonomialTable::<Lanes8>::with_max_exponent(4, 60);
        assert_ne!(basis.thresholds, symbolic.thresholds);

        let a = basis.intern(&[3, 1, 0, 5]).unwrap();
        let b = basis.intern(&[1, 4, 0, 2]).unwrap();
        check_cross_copy(&basis, &mut symbolic, a);
        check_cross_lcm(&mut basis, &mut symbolic, a, b);
        check_cross_quotient(&mut basis, &mut symbolic, b);
        check_cross_product(&basis, &mut symbolic, b);
    }

    fn check_cross_copy(
        basis: &MonomialTable<Lanes8>,
        symbolic: &mut MonomialTable<Lanes8>,
        mono: MonomialId,
    ) {
        let moved = intern_from(symbolic, basis, mono).unwrap();
        assert_eq!(exps_of(symbolic, moved), vec![3, 1, 0, 5]);
        assert_eq!(symbolic.degree(moved), basis.degree(mono));
        assert_eq!(symbolic.hash(moved), basis.hash(mono));
        assert_eq!(symbolic.mask(moved), mask_here(symbolic, moved));
        assert_ne!(symbolic.mask(moved), basis.mask(mono));
    }

    fn check_cross_lcm(
        basis: &mut MonomialTable<Lanes8>,
        symbolic: &mut MonomialTable<Lanes8>,
        left: MonomialId,
        right: MonomialId,
    ) {
        let inside = basis.lcm(left, right).unwrap();
        let across = lcm_into(symbolic, basis, left, right).unwrap();
        assert_eq!(exps_of(symbolic, across), exps_of(basis, inside));
    }

    fn check_cross_quotient(
        basis: &mut MonomialTable<Lanes8>,
        symbolic: &mut MonomialTable<Lanes8>,
        divisor: MonomialId,
    ) {
        let left = basis.intern(&[3, 1, 0, 5]).unwrap();
        let inside = basis.quotient(left, divisor).unwrap();
        let across = quotient_into(symbolic, basis, left, basis, divisor).unwrap();
        assert_eq!(inside.is_none(), across.is_none());
        let numerator = basis.intern(&[4, 5, 0, 7]).unwrap();
        let inside = basis.quotient(numerator, divisor).unwrap().unwrap();
        let across = quotient_into(symbolic, basis, numerator, basis, divisor)
            .unwrap()
            .unwrap();
        assert_eq!(exps_of(symbolic, across), exps_of(basis, inside));
    }

    fn check_cross_product(
        basis: &MonomialTable<Lanes8>,
        symbolic: &mut MonomialTable<Lanes8>,
        mono: MonomialId,
    ) {
        let multiplier = symbolic.intern(&[2, 0, 1, 0]).unwrap();
        let product = symbolic.mul_external(multiplier, basis, mono).unwrap();
        assert_eq!(exps_of(symbolic, product), vec![3, 4, 1, 2]);
        assert_eq!(
            symbolic.hash(product),
            symbolic.hash(multiplier).wrapping_add(basis.hash(mono))
        );
    }

    /// A deterministic pool of monomials of degree at most `max_degree`.
    fn monomial_pool(nvars: usize, max_degree: u32, count: usize) -> Vec<Vec<u32>> {
        let mut state = 0x243f_6a88_85a3_08d3u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        (0..count)
            .map(|_| {
                let mut exps = vec![0u32; nvars];
                let degree = (next() % (max_degree as u64 + 1)) as u32;
                for _ in 0..degree {
                    exps[(next() % nvars as u64) as usize] += 1;
                }
                exps
            })
            .collect()
    }

    #[test]
    #[ignore = "timing study: cargo test --release --lib -- --ignored --nocapture"]
    fn intern_and_compare_throughput() {
        use std::time::Instant;

        for (nvars, max_degree) in [(9usize, 12u32), (7, 10)] {
            let pool = monomial_pool(nvars, max_degree, 60_000);
            let rounds = 20;
            let mut bytes = 0;
            let mut sink = 0usize;

            let start = Instant::now();
            for _ in 0..rounds {
                let mut table = MonomialTable::<Lanes8>::with_max_exponent(nvars, max_degree);
                for exps in &pool {
                    sink += table.intern(exps).unwrap().index();
                }
                bytes = table.memory_bytes();
            }
            let intern_ns = start.elapsed().as_nanos() as f64 / (rounds * pool.len()) as f64;

            let mut table = MonomialTable::<Lanes8>::with_max_exponent(nvars, max_degree);
            let ids: Vec<MonomialId> = pool
                .iter()
                .map(|exps| table.intern(exps).unwrap())
                .collect();
            let start = Instant::now();
            let mut comparisons = 0usize;
            for _ in 0..rounds {
                let mut order = ids.clone();
                order.sort_by(|&x, &y| {
                    comparisons += 1;
                    table.cmp(x, y)
                });
                sink += order[0].index();
            }
            let sort_ns = start.elapsed().as_nanos() as f64 / comparisons as f64;

            println!(
                "nvars {nvars} degree {max_degree}: {} monomials, intern {intern_ns:.1} ns, \
                 compare {sort_ns:.1} ns, table {bytes} bytes (sink {sink})",
                table.len()
            );
        }
    }

    fn mask_here<L: Lanes>(table: &MonomialTable<L>, id: MonomialId) -> u32 {
        table.mask_of_words(table.words_of(id))
    }

    #[test]
    fn mul_external_survives_a_slab_reallocation() {
        let mut basis = MonomialTable::<Lanes8>::new(5);
        let mut symbolic = MonomialTable::<Lanes8>::new(5);
        let rhs = basis.intern(&[0, 1, 0, 0, 0]).unwrap();
        let mut expected = Vec::new();
        for k in 0..1000u32 {
            let lhs = symbolic.intern(&[k % 11, 0, k % 7, k % 3, k % 5]).unwrap();
            let product = symbolic.mul_external(lhs, &basis, rhs).unwrap();
            let mut want = exps_of(&symbolic, lhs);
            want[1] += 1;
            assert_eq!(exps_of(&symbolic, product), want);
            expected.push((product, want));
        }
        for (id, want) in expected {
            assert_eq!(exps_of(&symbolic, id), want);
        }
    }

    #[test]
    fn the_mask_never_rejects_a_true_divisor() {
        let mut table = MonomialTable::<Lanes8>::with_max_exponent(6, 9);
        let mut ids = Vec::new();
        for k in 0..400u32 {
            let exps = vec![k % 10, k % 4, (k * 3) % 7, k % 2, (k * 5) % 11, k % 3];
            ids.push((table.intern(&exps).unwrap(), exps));
        }
        for (ia, ea) in &ids {
            for (ib, eb) in &ids {
                let divides = reference(ea).divides(&reference(eb));
                assert_eq!(table.divides(*ia, *ib), divides);
                if divides {
                    assert_eq!(table.mask(*ia) & !table.mask(*ib), 0);
                }
            }
        }
    }
}
