//! Monomials, terms, and polynomials of a [`PolynomialRing`].

use smallvec::SmallVec;
use std::cmp::Ordering;
use std::fmt;
use std::mem::size_of;

use crate::ring::PolynomialRing;
use crate::ring::field::Felt;

/// The exponent vector of a monomial.
///
/// Past [`INLINE_VARIABLES`] variables, the vector spills to the heap.
pub(crate) type Exps = SmallVec<[u16; INLINE_VARIABLES]>;

/// The number of exponents a monomial holds inline.
///
/// Eleven exponents are the most that fit the space `SmallVec` reserves
/// for the inline array, so `Monomial` stays 40 bytes and `Term` 48. A
/// wider inline array grows both and measures slower on every benchmark
/// input.
const INLINE_VARIABLES: usize = 11;

/// The heap bytes one monomial's exponent vector costs beyond what
/// `size_of::<Term>()` already counts.
///
/// A ring of at most [`INLINE_VARIABLES`] variables holds every exponent
/// vector inline, so the heap cost is zero. Past that width, `Exps` spills
/// to the heap and allocates space for all `nvars` exponents there, so a
/// memory count that only reads `size_of::<Term>()` misses this every
/// time.
pub(crate) fn heap_exps_bytes(nvars: usize) -> usize {
    if nvars > INLINE_VARIABLES {
        nvars * size_of::<u16>()
    } else {
        0
    }
}

/// A monomial product past the width one exponent holds.
///
/// [`Monomial::checked_mul`] reports this. The engines turn it into
/// [`crate::compute::ComputeError::DegreeLimit`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ExponentOverflow;

/// A monomial as an exponent vector with its total degree.
///
/// Every monomial of one ring holds one exponent per variable. The
/// arithmetic checks that width in debug builds, because a mismatch is a
/// defect.
#[derive(Debug, PartialEq, Eq, Hash)]
pub(crate) struct Monomial {
    pub(crate) exps: Exps,
    pub(crate) deg: u32,
}

/// Copy an exponent vector.
///
/// Every monomial the engines build starts from a copy of one they hold,
/// so this is the innermost step of both backends. `SmallVec`'s `Clone`
/// runs a capacity check per element, and `from_slice` calls `memcpy`.
/// A ring inside the inline width copies a fixed-width buffer instead,
/// which the compiler unrolls.
#[inline]
fn copy_exps(exps: &[u16]) -> Exps {
    if exps.len() > INLINE_VARIABLES {
        return Exps::from_slice(exps);
    }
    let mut buf = [0u16; INLINE_VARIABLES];
    for (slot, &exp) in buf.iter_mut().zip(exps) {
        *slot = exp;
    }
    Exps::from_buf_and_len(buf, exps.len())
}

/// Build an exponent vector from two of the same width, entry by entry.
///
/// This is [`copy_exps`] with one operation folded into the copy, so a
/// monomial product or quotient passes over the exponents once.
#[inline]
fn combine_exps(a: &[u16], b: &[u16], op: impl Fn(u16, u16) -> u16) -> Exps {
    debug_assert_eq!(a.len(), b.len(), "monomials must share nvars");
    if a.len() > INLINE_VARIABLES {
        return a.iter().zip(b).map(|(&x, &y)| op(x, y)).collect();
    }
    let mut buf = [0u16; INLINE_VARIABLES];
    for (slot, (&x, &y)) in buf.iter_mut().zip(a.iter().zip(b)) {
        *slot = op(x, y);
    }
    Exps::from_buf_and_len(buf, a.len())
}

impl Clone for Monomial {
    #[inline]
    fn clone(&self) -> Self {
        Monomial {
            exps: copy_exps(&self.exps),
            deg: self.deg,
        }
    }
}

impl Monomial {
    pub(crate) fn one(nvars: usize) -> Self {
        Monomial {
            exps: Exps::from_elem(0u16, nvars),
            deg: 0,
        }
    }

    pub(crate) fn from_exps(exps: Exps) -> Self {
        let deg = exps
            .iter()
            .try_fold(0u32, |acc, &e| acc.checked_add(e as u32))
            // A ring holds at most 256 variables, each below 2^16.
            .expect("monomial degree overflow");
        Monomial { exps, deg }
    }

    pub(crate) fn nvars(&self) -> usize {
        self.exps.len()
    }

    pub(crate) fn divides(&self, other: &Self) -> bool {
        debug_assert_eq!(self.nvars(), other.nvars(), "monomials must share nvars");
        // A divisor has the smaller degree, and most candidates fail here.
        self.deg <= other.deg && self.exps.iter().zip(&other.exps).all(|(a, b)| a <= b)
    }

    /// Return `self * other`, or [`ExponentOverflow`] when an exponent of
    /// the product leaves the width a `u16` holds.
    ///
    /// The degree of the product is the sum of the two degrees, which fits
    /// a `u32` for the same reason [`Monomial::from_exps`] gives.
    pub(crate) fn checked_mul(&self, other: &Self) -> Result<Self, ExponentOverflow> {
        debug_assert_eq!(self.nvars(), other.nvars(), "monomials must share nvars");
        let deg = self.deg + other.deg;
        // Every exponent of the product is at most the product's degree, so
        // a degree inside a u16 rules out an exponent overflow and the
        // addition needs no per-entry check. Past that width one entry may
        // overflow, so that path checks each entry and reports it.
        if deg <= u16::MAX as u32 {
            Ok(Monomial {
                exps: combine_exps(&self.exps, &other.exps, u16::wrapping_add),
                deg,
            })
        } else {
            let mut exps = copy_exps(&self.exps);
            for (slot, &b) in exps.iter_mut().zip(&other.exps) {
                *slot = slot.checked_add(b).ok_or(ExponentOverflow)?;
            }
            Ok(Monomial { exps, deg })
        }
    }

    pub(crate) fn quotient(&self, divisor: &Self) -> Option<Self> {
        debug_assert_eq!(self.nvars(), divisor.nvars(), "monomials must share nvars");
        if !divisor.divides(self) {
            return None;
        }
        Some(Monomial {
            exps: combine_exps(&self.exps, &divisor.exps, |a, b| a - b),
            deg: self.deg - divisor.deg,
        })
    }

    pub(crate) fn lcm(&self, other: &Self) -> Self {
        debug_assert_eq!(self.nvars(), other.nvars(), "monomials must share nvars");
        Monomial::from_exps(combine_exps(&self.exps, &other.exps, u16::max))
    }

    /// Compare `self * (num / den)` with `other` without building the product.
    ///
    /// `den` must divide `num`. The result is what [`Ord`] gives for the
    /// product, so a caller that rejects the product pays no allocation.
    pub(crate) fn cmp_shifted(&self, num: &Self, den: &Self, other: &Self) -> Ordering {
        debug_assert_eq!(self.nvars(), num.nvars(), "monomials must share nvars");
        debug_assert_eq!(self.nvars(), den.nvars(), "monomials must share nvars");
        debug_assert_eq!(self.nvars(), other.nvars(), "monomials must share nvars");
        match (self.deg + num.deg - den.deg).cmp(&other.deg) {
            Ordering::Equal => {
                for index in (0..self.exps.len()).rev() {
                    let exp =
                        self.exps[index] as u32 + num.exps[index] as u32 - den.exps[index] as u32;
                    match exp.cmp(&(other.exps[index] as u32)) {
                        Ordering::Equal => continue,
                        ord => return ord.reverse(),
                    }
                }
                Ordering::Equal
            }
            ord => ord,
        }
    }

    /// The first sixteen exponents in four bits each, capped at seven.
    ///
    /// A divisor's exponents all sit at or below its multiple's, and the
    /// cap keeps that true, so [`key_divides`] rejects a candidate without
    /// reading either exponent vector. Variables past the sixteenth hold
    /// no field, which only makes the test accept more.
    pub(crate) fn divisor_key(&self) -> u64 {
        let mut key = 0u64;
        for (index, &exp) in self.exps.iter().take(16).enumerate() {
            key |= (exp.min(7) as u64) << (index * 4);
        }
        key
    }

    /// Report whether no variable appears in both monomials.
    pub(crate) fn is_coprime(&self, other: &Self) -> bool {
        self.exps
            .iter()
            .zip(&other.exps)
            .all(|(&a, &b)| a == 0 || b == 0)
    }
}

/// Report whether every field of `divisor` is at or below `multiple`.
///
/// Both keys come from [`Monomial::divisor_key`], so every field is at
/// most seven and the borrow of one field never reaches the next.
#[inline]
pub(crate) fn key_divides(divisor: u64, multiple: u64) -> bool {
    const HIGH: u64 = 0x8888_8888_8888_8888;
    ((multiple | HIGH) - divisor) & HIGH == HIGH
}

/// Four exponents as one integer, the later variable in the higher bits.
///
/// Comparing two of these as integers gives the reverse-lexicographic
/// order of the four exponents they hold.
#[inline]
fn pack_four(exps: &[u16]) -> u64 {
    let [a, b, c, d]: [u16; 4] = exps.try_into().expect("a chunk holds four exponents");
    a as u64 | (b as u64) << 16 | (c as u64) << 32 | (d as u64) << 48
}

/// Compare two exponent vectors of one width by the last entry that
/// differs, with the smaller entry the greater monomial.
///
/// This is the tie-break of grevlex, and the merge inside every reduction
/// step runs it once per term, so it takes four exponents per comparison
/// where the width allows.
#[inline]
fn reverse_lex(a: &[u16], b: &[u16]) -> Ordering {
    let mut a_chunks = a.rchunks_exact(4);
    let mut b_chunks = b.rchunks_exact(4);
    for (a_chunk, b_chunk) in a_chunks.by_ref().zip(b_chunks.by_ref()) {
        let (a_word, b_word) = (pack_four(a_chunk), pack_four(b_chunk));
        if a_word != b_word {
            return a_word.cmp(&b_word).reverse();
        }
    }
    for (&a, &b) in a_chunks.remainder().iter().zip(b_chunks.remainder()).rev() {
        match a.cmp(&b) {
            Ordering::Equal => continue,
            ord => return ord.reverse(),
        }
    }
    Ordering::Equal
}

/// Grevlex, with the larger monomial the greater value.
impl Ord for Monomial {
    fn cmp(&self, other: &Self) -> Ordering {
        debug_assert_eq!(self.nvars(), other.nvars(), "monomials must share nvars");
        match self.deg.cmp(&other.deg) {
            Ordering::Equal => reverse_lex(&self.exps, &other.exps),
            ord => ord,
        }
    }
}

impl PartialOrd for Monomial {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// One coefficient and one monomial.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Term {
    pub(crate) coeff: Felt,
    pub(crate) mono: Monomial,
}

/// A polynomial of a [`PolynomialRing`].
///
/// Build one with [`PolynomialRing::polynomial`] or
/// [`PolynomialRing::parse_polynomial`]. A polynomial carries its ring, so
/// the zero polynomial names its variables and its prime. The coefficients
/// are reduced, the monomials are distinct, and no coefficient is zero.
///
/// `Display` writes the syntax [`PolynomialRing::parse_polynomial`] reads,
/// with the terms in descending grevlex order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial {
    ring: PolynomialRing,
    /// Terms in ascending grevlex order, so the last term is the leading
    /// term.
    pub(crate) terms: Vec<Term>,
}

impl Polynomial {
    /// The ring this polynomial belongs to.
    pub fn ring(&self) -> &PolynomialRing {
        &self.ring
    }

    /// Report whether the polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// The terms, largest monomial first.
    ///
    /// Each item is the coefficient, a representative in `[1, p)`, and one
    /// exponent per variable.
    pub fn terms(&self) -> impl ExactSizeIterator<Item = (u64, &[u16])> {
        self.terms
            .iter()
            .rev()
            .map(|term| (term.coeff.value(), term.mono.exps.as_slice()))
    }

    /// The largest term, or `None` for the zero polynomial.
    pub fn leading_term(&self) -> Option<(u64, &[u16])> {
        self.lt()
            .map(|term| (term.coeff.value(), term.mono.exps.as_slice()))
    }

    /// The total degree, or `None` for the zero polynomial.
    ///
    /// Grevlex is graded, so the leading monomial carries the degree.
    pub fn degree(&self) -> Option<u32> {
        self.lt().map(|term| term.mono.deg)
    }

    /// Build a polynomial from terms that already hold the invariants.
    ///
    /// The terms run strictly ascending under grevlex, no coefficient is
    /// zero, and every monomial holds one exponent per variable.
    pub(crate) fn from_sorted_terms(ring: PolynomialRing, terms: Vec<Term>) -> Self {
        debug_assert!(
            terms.iter().all(|term| {
                !term.coeff.is_zero()
                    && term.coeff.value() < ring.modulus()
                    && term.mono.nvars() == ring.nvars()
            }) && terms.windows(2).all(|pair| pair[0].mono < pair[1].mono),
            "a polynomial must hold reduced coefficients and ascending distinct monomials"
        );
        Polynomial { ring, terms }
    }

    /// Build a polynomial from unsorted terms, merging like monomials.
    ///
    /// A term that reduces to zero drops out.
    pub(crate) fn from_terms(ring: PolynomialRing, mut terms: Vec<Term>) -> Self {
        let p = ring.modulus();
        terms.retain(|t| !t.coeff.is_zero());
        terms.sort_by(|a, b| a.mono.cmp(&b.mono));

        let mut out: Vec<Term> = Vec::with_capacity(terms.len());
        for term in terms {
            if let Some(last) = out.last_mut()
                && last.mono == term.mono
            {
                let c = last.coeff.add(term.coeff, p);
                if c.is_zero() {
                    out.pop();
                } else {
                    last.coeff = c;
                }
                continue;
            }
            out.push(term);
        }
        Polynomial { ring, terms: out }
    }

    /// The zero polynomial of the same ring.
    pub(crate) fn zero_like(&self) -> Self {
        Polynomial {
            ring: self.ring.clone(),
            terms: Vec::new(),
        }
    }

    /// Stop a step that mixes two rings, or that names the wrong prime.
    ///
    /// The arithmetic takes the modulus as an argument, so this is what
    /// ties the argument back to the ring the result claims.
    fn check_operand(&self, other: &Self, p: u64) {
        debug_assert_eq!(p, self.ring.modulus(), "the modulus must be the ring's");
        debug_assert!(self.ring == other.ring, "both operands must share a ring");
    }

    fn with_terms(&self, terms: Vec<Term>) -> Self {
        Polynomial {
            ring: self.ring.clone(),
            terms,
        }
    }

    pub(crate) fn lt(&self) -> Option<&Term> {
        self.terms.last()
    }

    pub(crate) fn lc(&self) -> Option<Felt> {
        self.lt().map(|t| t.coeff)
    }

    pub(crate) fn lm(&self) -> Option<&Monomial> {
        self.lt().map(|t| &t.mono)
    }

    pub(crate) fn pop_lt(&mut self) -> Option<Term> {
        self.terms.pop()
    }

    pub(crate) fn push_term(&mut self, term: Term, p: u64) {
        if term.coeff.is_zero() {
            return;
        }
        match self
            .terms
            .binary_search_by(|probe| probe.mono.cmp(&term.mono))
        {
            Ok(pos) => {
                let new_coeff = self.terms[pos].coeff.add(term.coeff, p);
                if new_coeff.is_zero() {
                    self.terms.remove(pos);
                } else {
                    self.terms[pos].coeff = new_coeff;
                }
            }
            Err(pos) => self.terms.insert(pos, term),
        }
    }

    pub(crate) fn add(&self, other: &Self, p: u64) -> Self {
        self.merge::<false>(other, p)
    }

    pub(crate) fn sub(&self, other: &Self, p: u64) -> Self {
        self.merge::<true>(other, p)
    }

    /// Merge two polynomials of one ring under grevlex.
    ///
    /// `SUB` negates every term of `other`, so one loop serves both
    /// `self + other` and `self - other`. The choice is a `const`, so the
    /// compiler resolves it and the loop carries no runtime branch for it.
    fn merge<const SUB: bool>(&self, other: &Self, p: u64) -> Self {
        self.check_operand(other, p);
        let mut out: Vec<Term> = Vec::with_capacity(self.terms.len() + other.terms.len());

        let mut i = 0usize;
        let mut j = 0usize;
        while i < self.terms.len() && j < other.terms.len() {
            let a = &self.terms[i];
            let b = &other.terms[j];
            match a.mono.cmp(&b.mono) {
                Ordering::Less => {
                    out.push(a.clone());
                    i += 1;
                }
                Ordering::Greater => {
                    out.push(Term {
                        coeff: if SUB { b.coeff.neg(p) } else { b.coeff },
                        mono: b.mono.clone(),
                    });
                    j += 1;
                }
                Ordering::Equal => {
                    let c = if SUB {
                        a.coeff.sub(b.coeff, p)
                    } else {
                        a.coeff.add(b.coeff, p)
                    };
                    if !c.is_zero() {
                        out.push(Term {
                            coeff: c,
                            mono: a.mono.clone(),
                        });
                    }
                    i += 1;
                    j += 1;
                }
            }
        }

        out.extend(self.terms[i..].iter().cloned());
        out.extend(other.terms[j..].iter().map(|term| Term {
            coeff: if SUB { term.coeff.neg(p) } else { term.coeff },
            mono: term.mono.clone(),
        }));
        self.with_terms(out)
    }

    /// Return `c * m * self`, or [`ExponentOverflow`] when a product leaves
    /// the width one exponent holds.
    pub(crate) fn scale_monomial(
        &self,
        c: Felt,
        m: &Monomial,
        p: u64,
    ) -> Result<Self, ExponentOverflow> {
        if c.is_zero() || self.is_zero() {
            return Ok(self.zero_like());
        }
        if c == Felt::one() {
            return self.shift_monomial(m);
        }
        let mut terms: Vec<Term> = Vec::with_capacity(self.terms.len());
        for t in &self.terms {
            terms.push(Term {
                coeff: t.coeff.mul(c, p),
                mono: t.mono.checked_mul(m)?,
            });
        }
        Ok(self.with_terms(terms))
    }

    /// Return `m * self`, or [`ExponentOverflow`] when a product leaves the
    /// width one exponent holds.
    ///
    /// Grevlex is a monomial order, so the terms stay in ascending order.
    pub(crate) fn shift_monomial(&self, m: &Monomial) -> Result<Self, ExponentOverflow> {
        let mut terms: Vec<Term> = Vec::with_capacity(self.terms.len());
        for t in &self.terms {
            terms.push(Term {
                coeff: t.coeff,
                mono: t.mono.checked_mul(m)?,
            });
        }
        Ok(self.with_terms(terms))
    }

    /// Return `self - c * m * other`, or [`ExponentOverflow`] when a
    /// product leaves the width one exponent holds.
    pub(crate) fn sub_scaled(
        &self,
        other: &Self,
        c: Felt,
        m: &Monomial,
        p: u64,
    ) -> Result<Self, ExponentOverflow> {
        self.check_operand(other, p);
        if other.is_zero() || c.is_zero() {
            return Ok(self.clone());
        }

        let mut out: Vec<Term> = Vec::with_capacity(self.terms.len() + other.terms.len());

        let mut i = 0usize;
        let mut j = 0usize;
        // The multiple of the term at `j`, held across the iterations that
        // advance `i` only.
        let mut b_mono: Option<Monomial> = None;

        while i < self.terms.len() && j < other.terms.len() {
            let a = &self.terms[i];
            let b0 = &other.terms[j];
            if b_mono.is_none() {
                b_mono = Some(b0.mono.checked_mul(m)?);
            }
            let ord = a.mono.cmp(b_mono.as_ref().expect("the multiple is held"));
            match ord {
                Ordering::Less => {
                    out.push(a.clone());
                    i += 1;
                }
                Ordering::Greater => {
                    let b_coeff = b0.coeff.mul(c, p).neg(p);
                    let mono = b_mono.take().expect("the multiple is held");
                    out.push(Term {
                        coeff: b_coeff,
                        mono,
                    });
                    j += 1;
                }
                Ordering::Equal => {
                    let b_coeff = b0.coeff.mul(c, p);
                    let new_coeff = a.coeff.sub(b_coeff, p);
                    let mono = b_mono.take().expect("the multiple is held");
                    if !new_coeff.is_zero() {
                        out.push(Term {
                            coeff: new_coeff,
                            mono,
                        });
                    }
                    i += 1;
                    j += 1;
                }
            }
        }

        out.extend(self.terms[i..].iter().cloned());
        append_scaled_tail(&mut out, &other.terms, j, b_mono, c, m, p)?;

        Ok(self.with_terms(out))
    }

    pub(crate) fn make_monic(&self, p: u64) -> Self {
        let Some(lc) = self.lc() else {
            return self.clone();
        };
        if lc == Felt::one() {
            return self.clone();
        }
        let inv = lc.inv(p);
        let terms: Vec<Term> = self
            .terms
            .iter()
            .map(|t| Term {
                coeff: t.coeff.mul(inv, p),
                mono: t.mono.clone(),
            })
            .collect();
        self.with_terms(terms)
    }

    #[cfg(test)]
    pub(crate) fn normal_form(
        &self,
        reducers: &[Polynomial],
        p: u64,
    ) -> Result<Self, ExponentOverflow> {
        let reducers: Vec<&Polynomial> = reducers.iter().collect();
        self.normal_form_refs(&reducers, p)
    }

    /// Reduce `self` by `reducers` and return the remainder, or
    /// [`ExponentOverflow`] when a reduction multiple leaves the width one
    /// exponent holds.
    pub(crate) fn normal_form_refs(
        &self,
        reducers: &[&Polynomial],
        p: u64,
    ) -> Result<Self, ExponentOverflow> {
        let mut poly = self.clone();
        // The remainder grows by the largest term left, so it is built
        // descending and turned around once.
        let mut remainder: Vec<Term> = Vec::new();

        while let Some(lt_p) = poly.lt().cloned() {
            let mut reduced = false;

            for g in reducers {
                let Some(lt_g) = g.lt() else { continue };
                if lt_g.mono.divides(&lt_p.mono) {
                    let m = lt_p
                        .mono
                        .quotient(&lt_g.mono)
                        .expect("divides() implies quotient()");
                    let scale = lt_p.coeff.div(lt_g.coeff, p);
                    poly = poly.sub_scaled(g, scale, &m, p)?;
                    reduced = true;
                    break;
                }
            }

            if !reduced {
                let term = poly.pop_lt().expect("lt_p was Some, so poly is non-empty");
                remainder.push(term);
            }
        }

        remainder.reverse();
        Ok(Polynomial::from_sorted_terms(self.ring.clone(), remainder))
    }

    /// Return the S-polynomial of `self` and `other`, or
    /// [`ExponentOverflow`] when a multiplier leaves the width one exponent
    /// holds.
    pub(crate) fn s_polynomial(&self, other: &Self, p: u64) -> Result<Self, ExponentOverflow> {
        self.check_operand(other, p);
        if self.is_zero() || other.is_zero() {
            return Ok(self.zero_like());
        }
        let lt_f = self.lt().expect("non-zero polynomial must have lt");
        let lt_g = other.lt().expect("non-zero polynomial must have lt");
        let lcm = lt_f.mono.lcm(&lt_g.mono);

        let m_f = lcm.quotient(&lt_f.mono).expect("lm(f) divides lcm");
        let m_g = lcm.quotient(&lt_g.mono).expect("lm(g) divides lcm");

        let f_scaled = self.scale_monomial(lt_g.coeff, &m_f, p)?;
        let g_scaled = other.scale_monomial(lt_f.coeff, &m_g, p)?;
        Ok(f_scaled.sub(&g_scaled, p))
    }
}

fn append_scaled_tail(
    out: &mut Vec<Term>,
    terms: &[Term],
    mut start: usize,
    held: Option<Monomial>,
    c: Felt,
    m: &Monomial,
    p: u64,
) -> Result<(), ExponentOverflow> {
    if let Some(mono) = held {
        push_negated_scaled(out, &terms[start], c, mono, p);
        start += 1;
    }
    for term in &terms[start..] {
        push_negated_scaled(out, term, c, term.mono.checked_mul(m)?, p);
    }
    Ok(())
}

fn push_negated_scaled(out: &mut Vec<Term>, term: &Term, c: Felt, mono: Monomial, p: u64) {
    out.push(Term {
        coeff: term.coeff.mul(c, p).neg(p),
        mono,
    });
}

fn write_term(
    f: &mut fmt::Formatter<'_>,
    variables: &[String],
    coeff: u64,
    exps: &[u16],
) -> fmt::Result {
    let mut written = false;
    if coeff != 1 || exps.iter().all(|&exp| exp == 0) {
        write!(f, "{coeff}")?;
        written = true;
    }
    for (name, &exp) in variables.iter().zip(exps) {
        written = write_variable(f, name, exp, written)?;
    }
    Ok(())
}

fn write_variable(
    f: &mut fmt::Formatter<'_>,
    name: &str,
    exp: u16,
    mut written: bool,
) -> Result<bool, fmt::Error> {
    if exp == 0 {
        return Ok(written);
    }
    if written {
        f.write_str("*")?;
    }
    written = true;
    match exp {
        1 => f.write_str(name)?,
        _ => write!(f, "{name}^{exp}")?,
    }
    Ok(written)
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return f.write_str("0");
        }
        let variables = self.ring.variables();
        for (index, (coeff, exps)) in self.terms().enumerate() {
            if index > 0 {
                f.write_str(" + ")?;
            }
            write_term(f, variables, coeff, exps)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::PolynomialRing;

    fn ring() -> PolynomialRing {
        PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime")
    }

    fn mono(exps: &[u16]) -> Monomial {
        Monomial::from_exps(exps.iter().copied().collect())
    }

    #[test]
    fn grevlex_orders_by_degree_then_reverse_lex() {
        assert!(mono(&[3, 0]) > mono(&[0, 2]));
        assert!(mono(&[2, 0]) > mono(&[1, 1]));
        assert!(mono(&[1, 1]) > mono(&[0, 2]));
        assert!(mono(&[2, 1]) > mono(&[1, 2]));
    }

    #[test]
    fn packed_comparison_agrees_with_the_entry_by_entry_order() {
        fn by_entry(a: &Monomial, b: &Monomial) -> Ordering {
            match a.deg.cmp(&b.deg) {
                Ordering::Equal => {
                    for (&x, &y) in a.exps.iter().zip(&b.exps).rev() {
                        match x.cmp(&y) {
                            Ordering::Equal => continue,
                            ord => return ord.reverse(),
                        }
                    }
                    Ordering::Equal
                }
                ord => ord,
            }
        }

        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        for nvars in 1..=13usize {
            for _ in 0..2000 {
                let draw = |state: &mut u64| {
                    let exps: Exps = (0..nvars)
                        .map(|_| {
                            *state ^= *state << 13;
                            *state ^= *state >> 7;
                            *state ^= *state << 17;
                            (*state % 4) as u16
                        })
                        .collect();
                    Monomial::from_exps(exps)
                };
                let a = draw(&mut state);
                let b = draw(&mut state);
                assert_eq!(a.cmp(&b), by_entry(&a, &b), "at {:?} {:?}", a.exps, b.exps);
            }
        }
    }

    #[test]
    fn the_divisor_key_never_rejects_a_divisor() {
        let mut state = 0x243f_6a88_85a3_08d3u64;
        for nvars in 1..=18usize {
            for _ in 0..2000 {
                let draw = |state: &mut u64, bound: u64| {
                    let exps: Exps = (0..nvars)
                        .map(|_| {
                            *state ^= *state << 13;
                            *state ^= *state >> 7;
                            *state ^= *state << 17;
                            (*state % bound) as u16
                        })
                        .collect();
                    Monomial::from_exps(exps)
                };
                let a = draw(&mut state, 12);
                let b = draw(&mut state, 12);
                if a.divides(&b) {
                    assert!(
                        key_divides(a.divisor_key(), b.divisor_key()),
                        "the key rejected a divisor: {:?} | {:?}",
                        a.exps,
                        b.exps
                    );
                }
            }
        }
    }

    #[test]
    fn monomial_division_and_lcm() {
        let a = mono(&[2, 1]);
        let b = mono(&[1, 1]);
        assert!(b.divides(&a));
        assert_eq!(a.quotient(&b), Some(mono(&[1, 0])));
        assert!(!a.divides(&b));
        assert_eq!(a.lcm(&mono(&[1, 3])), mono(&[2, 3]));
        assert!(mono(&[2, 0]).is_coprime(&mono(&[0, 3])));
        assert!(!mono(&[2, 1]).is_coprime(&mono(&[0, 3])));
    }

    #[test]
    fn construction_sorts_and_merges() {
        let ring = ring();
        let f = ring
            .polynomial([(3, [0, 2]), (2, [2, 0]), (6, [2, 0])])
            .expect("fits");
        assert_eq!(f.terms.len(), 2);
        assert_eq!(f.terms[0].mono, mono(&[0, 2]));
        assert_eq!(f.terms[1].coeff.value(), 1);
    }

    #[test]
    fn terms_run_descending() {
        let ring = ring();
        let f = ring.parse_polynomial("y^2 + x^2 + x").expect("parses");
        let monomials: Vec<Vec<u16>> = f.terms().map(|(_, exps)| exps.to_vec()).collect();
        assert_eq!(monomials, vec![vec![2, 0], vec![0, 2], vec![1, 0]]);
        assert_eq!(f.leading_term(), Some((1, [2u16, 0].as_slice())));
        assert_eq!(f.degree(), Some(2));
    }

    #[test]
    fn display_writes_the_parse_syntax() {
        let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"]).expect("32003 is prime");
        let f = ring.parse_polynomial("x^2*y - 3*z + 1").expect("parses");
        assert_eq!(f.to_string(), "x^2*y + 32000*z + 1");
        assert_eq!(ring.zero().to_string(), "0");
    }

    #[test]
    fn subtraction_cancels_the_leading_term() {
        let ring = ring();
        let f = ring.polynomial([(2, [2, 0]), (1, [0, 0])]).expect("fits");
        let g = ring.polynomial([(1, [1, 0]), (1, [0, 0])]).expect("fits");
        let out = f
            .sub_scaled(&g, Felt::new(2, 7), &mono(&[1, 0]), 7)
            .expect("the product fits");
        assert_eq!(out.leading_term(), Some((5, [1u16, 0].as_slice())));
    }

    #[test]
    fn normal_form_reduces_completely() {
        let ring = ring();
        let f = ring
            .polynomial([(1, [2, 0]), (2, [1, 0]), (3, [0, 0])])
            .expect("fits");
        let g = ring.polynomial([(1, [1, 0])]).expect("fits");
        let r = f.normal_form(&[g], 7).expect("the product fits");
        assert_eq!(r.terms.len(), 1);
        assert_eq!(r.leading_term(), Some((3, [0u16, 0].as_slice())));
    }

    #[test]
    fn an_s_polynomial_cancels_the_lcm_term() {
        let ring = ring();
        let f = ring.polynomial([(2, [2, 0]), (1, [0, 0])]).expect("fits");
        let g = ring.polynomial([(3, [1, 1])]).expect("fits");
        let s = f.s_polynomial(&g, 7).expect("the product fits");
        assert!(s.terms.iter().all(|t| t.mono != mono(&[2, 1])));
    }
}
