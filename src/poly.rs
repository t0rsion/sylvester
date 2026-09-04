//! Monomials, terms, and polynomials of a [`PolynomialRing`].

use smallvec::SmallVec;
use std::cmp::Ordering;
use std::fmt;
use std::mem::size_of;

use crate::ring::{Domain, DomainOps, PolynomialRing, PrimeField};

/// The exponent vector of a monomial.
///
/// A ring with more than [`INLINE_VARIABLES`] variables still works: the
/// vector spills to the heap.
pub(crate) type Exps = SmallVec<[u16; INLINE_VARIABLES]>;

/// The number of exponents a monomial holds inline.
///
/// Eleven exponents are the most that fit the space `SmallVec` reserves
/// anyway, so `Monomial` stays 40 bytes and `Term` 48. A wider inline
/// array grows both and measures slower on every benchmark input.
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

/// A monomial as an exponent vector with its total degree.
///
/// Every monomial of one ring holds one exponent per variable. The
/// arithmetic checks that width in debug builds, because a mismatch is a
/// defect and not a value.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Monomial {
    pub(crate) exps: Exps,
    pub(crate) deg: u32,
}

impl Monomial {
    pub(crate) fn one(nvars: usize) -> Self {
        Monomial {
            exps: std::iter::repeat_n(0u16, nvars).collect(),
            deg: 0,
        }
    }

    pub(crate) fn from_exps(exps: Exps) -> Self {
        let deg = exps
            .iter()
            .try_fold(0u32, |acc, &e| acc.checked_add(e as u32))
            // a ring holds at most 65535 variables and every
            // exponent fits a u16, so the sum fits a u32.
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

    pub(crate) fn mul(&self, other: &Self) -> Self {
        debug_assert_eq!(self.nvars(), other.nvars(), "monomials must share nvars");
        let exps = self
            .exps
            .iter()
            .zip(&other.exps)
            // The engines check every critical pair's lcm against
            // compute::DEGREE_LIMIT before they build a product from it, so
            // this addition never overflows a u16 on a path that check
            // guards. The expect documents that invariant; it is not a
            // check a caller can trigger.
            .map(|(&a, &b)| a.checked_add(b).expect("monomial exponent overflow"))
            .collect();
        Monomial {
            exps,
            deg: self.deg + other.deg,
        }
    }

    /// Return `self * other`, or report an exponent past the width of a
    /// `u16`.
    ///
    /// Division by a basis bounds no degree in advance: reducing
    /// `x^65535*y` by `x^65535 + y^65535` needs `y^65536`. Every step that
    /// runs without such a bound multiplies here.
    pub(crate) fn checked_mul(&self, other: &Self) -> Result<Self, ExponentOverflow> {
        debug_assert_eq!(self.nvars(), other.nvars(), "monomials must share nvars");
        let mut exps = Exps::with_capacity(self.nvars());
        for (&a, &b) in self.exps.iter().zip(&other.exps) {
            exps.push(a.checked_add(b).ok_or(ExponentOverflow)?);
        }
        Ok(Monomial {
            exps,
            deg: self.deg + other.deg,
        })
    }

    pub(crate) fn quotient(&self, divisor: &Self) -> Option<Self> {
        debug_assert_eq!(self.nvars(), divisor.nvars(), "monomials must share nvars");
        if !divisor.divides(self) {
            return None;
        }
        let exps = self
            .exps
            .iter()
            .zip(&divisor.exps)
            .map(|(&a, &b)| a - b)
            .collect();
        Some(Monomial {
            exps,
            deg: self.deg - divisor.deg,
        })
    }

    pub(crate) fn lcm(&self, other: &Self) -> Self {
        debug_assert_eq!(self.nvars(), other.nvars(), "monomials must share nvars");
        let exps: Exps = self
            .exps
            .iter()
            .zip(&other.exps)
            .map(|(&a, &b)| a.max(b))
            .collect();
        Monomial::from_exps(exps)
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

    /// One bit per variable the monomial uses, for the first 64 variables.
    ///
    /// A divisor uses a subset of the variables of its multiple, so a bit
    /// set here and clear in the multiple rejects the divisor. Variables
    /// past the 64th set no bit, which only makes the test accept more.
    pub(crate) fn var_mask(&self) -> u64 {
        let mut mask = 0u64;
        for (index, &exp) in self.exps.iter().take(64).enumerate() {
            if exp != 0 {
                mask |= 1 << index;
            }
        }
        mask
    }

    /// Report whether no variable appears in both monomials.
    pub(crate) fn is_coprime(&self, other: &Self) -> bool {
        self.exps
            .iter()
            .zip(&other.exps)
            .all(|(&a, &b)| a == 0 || b == 0)
    }
}

/// The report that a monomial product passes the width of one exponent.
///
/// The width is 65,535, the largest value a `u16` holds. Each caller maps
/// this to its own limit error.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ExponentOverflow;

/// Grevlex, with the larger monomial the greater value.
impl Ord for Monomial {
    fn cmp(&self, other: &Self) -> Ordering {
        debug_assert_eq!(self.nvars(), other.nvars(), "monomials must share nvars");
        match self.deg.cmp(&other.deg) {
            Ordering::Equal => {
                for (&a, &b) in self.exps.iter().zip(&other.exps).rev() {
                    match a.cmp(&b) {
                        Ordering::Equal => continue,
                        ord => return ord.reverse(),
                    }
                }
                Ordering::Equal
            }
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
///
/// The domain parameter defaults to [`PrimeField`], as it does on
/// [`Polynomial`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Term<D: Domain = PrimeField> {
    pub(crate) coeff: D::Coeff,
    pub(crate) mono: Monomial,
}

/// A polynomial of a [`PolynomialRing`].
///
/// Build one with [`PolynomialRing::polynomial`] or
/// [`PolynomialRing::parse_polynomial`]. A polynomial carries its ring, so
/// even the zero polynomial names its variables and its domain. The
/// coefficients are reduced, the monomials are distinct, and no coefficient
/// is zero.
///
/// The domain parameter defaults to [`PrimeField`], so `Polynomial` alone
/// names a polynomial over `F_p`.
///
/// `Display` writes the syntax [`PolynomialRing::parse_polynomial`] reads,
/// with the terms in descending grevlex order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial<D: Domain = PrimeField> {
    ring: PolynomialRing<D>,
    /// Terms in ascending grevlex order, so the last term is the leading
    /// term.
    pub(crate) terms: Vec<Term<D>>,
}

impl<D: Domain> Polynomial<D> {
    /// The ring this polynomial belongs to.
    pub fn ring(&self) -> &PolynomialRing<D> {
        &self.ring
    }

    /// Report whether the polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// The terms, largest monomial first.
    ///
    /// Each item is the coefficient of the domain and one exponent per
    /// variable. Over [`PrimeField`] the coefficient is a [`Felt`], whose
    /// [`value`] is the representative in `[1, p)`.
    ///
    /// [`Felt`]: crate::Felt
    /// [`value`]: crate::Felt::value
    pub fn terms(&self) -> impl ExactSizeIterator<Item = (&D::Coeff, &[u16])> {
        self.terms
            .iter()
            .rev()
            .map(|term| (&term.coeff, term.mono.exps.as_slice()))
    }

    /// The largest term, or `None` for the zero polynomial.
    pub fn leading_term(&self) -> Option<(&D::Coeff, &[u16])> {
        self.lt()
            .map(|term| (&term.coeff, term.mono.exps.as_slice()))
    }

    /// The total degree, or `None` for the zero polynomial.
    ///
    /// Grevlex is graded, so the leading monomial carries the degree.
    pub fn degree(&self) -> Option<u32> {
        self.lt().map(|term| term.mono.deg)
    }

    /// The heap bytes this polynomial holds beyond its own value.
    ///
    /// One term costs its own size, the exponents that spill to the heap
    /// past [`INLINE_VARIABLES`] variables, and the heap bytes of its
    /// coefficient. A [`Felt`] holds none, so the count over a prime field
    /// is the term bytes alone. The classic backend and the interreduction
    /// size a polynomial with this. The certificate writer estimates the
    /// terms from its own per-term figure and adds
    /// [`Polynomial::coefficient_bytes`].
    ///
    /// [`Felt`]: crate::Felt
    pub(crate) fn heap_bytes(&self) -> usize {
        let per_term = size_of::<Term<D>>() + heap_exps_bytes(self.ring.nvars());
        let terms = self.terms.len().saturating_mul(per_term);
        terms.saturating_add(self.coefficient_bytes())
    }

    /// The heap bytes the coefficients of this polynomial hold.
    ///
    /// It is 0 over a prime field, where a coefficient is one `u64`.
    pub(crate) fn coefficient_bytes(&self) -> usize {
        let ops = self.ring.ops();
        self.terms.iter().fold(0usize, |bytes, term| {
            bytes.saturating_add(ops.heap_bytes(&term.coeff))
        })
    }

    /// Build a polynomial from terms that already hold the invariants.
    ///
    /// The terms run strictly ascending under grevlex, no coefficient is
    /// zero, and every monomial holds one exponent per variable.
    pub(crate) fn from_sorted_terms(ring: PolynomialRing<D>, terms: Vec<Term<D>>) -> Self {
        debug_assert!(
            terms.iter().all(|term| term.mono.nvars() == ring.nvars())
                && terms.windows(2).all(|pair| pair[0].mono < pair[1].mono),
            "a polynomial must hold ascending distinct monomials of the ring's width"
        );
        Polynomial { ring, terms }
    }

    /// Sort and merge terms whose coefficients are already nonzero.
    ///
    /// A merge that cancels drops the monomial, so the result holds no
    /// zero coefficient either.
    pub(crate) fn from_terms(ring: PolynomialRing<D>, mut terms: Vec<Term<D>>) -> Self {
        let ops = ring.ops();
        terms.sort_by(|a, b| a.mono.cmp(&b.mono));

        let mut out: Vec<Term<D>> = Vec::with_capacity(terms.len());
        for term in terms {
            if let Some(last) = out.last_mut()
                && last.mono == term.mono
            {
                match ops.add(&last.coeff, &term.coeff) {
                    Some(coeff) => last.coeff = coeff,
                    None => {
                        out.pop();
                    }
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

    /// Stop a step that mixes two rings, or that names the wrong domain
    /// arithmetic.
    ///
    /// The arithmetic travels as an argument, so this is what ties the
    /// argument back to the ring the result claims.
    fn check_operand(&self, other: &Self, ops: &D::Ops) {
        debug_assert!(ops == self.ring.ops(), "the arithmetic must be the ring's");
        debug_assert!(self.ring == other.ring, "both operands must share a ring");
    }

    fn with_terms(&self, terms: Vec<Term<D>>) -> Self {
        Polynomial {
            ring: self.ring.clone(),
            terms,
        }
    }

    pub(crate) fn lt(&self) -> Option<&Term<D>> {
        self.terms.last()
    }

    pub(crate) fn lc(&self) -> Option<&D::Coeff> {
        self.lt().map(|t| &t.coeff)
    }

    pub(crate) fn lm(&self) -> Option<&Monomial> {
        self.lt().map(|t| &t.mono)
    }

    pub(crate) fn pop_lt(&mut self) -> Option<Term<D>> {
        self.terms.pop()
    }

    /// Add one term whose coefficient is nonzero.
    pub(crate) fn push_term(&mut self, term: Term<D>, ops: &D::Ops) {
        if self.terms.is_empty() {
            self.terms.push(term);
            return;
        }

        match self
            .terms
            .binary_search_by(|probe| probe.mono.cmp(&term.mono))
        {
            Ok(pos) => match ops.add(&self.terms[pos].coeff, &term.coeff) {
                Some(coeff) => self.terms[pos].coeff = coeff,
                None => {
                    self.terms.remove(pos);
                }
            },
            Err(pos) => self.terms.insert(pos, term),
        }
    }

    pub(crate) fn add(&self, other: &Self, ops: &D::Ops) -> Self {
        self.check_operand(other, ops);
        let mut out: Vec<Term<D>> = Vec::with_capacity(self.terms.len() + other.terms.len());

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
                    out.push(b.clone());
                    j += 1;
                }
                Ordering::Equal => {
                    if let Some(coeff) = ops.add(&a.coeff, &b.coeff) {
                        out.push(Term {
                            coeff,
                            mono: a.mono.clone(),
                        });
                    }
                    i += 1;
                    j += 1;
                }
            }
        }

        out.extend(self.terms[i..].iter().cloned());
        out.extend(other.terms[j..].iter().cloned());
        self.with_terms(out)
    }

    pub(crate) fn sub(&self, other: &Self, ops: &D::Ops) -> Self {
        self.check_operand(other, ops);
        let mut out: Vec<Term<D>> = Vec::with_capacity(self.terms.len() + other.terms.len());

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
                        coeff: ops.neg(&b.coeff),
                        mono: b.mono.clone(),
                    });
                    j += 1;
                }
                Ordering::Equal => {
                    if let Some(coeff) = ops.sub(&a.coeff, &b.coeff) {
                        out.push(Term {
                            coeff,
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
            coeff: ops.neg(&term.coeff),
            mono: term.mono.clone(),
        }));

        self.with_terms(out)
    }

    /// Return `c * m * self`, for a nonzero `c`.
    ///
    /// The multiplication is unchecked. The engines bound the degree of a
    /// critical pair before they build any product from it. A caller
    /// without that bound takes [`Polynomial::scale_monomial_checked`].
    pub(crate) fn scale_monomial(&self, c: &D::Coeff, m: &Monomial, ops: &D::Ops) -> Self {
        if self.is_zero() {
            return self.zero_like();
        }
        if ops.is_one(c) {
            return self.shift_monomial(m);
        }
        let terms: Vec<Term<D>> = self
            .terms
            .iter()
            .map(|t| Term {
                coeff: ops.mul(&t.coeff, c),
                mono: t.mono.mul(m),
            })
            .collect();
        self.with_terms(terms)
    }

    /// Return `m * self`.
    ///
    /// The coefficients are copied. Grevlex is a monomial order, so the
    /// terms stay in ascending order.
    pub(crate) fn shift_monomial(&self, m: &Monomial) -> Self {
        let terms: Vec<Term<D>> = self
            .terms
            .iter()
            .map(|t| Term {
                coeff: t.coeff.clone(),
                mono: t.mono.mul(m),
            })
            .collect();
        self.with_terms(terms)
    }

    /// Return `c * m * self`, or report an exponent past the width.
    ///
    /// This is [`Polynomial::scale_monomial`] with every product checked.
    pub(crate) fn scale_monomial_checked(
        &self,
        c: &D::Coeff,
        m: &Monomial,
        ops: &D::Ops,
    ) -> Result<Self, ExponentOverflow> {
        if self.is_zero() {
            return Ok(self.zero_like());
        }
        if ops.is_one(c) {
            return self.shift_monomial_checked(m);
        }
        let mut terms: Vec<Term<D>> = Vec::with_capacity(self.terms.len());
        for t in &self.terms {
            terms.push(Term {
                coeff: ops.mul(&t.coeff, c),
                mono: t.mono.checked_mul(m)?,
            });
        }
        Ok(self.with_terms(terms))
    }

    /// Return `m * self`, or report an exponent past the width.
    fn shift_monomial_checked(&self, m: &Monomial) -> Result<Self, ExponentOverflow> {
        let mut terms: Vec<Term<D>> = Vec::with_capacity(self.terms.len());
        for t in &self.terms {
            terms.push(Term {
                coeff: t.coeff.clone(),
                mono: t.mono.checked_mul(m)?,
            });
        }
        Ok(self.with_terms(terms))
    }

    /// Return `self - c * m * other`, for a nonzero `c`.
    ///
    /// The multiplication is unchecked, as in
    /// [`Polynomial::scale_monomial`]. A caller without a degree bound
    /// takes [`Polynomial::sub_scaled_checked`].
    pub(crate) fn sub_scaled(
        &self,
        other: &Self,
        c: &D::Coeff,
        m: &Monomial,
        ops: &D::Ops,
    ) -> Self {
        self.check_operand(other, ops);
        if other.is_zero() {
            return self.clone();
        }

        let mut out: Vec<Term<D>> = Vec::with_capacity(self.terms.len() + other.terms.len());

        let mut i = 0usize;
        let mut j = 0usize;
        // The multiple of the term at `j`, held across the iterations that
        // advance `i` only.
        let mut b_mono: Option<Monomial> = None;

        while i < self.terms.len() && j < other.terms.len() {
            let a = &self.terms[i];
            let b0 = &other.terms[j];
            if b_mono.is_none() {
                b_mono = Some(b0.mono.mul(m));
            }
            let ord = a.mono.cmp(b_mono.as_ref().expect("the multiple is held"));
            match ord {
                Ordering::Less => {
                    out.push(a.clone());
                    i += 1;
                }
                Ordering::Greater => {
                    let coeff = ops.neg(&ops.mul(&b0.coeff, c));
                    let mono = b_mono.take().expect("the multiple is held");
                    out.push(Term { coeff, mono });
                    j += 1;
                }
                Ordering::Equal => {
                    let scaled = ops.mul(&b0.coeff, c);
                    // The monomials are equal, so the held multiple is
                    // a.mono.
                    let mono = b_mono.take().expect("the multiple is held");
                    if let Some(coeff) = ops.sub(&a.coeff, &scaled) {
                        out.push(Term { coeff, mono });
                    }
                    i += 1;
                    j += 1;
                }
            }
        }

        out.extend(self.terms[i..].iter().cloned());
        while j < other.terms.len() {
            let b0 = &other.terms[j];
            let coeff = ops.neg(&ops.mul(&b0.coeff, c));
            let mono = b_mono.take().unwrap_or_else(|| b0.mono.mul(m));
            out.push(Term { coeff, mono });
            j += 1;
        }

        self.with_terms(out)
    }

    /// The bytes `self - c * m * other` holds at most, before the
    /// subtraction runs.
    ///
    /// The result holds at most one term per term of the two operands.
    /// Over `Q` the coefficient of a term of `other` grows by the size of
    /// `c`, which is why `c` is read here. The count is an estimate by the
    /// meter of [`Polynomial::heap_bytes`]: a subtraction that carries can
    /// add one limb per term.
    pub(crate) fn sub_scaled_bytes(&self, other: &Self, c: &D::Coeff, ops: &D::Ops) -> usize {
        let per_term = size_of::<Term<D>>() + heap_exps_bytes(self.ring.nvars());
        let terms = self
            .terms
            .len()
            .saturating_add(other.terms.len())
            .saturating_mul(per_term);
        let scale = ops.heap_bytes(c);
        other.terms.iter().fold(
            terms.saturating_add(self.coefficient_bytes()),
            |bytes, term| {
                bytes
                    .saturating_add(ops.heap_bytes(&term.coeff))
                    .saturating_add(scale)
            },
        )
    }

    /// Return `self - c * m * other`, or report an exponent past the
    /// width.
    ///
    /// This is [`Polynomial::sub_scaled`] with every product checked. The
    /// multiplier reaches every term of `other`, so a product past the
    /// width can appear in the tail while both leading monomials fit.
    pub(crate) fn sub_scaled_checked(
        &self,
        other: &Self,
        c: &D::Coeff,
        m: &Monomial,
        ops: &D::Ops,
    ) -> Result<Self, ExponentOverflow> {
        self.check_operand(other, ops);
        if other.is_zero() {
            return Ok(self.clone());
        }

        let mut out: Vec<Term<D>> = Vec::with_capacity(self.terms.len() + other.terms.len());

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
                    let coeff = ops.neg(&ops.mul(&b0.coeff, c));
                    let mono = b_mono.take().expect("the multiple is held");
                    out.push(Term { coeff, mono });
                    j += 1;
                }
                Ordering::Equal => {
                    let scaled = ops.mul(&b0.coeff, c);
                    // The monomials are equal, so the held multiple is
                    // a.mono.
                    let mono = b_mono.take().expect("the multiple is held");
                    if let Some(coeff) = ops.sub(&a.coeff, &scaled) {
                        out.push(Term { coeff, mono });
                    }
                    i += 1;
                    j += 1;
                }
            }
        }

        out.extend(self.terms[i..].iter().cloned());
        append_scaled_tail(&mut out, &other.terms, j, b_mono, c, m, ops)?;

        Ok(self.with_terms(out))
    }

    pub(crate) fn make_monic(&self, ops: &D::Ops) -> Self {
        let Some(lc) = self.lc() else {
            return self.clone();
        };
        if ops.is_one(lc) {
            return self.clone();
        }
        let inv = ops.inv(lc);
        let terms: Vec<Term<D>> = self
            .terms
            .iter()
            .map(|t| Term {
                coeff: ops.mul(&t.coeff, &inv),
                mono: t.mono.clone(),
            })
            .collect();
        self.with_terms(terms)
    }

    #[cfg(test)]
    pub(crate) fn normal_form(&self, reducers: &[Polynomial<D>], ops: &D::Ops) -> Self {
        let reducers: Vec<&Polynomial<D>> = reducers.iter().collect();
        self.normal_form_refs(&reducers, ops)
    }

    pub(crate) fn normal_form_refs(&self, reducers: &[&Polynomial<D>], ops: &D::Ops) -> Self {
        let mut poly = self.clone();
        // The remainder grows by the largest term left, so it is built
        // descending and turned around once.
        let mut remainder: Vec<Term<D>> = Vec::new();

        while let Some(lt_p) = poly.lt().cloned() {
            let mut reduced = false;

            for g in reducers {
                let Some(lt_g) = g.lt() else { continue };
                if lt_g.mono.divides(&lt_p.mono) {
                    let m = lt_p
                        .mono
                        .quotient(&lt_g.mono)
                        // divides() implies a quotient exists.
                        .expect("divides() implies quotient()");
                    let scale = divide_coefficients::<D>(ops, &lt_p.coeff, &lt_g.coeff);
                    poly = poly.sub_scaled(g, &scale, &m, ops);
                    reduced = true;
                    break;
                }
            }

            if !reduced {
                let term = poly
                    .pop_lt()
                    // lt_p was Some, so poly is non-empty here.
                    .expect("polynomial should not be empty");
                remainder.push(term);
            }
        }

        remainder.reverse();
        Polynomial::from_sorted_terms(self.ring.clone(), remainder)
    }

    /// Return the S-polynomial of `self` and `other`.
    ///
    /// The multiplication is unchecked, as in
    /// [`Polynomial::scale_monomial`]. A caller without a degree bound
    /// takes [`Polynomial::s_polynomial_checked`].
    pub(crate) fn s_polynomial(&self, other: &Self, ops: &D::Ops) -> Self {
        self.check_operand(other, ops);
        if self.is_zero() || other.is_zero() {
            return self.zero_like();
        }
        // both polynomials are non-zero here.
        let lt_f = self.lt().expect("non-zero polynomial must have lt");
        // both polynomials are non-zero here.
        let lt_g = other.lt().expect("non-zero polynomial must have lt");
        let lcm = lt_f.mono.lcm(&lt_g.mono);

        // lcm is a multiple of lt_f.mono.
        let m_f = lcm.quotient(&lt_f.mono).expect("lm(f) divides lcm");
        // lcm is a multiple of lt_g.mono.
        let m_g = lcm.quotient(&lt_g.mono).expect("lm(g) divides lcm");

        let f_scaled = self.scale_monomial(&lt_g.coeff, &m_f, ops);
        let g_scaled = other.scale_monomial(&lt_f.coeff, &m_g, ops);
        f_scaled.sub(&g_scaled, ops)
    }
    /// Return the S-polynomial of `self` and `other`, or report an
    /// exponent past the width.
    ///
    /// The least common multiple of the two leading monomials fits by
    /// construction. The tail of each side is multiplied by a quotient of
    /// it, and that product can pass the width.
    pub(crate) fn s_polynomial_checked(
        &self,
        other: &Self,
        ops: &D::Ops,
    ) -> Result<Self, ExponentOverflow> {
        self.check_operand(other, ops);
        if self.is_zero() || other.is_zero() {
            return Ok(self.zero_like());
        }
        // both polynomials are non-zero here.
        let lt_f = self.lt().expect("non-zero polynomial must have lt");
        // both polynomials are non-zero here.
        let lt_g = other.lt().expect("non-zero polynomial must have lt");
        let lcm = lt_f.mono.lcm(&lt_g.mono);

        // lcm is a multiple of lt_f.mono.
        let m_f = lcm.quotient(&lt_f.mono).expect("lm(f) divides lcm");
        // lcm is a multiple of lt_g.mono.
        let m_g = lcm.quotient(&lt_g.mono).expect("lm(g) divides lcm");

        let f_scaled = self.scale_monomial_checked(&lt_g.coeff, &m_f, ops)?;
        let g_scaled = other.scale_monomial_checked(&lt_f.coeff, &m_g, ops)?;
        Ok(f_scaled.sub(&g_scaled, ops))
    }
}

/// Return `a / b` in the domain.
///
/// Every basis element is monic, so a divisor of 1 is the common case and
/// it needs no inversion.
pub(crate) fn divide_coefficients<D: Domain>(ops: &D::Ops, a: &D::Coeff, b: &D::Coeff) -> D::Coeff {
    if ops.is_one(b) {
        a.clone()
    } else {
        ops.mul(a, &ops.inv(b))
    }
}

fn append_scaled_tail<D: Domain>(
    out: &mut Vec<Term<D>>,
    terms: &[Term<D>],
    mut start: usize,
    held: Option<Monomial>,
    c: &D::Coeff,
    m: &Monomial,
    ops: &D::Ops,
) -> Result<(), ExponentOverflow> {
    if let Some(mono) = held {
        push_negated_scaled(out, &terms[start], c, mono, ops);
        start += 1;
    }
    for term in &terms[start..] {
        push_negated_scaled(out, term, c, term.mono.checked_mul(m)?, ops);
    }
    Ok(())
}

fn push_negated_scaled<D: Domain>(
    out: &mut Vec<Term<D>>,
    term: &Term<D>,
    c: &D::Coeff,
    mono: Monomial,
    ops: &D::Ops,
) {
    out.push(Term {
        coeff: ops.neg(&ops.mul(&term.coeff, c)),
        mono,
    });
}

fn write_term<D: Domain>(
    f: &mut fmt::Formatter<'_>,
    variables: &[String],
    ops: &D::Ops,
    coeff: &D::Coeff,
    exps: &[u16],
) -> fmt::Result {
    let negative = ops.is_negative(coeff);
    let magnitude;
    let value = if negative {
        magnitude = ops.neg(coeff);
        &magnitude
    } else {
        coeff
    };
    let mut written = false;
    if !ops.is_one(value) || exps.iter().all(|&exp| exp == 0) {
        ops.write(value, f)?;
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

impl<D: Domain> fmt::Display for Polynomial<D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return f.write_str("0");
        }
        let ops = self.ring.ops();
        let variables = self.ring.variables();
        for (index, (coeff, exps)) in self.terms().enumerate() {
            let negative = ops.is_negative(coeff);
            if index > 0 {
                f.write_str(if negative { " - " } else { " + " })?;
            } else if negative {
                f.write_str("-")?;
            }
            write_term::<D>(f, variables, ops, coeff, exps)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::{Felt, PolynomialRing, Rationals};

    fn ring() -> PolynomialRing {
        PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime")
    }

    fn mono(exps: &[u16]) -> Monomial {
        Monomial::from_exps(exps.iter().copied().collect())
    }

    /// The leading term as a residue and its exponents.
    fn lead(poly: &Polynomial) -> Option<(u64, &[u16])> {
        poly.leading_term()
            .map(|(coeff, exps)| (coeff.value(), exps))
    }

    #[test]
    fn grevlex_orders_by_degree_then_reverse_lex() {
        assert!(mono(&[3, 0]) > mono(&[0, 2]));
        assert!(mono(&[2, 0]) > mono(&[1, 1]));
        assert!(mono(&[1, 1]) > mono(&[0, 2]));
        assert!(mono(&[2, 1]) > mono(&[1, 2]));
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

    /// The engines bound the degree of a pair before they multiply, so
    /// `Monomial::mul` panics on an overflow it is not meant to see.
    /// Division by a basis takes `checked_mul` and reports it.
    #[test]
    fn a_checked_multiply_reports_an_exponent_past_the_width() {
        let half = mono(&[40000, 1]);
        let product = half.checked_mul(&mono(&[25535, 2])).expect("65535 fits");
        assert_eq!(product.exps.as_slice(), [65535, 3]);
        assert_eq!(product.deg, 65538);
        assert_eq!(half.checked_mul(&mono(&[25536, 0])), Err(ExponentOverflow));
        assert_eq!(
            mono(&[0, 65535]).checked_mul(&mono(&[0, 1])),
            Err(ExponentOverflow)
        );
    }

    /// The multiplier reaches every term of the reducer, so the tail can
    /// pass the width while both leading monomials fit.
    #[test]
    fn a_checked_subtraction_reports_the_tail_that_passes_the_width() {
        let ring = ring();
        let f = ring.polynomial([(1, [65535, 1])]).expect("fits");
        let g = ring
            .polynomial([(1, [65535, 0]), (1, [0, 65535])])
            .expect("fits");
        let one = Felt::new(1, 7);
        assert_eq!(
            f.sub_scaled_checked(&g, &one, &mono(&[0, 1]), ring.ops()),
            Err(ExponentOverflow)
        );
        assert_eq!(
            f.s_polynomial_checked(&g, ring.ops()),
            Err(ExponentOverflow)
        );
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
        assert_eq!(lead(&f), Some((1, [2u16, 0].as_slice())));
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
        let out = f.sub_scaled(&g, &Felt::new(2, 7), &mono(&[1, 0]), ring.ops());
        assert_eq!(lead(&out), Some((5, [1u16, 0].as_slice())));
    }

    #[test]
    fn normal_form_reduces_completely() {
        let ring = ring();
        let f = ring
            .polynomial([(1, [2, 0]), (2, [1, 0]), (3, [0, 0])])
            .expect("fits");
        let g = ring.polynomial([(1, [1, 0])]).expect("fits");
        let r = f.normal_form(&[g], ring.ops());
        assert_eq!(r.terms.len(), 1);
        assert_eq!(lead(&r), Some((3, [0u16, 0].as_slice())));
    }

    #[test]
    fn an_s_polynomial_cancels_the_lcm_term() {
        let ring = ring();
        let f = ring.polynomial([(2, [2, 0]), (1, [0, 0])]).expect("fits");
        let g = ring.polynomial([(3, [1, 1])]).expect("fits");
        let s = f.s_polynomial(&g, ring.ops());
        assert!(s.terms.iter().all(|t| t.mono != mono(&[2, 1])));
    }

    /// The memory meters of the crate size a polynomial with
    /// [`Polynomial::heap_bytes`], so every one of them counts the heap
    /// bytes of a coefficient. Over a prime field a coefficient holds
    /// none, and over the rationals it holds the numerator and the
    /// denominator.
    #[test]
    fn the_memory_meter_counts_the_coefficient_bytes() {
        let rationals: PolynomialRing<Rationals> =
            PolynomialRing::rationals(["x"]).expect("the name holds");
        let small = rationals
            .parse_polynomial("x + 1")
            .expect("the syntax holds");
        let large = rationals
            .parse_polynomial("340282366920938463463374607431768211457/3*x + 1")
            .expect("the syntax holds");
        assert!(small.coefficient_bytes() > 0);
        assert!(
            large.coefficient_bytes() > small.coefficient_bytes(),
            "a wider numerator costs more"
        );

        let terms = large.terms.len() * size_of::<Term<Rationals>>();
        assert_eq!(large.heap_bytes(), terms + large.coefficient_bytes());

        let prime = ring().parse_polynomial("x + 1").expect("the syntax holds");
        assert_eq!(prime.coefficient_bytes(), 0, "a Felt is one u64");
        assert_eq!(prime.heap_bytes(), prime.terms.len() * size_of::<Term>());
    }
}
