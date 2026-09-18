//! Arithmetic for the verifier: prime field, monomials, polynomials.
//!
//! Every operation that builds a polynomial takes a [`Budget`]. It charges
//! the buffers it is about to reserve, plus the `live` terms the caller
//! holds beside them, before it reserves anything. Every loop over the
//! terms of an argument polls the deadline through the same budget, so one
//! multiplication cannot run past it.

use std::cmp::Ordering;

use super::error::VerifyError;
use super::limits::Budget;

/// An exponent of one variable.
///
/// A decoded exponent is at most 65535. The type is wider so a product of
/// two decoded monomials cannot overflow.
pub type Exp = u32;

pub(crate) fn add_mod(a: u64, b: u64, modulus: u64) -> u64 {
    let sum = a + b;
    if sum >= modulus { sum - modulus } else { sum }
}

pub(crate) fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

pub(crate) fn pow_mod(base: u64, exponent: u64, modulus: u64) -> u64 {
    let mut result = 1 % modulus;
    let mut acc = base % modulus;
    let mut left = exponent;
    while left > 0 {
        if left & 1 == 1 {
            result = mul_mod(result, acc, modulus);
        }
        acc = mul_mod(acc, acc, modulus);
        left >>= 1;
    }
    result
}

/// Return the inverse of `a` modulo a prime `modulus`.
pub(crate) fn inv_mod(a: u64, modulus: u64) -> u64 {
    pow_mod(a, modulus - 2, modulus)
}

const WITNESSES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// Report whether `n` is prime.
///
/// The Miller-Rabin test with these witnesses is exact for every `n` below
/// 3.3e24, so it decides every value this crate accepts.
pub(crate) fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for small in WITNESSES {
        if n.is_multiple_of(small) {
            return n == small;
        }
    }
    let mut d = n - 1;
    let mut shift = 0u32;
    while d.is_multiple_of(2) {
        d /= 2;
        shift += 1;
    }
    'witness: for base in WITNESSES {
        let mut x = pow_mod(base, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..shift {
            x = mul_mod(x, x, n);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// A monomial as an exponent vector over a fixed variable order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mono {
    pub(crate) exps: Vec<Exp>,
    pub(crate) deg: u64,
}

impl Mono {
    pub(crate) fn new(exps: Vec<Exp>) -> Self {
        let deg = exps.iter().map(|&e| e as u64).sum();
        Mono { exps, deg }
    }

    /// The exponent vector, one entry per variable.
    pub fn exps(&self) -> &[Exp] {
        &self.exps
    }

    /// The total degree.
    pub fn degree(&self) -> u64 {
        self.deg
    }

    /// The number of variables.
    pub fn nvars(&self) -> usize {
        self.exps.len()
    }

    /// Report whether this monomial divides `other`.
    pub fn divides(&self, other: &Mono) -> bool {
        self.exps.len() == other.exps.len()
            && self.exps.iter().zip(&other.exps).all(|(a, b)| a <= b)
    }

    /// Report whether the two monomials share no variable.
    pub fn is_coprime(&self, other: &Mono) -> bool {
        self.exps
            .iter()
            .zip(&other.exps)
            .all(|(&a, &b)| a == 0 || b == 0)
    }

    /// The product of two monomials.
    ///
    /// Every factor comes from the certificate, so every exponent is at
    /// most 65535 and the sum stays inside [`Exp`].
    pub(crate) fn mul(&self, other: &Mono) -> Mono {
        let exps = self
            .exps
            .iter()
            .zip(&other.exps)
            .map(|(&a, &b)| a + b)
            .collect();
        Mono {
            exps,
            deg: self.deg + other.deg,
        }
    }

    /// The least common multiple of the two monomials.
    pub(crate) fn lcm(&self, other: &Mono) -> Mono {
        Mono::new(
            self.exps
                .iter()
                .zip(&other.exps)
                .map(|(&a, &b)| a.max(b))
                .collect(),
        )
    }

    /// Divide a least common multiple by one of the monomials it covers.
    ///
    /// Every exponent of the multiple is at least the matching exponent of
    /// the divisor, so the subtraction cannot go below zero.
    fn divide_multiple(&self, divisor: &Mono) -> Mono {
        Mono::new(
            self.exps
                .iter()
                .zip(&divisor.exps)
                .map(|(&a, &b)| a - b)
                .collect(),
        )
    }
}

/// Compare two monomials under grevlex.
///
/// The larger degree is greater. At equal degree the exponent vectors are
/// scanned from the last variable back, and the smaller exponent at the
/// last difference is greater.
pub(crate) fn cmp_grevlex(a: &Mono, b: &Mono) -> Ordering {
    match a.deg.cmp(&b.deg) {
        Ordering::Equal => {
            for (&x, &y) in a.exps.iter().zip(&b.exps).rev() {
                match x.cmp(&y) {
                    Ordering::Equal => continue,
                    other => return other.reverse(),
                }
            }
            Ordering::Equal
        }
        other => other,
    }
}

/// A term: a coefficient in [1, p-1] and a monomial.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Term {
    pub(crate) coeff: u64,
    pub(crate) mono: Mono,
}

impl Term {
    pub(crate) fn new(coeff: u64, mono: Mono) -> Self {
        Term { coeff, mono }
    }

    /// The coefficient.
    pub fn coeff(&self) -> u64 {
        self.coeff
    }

    /// The monomial.
    pub fn mono(&self) -> &Mono {
        &self.mono
    }
}

/// A polynomial as a term list, sorted strictly descending under grevlex.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Poly {
    pub(crate) terms: Vec<Term>,
}

impl Poly {
    /// Wrap a term list.
    ///
    /// The caller keeps the descending order. Terms read from a
    /// certificate are checked separately.
    pub(crate) fn new(terms: Vec<Term>) -> Self {
        Poly { terms }
    }

    pub(crate) fn zero() -> Self {
        Poly { terms: Vec::new() }
    }

    /// Report whether the polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// The terms, sorted strictly descending under grevlex.
    pub fn terms(&self) -> &[Term] {
        &self.terms
    }

    /// The leading monomial, or `None` for the zero polynomial.
    pub fn lm(&self) -> Option<&Mono> {
        self.terms.first().map(|t| &t.mono)
    }

    /// The leading coefficient, or `None` for the zero polynomial.
    pub fn lc(&self) -> Option<u64> {
        self.terms.first().map(|t| t.coeff)
    }
}

/// Add two polynomials.
///
/// The arguments move into the sum through the merge the product uses, so
/// no exponent vector is copied. The merge charges them.
pub(crate) fn add(
    a: Poly,
    b: Poly,
    modulus: u64,
    budget: &mut Budget,
    live: usize,
) -> Result<Poly, VerifyError> {
    Ok(Poly {
        terms: merge(a.terms, b.terms, modulus, budget, live)?,
    })
}

/// Subtract `b` from `a`.
pub(crate) fn sub(
    a: Poly,
    b: Poly,
    modulus: u64,
    budget: &mut Budget,
    live: usize,
) -> Result<Poly, VerifyError> {
    let held = a.terms.len();
    let negated = scale(b, modulus - 1, modulus, budget, live + held)?;
    add(a, negated, modulus, budget, live)
}

/// Multiply a polynomial by a field element.
///
/// The result holds at most as many terms as the argument, and it reuses
/// the argument's term list rather than copying it.
pub(crate) fn scale(
    mut a: Poly,
    factor: u64,
    modulus: u64,
    budget: &mut Budget,
    live: usize,
) -> Result<Poly, VerifyError> {
    if factor == 0 {
        return Ok(Poly::zero());
    }
    budget.check_live(&[live, a.terms.len()])?;
    for term in &mut a.terms {
        budget.step()?;
        term.coeff = mul_mod(term.coeff, factor, modulus);
    }
    Ok(a)
}

/// Build a term list by mapping each term of `src`, dropping the `None`
/// results. Charges the source and the result buffer, next to the `live`
/// terms the caller holds, before it reserves, and polls once per term.
fn map_terms(
    src: &[Term],
    budget: &mut Budget,
    live: usize,
    mut f: impl FnMut(&Term) -> Option<Term>,
) -> Result<Vec<Term>, VerifyError> {
    budget.check_live(&[live, src.len()])?;
    let mut out = Vec::with_capacity(src.len());
    for term in src {
        budget.step()?;
        if let Some(mapped) = f(term) {
            out.push(mapped);
        }
    }
    Ok(out)
}

/// Multiply a polynomial by a single term.
///
/// A single term cannot cancel anything, so the order is preserved and the
/// result holds as many terms as the argument.
pub(crate) fn term_mul(
    a: &Poly,
    coeff: u64,
    mono: &Mono,
    modulus: u64,
    budget: &mut Budget,
    live: usize,
) -> Result<Poly, VerifyError> {
    if coeff == 0 {
        return Ok(Poly::zero());
    }
    let terms = map_terms(&a.terms, budget, live, |term| {
        Some(Term::new(
            mul_mod(term.coeff, coeff, modulus),
            term.mono.mul(mono),
        ))
    })?;
    Ok(Poly { terms })
}

/// Multiply two polynomials.
///
/// A term list times one term stays sorted, so the product needs no sort:
/// the multiplication halves the right factor, multiplies each half, and
/// merges the two results. The `live` terms the caller holds sit beside the
/// lists the recursion builds. The preflight charges one term per pair of
/// factors next to `live`, which rejects a product too large to hold before
/// the recursion runs; the recursion then charges every list it holds live,
/// so the multiplication never allocates past the cap.
pub(crate) fn mul(
    a: &Poly,
    b: &Poly,
    modulus: u64,
    budget: &mut Budget,
    live: usize,
) -> Result<Poly, VerifyError> {
    if a.is_zero() || b.is_zero() {
        return Ok(Poly::zero());
    }
    budget.check_deadline()?;
    budget.check_product(a.terms.len(), b.terms.len(), live)?;
    let mut terms = product(&a.terms, &b.terms, modulus, budget, live)?;
    if terms.len() * 2 <= terms.capacity() {
        terms.shrink_to_fit();
    }
    Ok(Poly { terms })
}

/// The product of two term lists, sorted strictly descending.
///
/// The recursion halves `right`, so it nests at most as deep as the base-2
/// logarithm of the term count. `live` is the term count of the buffers the
/// caller holds beside the list this call builds, so a merge one level up
/// charges the list built here as part of the buffers held live.
fn product(
    left: &[Term],
    right: &[Term],
    modulus: u64,
    budget: &mut Budget,
    live: usize,
) -> Result<Vec<Term>, VerifyError> {
    match right {
        [] => Ok(Vec::new()),
        [factor] => map_terms(left, budget, live, |term| {
            let coeff = mul_mod(term.coeff, factor.coeff, modulus);
            (coeff != 0).then(|| Term::new(coeff, term.mono.mul(&factor.mono)))
        }),
        _ => {
            let half = right.len() / 2;
            let low = product(left, &right[..half], modulus, budget, live)?;
            let high = product(left, &right[half..], modulus, budget, live + low.len())?;
            merge(low, high, modulus, budget, live)
        }
    }
}

/// Merge two descending term lists into one.
///
/// Equal monomials add, and a sum that vanishes drops out. The merge moves
/// the terms of both lists rather than copying their exponent vectors. The
/// budget charges the two lists it holds and the buffer it fills, next to
/// the `live` terms the caller holds.
fn merge(
    mut low: Vec<Term>,
    mut high: Vec<Term>,
    modulus: u64,
    budget: &mut Budget,
    live: usize,
) -> Result<Vec<Term>, VerifyError> {
    budget.check_live(&[live, low.len(), high.len(), low.len() + high.len()])?;
    let mut out = Vec::with_capacity(low.len() + high.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < low.len() && j < high.len() {
        budget.step()?;
        match cmp_grevlex(&low[i].mono, &high[j].mono) {
            Ordering::Greater => {
                out.push(take(&mut low[i]));
                i += 1;
            }
            Ordering::Less => {
                out.push(take(&mut high[j]));
                j += 1;
            }
            Ordering::Equal => {
                let coeff = add_mod(low[i].coeff, high[j].coeff, modulus);
                let mut term = take(&mut low[i]);
                if coeff != 0 {
                    term.coeff = coeff;
                    out.push(term);
                }
                i += 1;
                j += 1;
            }
        }
    }
    for term in low.drain(i..).chain(high.drain(j..)) {
        budget.step()?;
        out.push(term);
    }
    Ok(out)
}

/// Move a term out of a list, leaving an empty monomial in its place.
fn take(term: &mut Term) -> Term {
    Term {
        coeff: term.coeff,
        mono: Mono {
            exps: std::mem::take(&mut term.mono.exps),
            deg: term.mono.deg,
        },
    }
}

/// The S-polynomial of two nonzero polynomials.
///
/// Returns the zero polynomial if either argument is zero.
pub(crate) fn spoly(
    f: &Poly,
    g: &Poly,
    modulus: u64,
    budget: &mut Budget,
    live: usize,
) -> Result<Poly, VerifyError> {
    let (Some(lmf), Some(lmg)) = (f.lm(), g.lm()) else {
        return Ok(Poly::zero());
    };
    let (Some(lcf), Some(lcg)) = (f.lc(), g.lc()) else {
        return Ok(Poly::zero());
    };
    let lcm = lmf.lcm(lmg);
    let left = term_mul(
        f,
        inv_mod(lcf, modulus),
        &lcm.divide_multiple(lmf),
        modulus,
        budget,
        live,
    )?;
    let right = term_mul(
        g,
        inv_mod(lcg, modulus),
        &lcm.divide_multiple(lmg),
        modulus,
        budget,
        live + left.terms.len(),
    )?;
    sub(left, right, modulus, budget, live)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verify::error::Cap;
    use crate::verify::limits::{Limits, term_bytes};

    fn budget() -> Budget {
        Limits::default().budget(2)
    }

    /// A budget over two variables that admits `terms` terms, and no more.
    fn budget_for(terms: usize) -> Budget {
        Limits {
            max_intermediate_bytes: terms * term_bytes(2),
            ..Limits::default()
        }
        .budget(2)
    }

    fn mono(exps: &[u32]) -> Mono {
        Mono::new(exps.to_vec())
    }

    fn poly(terms: &[(u64, &[u32])]) -> Poly {
        Poly::new(
            terms
                .iter()
                .map(|(c, e)| Term::new(*c, mono(e)))
                .collect::<Vec<_>>(),
        )
    }

    #[test]
    fn a_merge_tail_observes_its_deadline() {
        let mut budget = Limits {
            deadline: Some(std::time::Instant::now()),
            ..Limits::default()
        }
        .budget(2);
        let tail = (0..crate::verify::limits::WORK_STRIDE)
            .rev()
            .map(|exponent| Term::new(1, mono(&[exponent as u32, 0])))
            .collect();
        assert_eq!(
            merge(Vec::new(), tail, 7, &mut budget, 0),
            Err(VerifyError::DeadlineExceeded)
        );
    }

    #[test]
    fn grevlex_prefers_degree() {
        assert_eq!(
            cmp_grevlex(&mono(&[3, 0]), &mono(&[0, 2])),
            Ordering::Greater
        );
        assert_eq!(
            cmp_grevlex(&mono(&[0, 3]), &mono(&[2, 0])),
            Ordering::Greater
        );
    }

    #[test]
    fn grevlex_breaks_ties_by_reversed_scan() {
        assert_eq!(
            cmp_grevlex(&mono(&[2, 0]), &mono(&[1, 1])),
            Ordering::Greater
        );
        assert_eq!(
            cmp_grevlex(&mono(&[1, 1]), &mono(&[0, 2])),
            Ordering::Greater
        );
        assert_eq!(
            cmp_grevlex(&mono(&[1, 0]), &mono(&[0, 1])),
            Ordering::Greater
        );
        assert_eq!(cmp_grevlex(&mono(&[1, 1]), &mono(&[1, 1])), Ordering::Equal);
    }

    #[test]
    fn grevlex_puts_the_last_variable_last() {
        assert_eq!(
            cmp_grevlex(&mono(&[1, 1, 1]), &mono(&[3, 0, 0])),
            Ordering::Less
        );
        assert_eq!(
            cmp_grevlex(&mono(&[1, 1, 1]), &mono(&[0, 0, 3])),
            Ordering::Greater
        );
    }

    #[test]
    fn primality_matches_trial_division() {
        for n in 0u64..2000 {
            let trial = n >= 2 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0);
            assert_eq!(is_prime(n), trial, "n = {n}");
        }
        assert!(is_prime(2147483647));
        assert!(!is_prime(2147483649));
        assert!(!is_prime(32003 * 3));
    }

    #[test]
    fn inverse_returns_the_field_inverse() {
        let p = 32003;
        for a in 1..500u64 {
            assert_eq!(mul_mod(a, inv_mod(a, p), p), 1);
        }
    }

    #[test]
    fn product_merges_and_drops_cancelling_terms() {
        let p = 7;
        let a = poly(&[(1, &[1, 0]), (1, &[0, 1])]);
        let b = poly(&[(1, &[1, 0]), (6, &[0, 1])]);
        let product = mul(&a, &b, p, &mut budget(), 0).expect("the budget holds");
        assert_eq!(product, poly(&[(1, &[2, 0]), (6, &[0, 2])]));
    }

    #[test]
    fn product_of_a_polynomial_with_its_negative_is_zero() {
        let p = 7;
        let a = poly(&[(1, &[1, 1]), (6, &[0, 0])]);
        assert!(
            sub(a.clone(), a, p, &mut budget(), 0)
                .expect("the budget holds")
                .is_zero()
        );
    }

    #[test]
    fn a_product_above_the_budget_is_rejected_before_it_allocates() {
        let p = 7;
        let a = poly(&[(1, &[1, 0]), (1, &[0, 1])]);
        let mut tight = budget_for(3);
        assert_eq!(
            mul(&a, &a, p, &mut tight, 0),
            Err(VerifyError::CapExceeded {
                cap: Cap::IntermediateBytes,
                limit: 3 * term_bytes(2)
            })
        );
    }

    #[test]
    fn a_product_the_caller_cannot_hold_beside_its_own_buffers_is_rejected() {
        // The two factors hold two terms each. The recursion holds the two
        // half products, two terms each, and the merged buffer, four terms,
        // live at the top merge: eight terms. That fits a budget for eight
        // terms and runs to the three-term product. One term the caller
        // holds beside them does not fit, and the check covers both.
        let p = 7;
        let a = poly(&[(1, &[1, 0]), (1, &[0, 1])]);
        assert_eq!(
            mul(&a, &a, p, &mut budget_for(8), 0).map(|r| r.terms().len()),
            Ok(3)
        );
        assert_eq!(
            mul(&a, &a, p, &mut budget_for(8), 1),
            Err(VerifyError::CapExceeded {
                cap: Cap::IntermediateBytes,
                limit: 8 * term_bytes(2)
            })
        );
    }

    #[test]
    fn a_product_that_collapses_gives_back_the_buffer_it_reserved() {
        // Every pair of factors meets in one monomial, so the merge drops
        // all but one term. The result must not hold the capacity of the
        // buffer the merge ran over.
        let p = 32003;
        let a = poly(&[(1, &[1, 0]), (1, &[0, 1])]);
        let b = poly(&[(1, &[1, 0]), (32002, &[0, 1])]);
        let product = mul(&a, &b, p, &mut budget(), 0).expect("the budget holds");
        assert_eq!(product, poly(&[(1, &[2, 0]), (32002, &[0, 2])]));
        assert!(product.terms.capacity() < 4);
    }

    #[test]
    fn a_sum_above_the_budget_is_rejected_before_it_allocates() {
        let p = 7;
        let a = poly(&[(1, &[1, 0]), (1, &[0, 1])]);
        let mut tight = budget_for(3);
        assert_eq!(
            add(a.clone(), a, p, &mut tight, 0),
            Err(VerifyError::CapExceeded {
                cap: Cap::IntermediateBytes,
                limit: 3 * term_bytes(2)
            })
        );
    }

    #[test]
    fn a_sum_is_charged_for_its_arguments_as_well_as_for_itself() {
        // The sum holds two terms, and the two arguments still live hold
        // two each. A budget for four terms is not enough for the three
        // buffers together.
        let p = 7;
        let a = poly(&[(1, &[1, 0]), (1, &[0, 1])]);
        assert!(add(a.clone(), a.clone(), p, &mut budget_for(4), 0).is_err());
        assert!(add(a.clone(), a, p, &mut budget_for(8), 0).is_ok());
    }

    #[test]
    fn spoly_cancels_the_leading_terms() {
        let p = 7;
        let g0 = poly(&[(1, &[2, 0]), (6, &[0, 1])]);
        let g1 = poly(&[(1, &[1, 1]), (6, &[0, 0])]);
        assert_eq!(
            spoly(&g0, &g1, p, &mut budget(), 0),
            Ok(poly(&[(6, &[0, 2]), (1, &[1, 0])]))
        );
    }

    #[test]
    fn spoly_of_a_pair_with_equal_multiples_is_zero() {
        let p = 7;
        let g0 = poly(&[(1, &[2, 0])]);
        let g1 = poly(&[(1, &[1, 1])]);
        assert!(
            spoly(&g0, &g1, p, &mut budget(), 0)
                .expect("the budget holds")
                .is_zero()
        );
    }
}
