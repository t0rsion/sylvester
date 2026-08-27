//! Arithmetic for the v2 verifier: prime field, sparse monomials, and
//! polynomials.
//!
//! Every operation charges its work before it runs and checks the cap on
//! one arithmetic result before it allocates.
//!
//! Section 7.1 of the contract fixes the operations: scale, monomial
//! multiply, add, and subtract.

use std::cmp::Ordering;

use crate::verify::error::VerifyError;
use crate::verify::limits::Meter;

/// The largest exponent the contract allows.
///
/// The bound holds for every monomial the verifier forms, not only for the
/// monomials it decodes.
pub(super) const MAX_EXP: u64 = 65535;

/// The number of terms between two deadline polls inside one operation.
const POLL_TERMS: usize = 1 << 20;

/// An exponent of one variable.
///
/// A decoded exponent is at most 65535. The type is wider, so the sum of
/// two exponents cannot wrap before the bound test.
pub(super) type Exp = u32;

/// A monomial as a sparse support.
///
/// The entries hold one (variable, exponent) pair per variable with a
/// nonzero exponent, sorted by variable index ascending.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct Mono {
    entries: Vec<(u32, Exp)>,
    degree: u64,
}

impl Mono {
    /// Build a monomial from its support.
    ///
    /// The caller keeps the two support rules: the variable indices
    /// increase, and every exponent is in [1, 65535].
    pub(super) fn from_support(entries: Vec<(u32, Exp)>) -> Self {
        let degree = entries.iter().map(|&(_, e)| u64::from(e)).sum();
        Mono { entries, degree }
    }

    /// The number of (variable, exponent) pairs.
    pub(super) fn support(&self) -> usize {
        self.entries.len()
    }

    /// Report whether every exponent is zero.
    pub(super) fn is_identity(&self) -> bool {
        self.entries.is_empty()
    }

    /// Compare two monomials under `grevlex-v1`.
    ///
    /// The larger degree is greater. At equal degree the variables are
    /// scanned from the last one back, and the smaller exponent at the
    /// first difference is greater.
    pub(super) fn cmp_grevlex(&self, other: &Mono) -> Ordering {
        match self.degree.cmp(&other.degree) {
            Ordering::Equal => {}
            other => return other,
        }
        let (mut i, mut j) = (self.entries.len(), other.entries.len());
        loop {
            match (i, j) {
                (0, 0) => return Ordering::Equal,
                (0, _) => return Ordering::Greater,
                (_, 0) => return Ordering::Less,
                _ => {}
            }
            let (left_var, left_exp) = self.entries[i - 1];
            let (right_var, right_exp) = other.entries[j - 1];
            match left_var.cmp(&right_var) {
                Ordering::Greater => return Ordering::Less,
                Ordering::Less => return Ordering::Greater,
                Ordering::Equal => match left_exp.cmp(&right_exp) {
                    Ordering::Equal => {
                        i -= 1;
                        j -= 1;
                    }
                    Ordering::Less => return Ordering::Greater,
                    Ordering::Greater => return Ordering::Less,
                },
            }
        }
    }

    /// Report whether this monomial divides `other`.
    pub(super) fn divides(&self, other: &Mono) -> bool {
        let mut j = 0;
        for &(var, exp) in &self.entries {
            while j < other.entries.len() && other.entries[j].0 < var {
                j += 1;
            }
            match other.entries.get(j) {
                Some(&(other_var, other_exp)) if other_var == var && other_exp >= exp => j += 1,
                _ => return false,
            }
        }
        true
    }

    /// Report whether the two monomials share a variable.
    pub(super) fn shares_variable(&self, other: &Mono) -> bool {
        let (mut i, mut j) = (0, 0);
        while i < self.entries.len() && j < other.entries.len() {
            match self.entries[i].0.cmp(&other.entries[j].0) {
                Ordering::Less => i += 1,
                Ordering::Greater => j += 1,
                Ordering::Equal => return true,
            }
        }
        false
    }

    /// The product of two monomials.
    ///
    /// The product is a rejection when an exponent passes 65535. Two
    /// decoded exponents sum to at most 131070, so the sum is computed in
    /// [`Exp`] and then tested.
    pub(super) fn mul(&self, other: &Mono) -> Result<Mono, VerifyError> {
        let mut entries = Vec::with_capacity(self.entries.len() + other.entries.len());
        let (mut i, mut j) = (0, 0);
        while i < self.entries.len() || j < other.entries.len() {
            let left = self.entries.get(i).copied();
            let right = other.entries.get(j).copied();
            let entry = match (left, right) {
                (Some(a), Some(b)) if a.0 == b.0 => {
                    i += 1;
                    j += 1;
                    (a.0, a.1 + b.1)
                }
                (Some(a), Some(b)) if a.0 < b.0 => {
                    i += 1;
                    a
                }
                (Some(a), None) => {
                    i += 1;
                    a
                }
                (_, Some(b)) => {
                    j += 1;
                    b
                }
                (None, None) => unreachable!("the loop runs while one side has entries"),
            };
            if u64::from(entry.1) > MAX_EXP {
                return Err(VerifyError::ExponentOverflow { max: MAX_EXP });
            }
            entries.push(entry);
        }
        Ok(Mono::from_support(entries))
    }

    /// The least common multiple of the two monomials.
    pub(super) fn lcm(&self, other: &Mono) -> Mono {
        let mut entries = Vec::with_capacity(self.entries.len() + other.entries.len());
        let (mut i, mut j) = (0, 0);
        while i < self.entries.len() || j < other.entries.len() {
            let left = self.entries.get(i).copied();
            let right = other.entries.get(j).copied();
            match (left, right) {
                (Some(a), Some(b)) if a.0 == b.0 => {
                    entries.push((a.0, a.1.max(b.1)));
                    i += 1;
                    j += 1;
                }
                (Some(a), Some(b)) if a.0 < b.0 => {
                    entries.push(a);
                    i += 1;
                }
                (Some(a), None) => {
                    entries.push(a);
                    i += 1;
                }
                (_, Some(b)) => {
                    entries.push(b);
                    j += 1;
                }
                (None, None) => unreachable!("the loop runs while one side has entries"),
            }
        }
        Mono::from_support(entries)
    }

    /// Divide this monomial by one that divides it.
    ///
    /// The caller checks the divisibility first, so the subtraction cannot
    /// go below zero.
    pub(super) fn divide(&self, divisor: &Mono) -> Mono {
        let mut entries = Vec::with_capacity(self.entries.len());
        let mut j = 0;
        for &(var, exp) in &self.entries {
            while j < divisor.entries.len() && divisor.entries[j].0 < var {
                j += 1;
            }
            let taken = match divisor.entries.get(j) {
                Some(&(divisor_var, divisor_exp)) if divisor_var == var => {
                    j += 1;
                    divisor_exp
                }
                _ => 0,
            };
            if exp > taken {
                entries.push((var, exp - taken));
            }
        }
        Mono::from_support(entries)
    }

    /// The exponent vector, one entry per variable.
    pub(super) fn dense(&self, nvars: usize) -> Vec<Exp> {
        let mut exps = vec![0; nvars];
        for &(var, exp) in &self.entries {
            exps[var as usize] = exp;
        }
        exps
    }
}

/// Charge one monomial operation on two operands.
pub(super) fn charge_monos(meter: &mut Meter, a: &Mono, b: &Mono) -> Result<(), VerifyError> {
    meter.charge((a.support() + b.support()).max(1) as u64)
}

/// A term: a coefficient in [1, p-1] and a monomial.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Term {
    pub(super) coeff: u64,
    pub(super) mono: Mono,
}

/// A polynomial as a term list, strictly descending under grevlex.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct Poly {
    terms: Vec<Term>,
}

impl Poly {
    /// Wrap a term list.
    ///
    /// The caller keeps the descending order.
    pub(super) fn new(terms: Vec<Term>) -> Self {
        Poly { terms }
    }

    /// Report whether the polynomial is zero.
    pub(super) fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// The number of terms.
    pub(super) fn term_count(&self) -> usize {
        self.terms.len()
    }

    /// The terms, strictly descending under grevlex.
    pub(super) fn terms(&self) -> &[Term] {
        &self.terms
    }

    /// The leading monomial, or `None` for the zero polynomial.
    pub(super) fn lm(&self) -> Option<&Mono> {
        self.terms.first().map(|term| &term.mono)
    }

    /// The leading coefficient, or `None` for the zero polynomial.
    pub(super) fn lc(&self) -> Option<u64> {
        self.terms.first().map(|term| term.coeff)
    }

    /// The terms as a coefficient and a dense exponent vector each.
    pub(super) fn into_dense(self, nvars: usize) -> Vec<(u64, Vec<Exp>)> {
        self.terms
            .into_iter()
            .map(|term| (term.coeff, term.mono.dense(nvars)))
            .collect()
    }
}

fn poll_terms(meter: &Meter, index: usize) -> Result<(), VerifyError> {
    if index > 0 && index.is_multiple_of(POLL_TERMS) {
        return meter.poll();
    }
    Ok(())
}

/// Multiply every coefficient by a factor in [1, p-1].
///
/// A nonzero factor kills no term, so the term count does not change.
pub(super) fn scale(
    a: &Poly,
    factor: u64,
    modulus: u64,
    meter: &mut Meter,
) -> Result<Poly, VerifyError> {
    meter.charge((a.term_count() as u64).max(1))?;
    meter.check_intermediate(a.term_count())?;
    let mut terms = Vec::with_capacity(a.term_count());
    for (index, term) in a.terms().iter().enumerate() {
        poll_terms(meter, index)?;
        terms.push(Term {
            coeff: mul_mod(term.coeff, factor, modulus),
            mono: term.mono.clone(),
        });
    }
    Ok(Poly::new(terms))
}

/// Multiply every monomial by `mono`.
///
/// Multiplication by a monomial is strictly increasing under grevlex, so
/// the order holds and the term count does not change. Each product
/// charges the support sizes of its two operands, which section 8.1
/// requires.
pub(super) fn mono_mul(a: &Poly, mono: &Mono, meter: &mut Meter) -> Result<Poly, VerifyError> {
    meter.charge((a.term_count() as u64).max(1))?;
    meter.check_intermediate(a.term_count())?;
    let mut terms = Vec::with_capacity(a.term_count());
    for (index, term) in a.terms().iter().enumerate() {
        poll_terms(meter, index)?;
        charge_monos(meter, &term.mono, mono)?;
        terms.push(Term {
            coeff: term.coeff,
            mono: term.mono.mul(mono)?,
        });
    }
    Ok(Poly::new(terms))
}

/// Add two polynomials.
///
/// Each comparison charges the support sizes of its two operands, which
/// section 8.1 requires. A term list of wide monomials therefore costs more
/// than its term count.
pub(super) fn add(
    a: &Poly,
    b: &Poly,
    modulus: u64,
    meter: &mut Meter,
) -> Result<Poly, VerifyError> {
    let room = a.term_count() + b.term_count();
    meter.charge((room as u64).max(1))?;
    meter.check_intermediate(room)?;
    let mut terms = Vec::with_capacity(room);
    let (mut i, mut j) = (0, 0);
    while i < a.term_count() && j < b.term_count() {
        poll_terms(meter, i + j)?;
        charge_monos(meter, &a.terms()[i].mono, &b.terms()[j].mono)?;
        match a.terms()[i].mono.cmp_grevlex(&b.terms()[j].mono) {
            Ordering::Greater => {
                terms.push(a.terms()[i].clone());
                i += 1;
            }
            Ordering::Less => {
                terms.push(b.terms()[j].clone());
                j += 1;
            }
            Ordering::Equal => {
                let coeff = add_mod(a.terms()[i].coeff, b.terms()[j].coeff, modulus);
                if coeff != 0 {
                    terms.push(Term {
                        coeff,
                        mono: a.terms()[i].mono.clone(),
                    });
                }
                i += 1;
                j += 1;
            }
        }
    }
    terms.extend_from_slice(&a.terms()[i..]);
    terms.extend_from_slice(&b.terms()[j..]);
    Ok(Poly::new(terms))
}

/// Subtract `b` from `a`, by adding `b` scaled by `p - 1`.
pub(super) fn sub(
    a: &Poly,
    b: &Poly,
    modulus: u64,
    meter: &mut Meter,
) -> Result<Poly, VerifyError> {
    let negated = scale(b, modulus - 1, modulus, meter)?;
    add(a, &negated, modulus, meter)
}

pub(super) fn add_mod(a: u64, b: u64, modulus: u64) -> u64 {
    let sum = a + b;
    if sum >= modulus { sum - modulus } else { sum }
}

pub(super) fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

fn pow_mod(base: u64, exponent: u64, modulus: u64) -> u64 {
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

/// The Miller-Rabin witnesses section 8.1 names.
const WITNESSES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// Report whether `n` is prime.
///
/// The test is deterministic over these witnesses for every `n` below
/// 3.3e24, so it decides every modulus the header permits. A modulus that
/// one witness divides is prime only when it equals that witness.
pub(super) fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for small in WITNESSES {
        if n.is_multiple_of(small) {
            return n == small;
        }
    }
    let mut odd = n - 1;
    let mut shift = 0u32;
    while odd.is_multiple_of(2) {
        odd /= 2;
        shift += 1;
    }
    'witness: for base in WITNESSES {
        let mut x = pow_mod(base, odd, n);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verify::limits::Limits;

    fn mono(entries: &[(u32, Exp)]) -> Mono {
        Mono::from_support(entries.to_vec())
    }

    fn meter() -> Meter {
        Meter::new(&Limits::default())
    }

    #[test]
    fn grevlex_prefers_the_larger_degree() {
        assert_eq!(
            mono(&[(0, 3)]).cmp_grevlex(&mono(&[(1, 2)])),
            Ordering::Greater
        );
    }

    #[test]
    fn grevlex_puts_the_first_variable_above_the_last() {
        assert_eq!(
            mono(&[(0, 1)]).cmp_grevlex(&mono(&[(1, 1)])),
            Ordering::Greater
        );
        assert_eq!(
            mono(&[(0, 2)]).cmp_grevlex(&mono(&[(0, 1), (1, 1)])),
            Ordering::Greater
        );
        assert_eq!(
            mono(&[(0, 1), (1, 1), (2, 1)]).cmp_grevlex(&mono(&[(2, 3)])),
            Ordering::Greater
        );
        assert_eq!(
            mono(&[(0, 1), (1, 1), (2, 1)]).cmp_grevlex(&mono(&[(0, 3)])),
            Ordering::Less
        );
    }

    #[test]
    fn the_identity_monomial_is_below_every_other() {
        assert_eq!(
            Mono::from_support(Vec::new()).cmp_grevlex(&mono(&[(0, 1)])),
            Ordering::Less
        );
        assert_eq!(
            Mono::from_support(Vec::new()).cmp_grevlex(&Mono::from_support(Vec::new())),
            Ordering::Equal
        );
    }

    #[test]
    fn divisibility_reads_the_sparse_support() {
        assert!(mono(&[(0, 1)]).divides(&mono(&[(0, 2), (1, 1)])));
        assert!(!mono(&[(0, 3)]).divides(&mono(&[(0, 2), (1, 1)])));
        assert!(!mono(&[(2, 1)]).divides(&mono(&[(0, 2), (1, 1)])));
        assert!(Mono::from_support(Vec::new()).divides(&mono(&[(0, 2)])));
    }

    #[test]
    fn a_product_past_the_exponent_bound_is_a_rejection() {
        let big = mono(&[(0, 65535)]);
        assert_eq!(
            big.mul(&mono(&[(0, 1)])),
            Err(VerifyError::ExponentOverflow { max: MAX_EXP })
        );
        assert_eq!(big.mul(&Mono::from_support(Vec::new())), Ok(big.clone()));
    }

    #[test]
    fn lcm_and_divide_invert_each_other() {
        let a = mono(&[(0, 2), (1, 1)]);
        let b = mono(&[(1, 3)]);
        let lcm = a.lcm(&b);
        assert_eq!(lcm, mono(&[(0, 2), (1, 3)]));
        assert_eq!(lcm.divide(&a), mono(&[(1, 2)]));
        assert_eq!(lcm.divide(&b), mono(&[(0, 2)]));
    }

    #[test]
    fn shared_variables_decide_the_product_criterion() {
        assert!(!mono(&[(0, 1)]).shares_variable(&mono(&[(1, 1)])));
        assert!(mono(&[(0, 1), (1, 1)]).shares_variable(&mono(&[(1, 1)])));
        assert!(!Mono::from_support(Vec::new()).shares_variable(&mono(&[(0, 1)])));
    }

    #[test]
    fn subtraction_cancels_equal_polynomials() {
        let f = Poly::new(vec![
            Term {
                coeff: 1,
                mono: mono(&[(0, 2)]),
            },
            Term {
                coeff: 3,
                mono: Mono::from_support(Vec::new()),
            },
        ]);
        let mut meter = meter();
        assert!(sub(&f, &f, 7, &mut meter).expect("the caps hold").is_zero());
    }

    #[test]
    fn primality_matches_trial_division() {
        for n in 0u64..2000 {
            let trial = n >= 2 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0);
            assert_eq!(is_prime(n), trial, "n = {n}");
        }
        assert!(is_prime(2147483647));
        assert!(!is_prime(1373653));
    }
}
