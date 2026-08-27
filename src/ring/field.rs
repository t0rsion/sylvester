//! Arithmetic in the prime field F_p.
//!
//! Every value is held as its representative in `[0, p)`. The modulus
//! travels as an argument, so one type serves every prime the crate
//! supports.

/// One element of the prime field.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub(crate) struct Felt(u64);

impl Felt {
    /// Reduce a signed integer into the field.
    #[inline]
    pub(crate) fn new(value: i64, p: u64) -> Self {
        debug_assert!(p > 1, "the modulus must be above 1");
        Felt((value as i128).rem_euclid(p as i128) as u64)
    }

    /// Wrap a value that is already reduced.
    #[inline]
    pub(crate) fn from_residue(value: u64) -> Self {
        Felt(value)
    }

    /// The representative in `[0, p)`.
    #[inline]
    pub(crate) fn value(self) -> u64 {
        self.0
    }

    #[inline]
    pub(crate) fn one() -> Self {
        Felt(1)
    }

    #[inline]
    pub(crate) fn is_zero(self) -> bool {
        self.0 == 0
    }

    /// Add in the field.
    ///
    /// A ring holds a modulus of at most 2^31 - 1, so the sum fits a `u64`.
    #[inline]
    pub(crate) fn add(self, other: Self, p: u64) -> Self {
        debug_assert!(p < 1 << 32, "a ring modulus is below 2^31");
        let sum = self.0 + other.0;
        Felt(if sum >= p { sum - p } else { sum })
    }

    /// Subtract in the field.
    ///
    /// A ring holds a modulus of at most 2^31 - 1, so the difference fits
    /// a `u64`.
    #[inline]
    pub(crate) fn sub(self, other: Self, p: u64) -> Self {
        debug_assert!(p < 1 << 32, "a ring modulus is below 2^31");
        if self.0 >= other.0 {
            Felt(self.0 - other.0)
        } else {
            Felt(self.0 + (p - other.0))
        }
    }

    #[inline]
    pub(crate) fn neg(self, p: u64) -> Self {
        if self.is_zero() {
            self
        } else {
            Felt(p - self.0)
        }
    }

    /// Multiply in the field.
    ///
    /// A ring holds a modulus of at most 2^31 - 1, so the product fits a
    /// `u64`.
    #[inline]
    pub(crate) fn mul(self, other: Self, p: u64) -> Self {
        debug_assert!(self.0 < p && other.0 < p, "both factors must be below p");
        debug_assert!(p < 1 << 32, "a ring modulus is below 2^31");
        Felt(self.0 * other.0 % p)
    }

    /// The multiplicative inverse, by the extended Euclidean algorithm.
    ///
    /// Zero has no inverse. One is its own inverse, which every monic
    /// polynomial hits.
    pub(crate) fn inv(self, p: u64) -> Self {
        assert!(!self.is_zero(), "zero has no inverse in a field");
        if self.0 == 1 {
            return self;
        }
        let mut a = self.0 as i128;
        let mut b = p as i128;
        let mut x0: i128 = 1;
        let mut x1: i128 = 0;

        while b != 0 {
            let q = a / b;
            let tmp = b;
            b = a - q * b;
            a = tmp;

            let tmp = x1;
            x1 = x0 - q * x1;
            x0 = tmp;
        }

        Felt(x0.rem_euclid(p as i128) as u64)
    }

    /// Divide in the field.
    ///
    /// A divisor of one is common, because every basis element is monic.
    #[inline]
    pub(crate) fn div(self, other: Self, p: u64) -> Self {
        if other.0 == 1 {
            return self;
        }
        self.mul(other.inv(p), p)
    }
}

impl std::fmt::Display for Felt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Multiply two values below `n`, modulo `n`, for `n` up to `u64::MAX`.
///
/// The primality test below checks a candidate modulus before a ring
/// bounds it, so it widens through a `u128`. The field multiply in
/// [`Felt::mul`] holds only for a modulus below 2^32 and cannot serve
/// here.
#[inline]
fn wide_mul_residue(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 * b as u128) % n as u128) as u64
}

/// Report whether `n` is prime.
///
/// The test is a deterministic Miller-Rabin over the bases that decide
/// every 64-bit integer.
pub(crate) fn is_prime(n: u64) -> bool {
    let Some(screened) = screen_small_primes(n) else {
        let (d, r) = factor_twos(n - 1);
        return PRIMALITY_BASES
            .into_iter()
            .all(|base| passes_miller_rabin(base, d, r, n));
    };
    screened
}

const PRIMALITY_BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

fn screen_small_primes(n: u64) -> Option<bool> {
    if n < 2 {
        return Some(false);
    }
    for base in PRIMALITY_BASES {
        if n == base {
            return Some(true);
        }
        if n.is_multiple_of(base) {
            return Some(false);
        }
    }
    None
}

fn factor_twos(mut d: u64) -> (u64, u32) {
    let mut r = 0u32;
    while d.is_multiple_of(2) {
        d /= 2;
        r += 1;
    }
    (d, r)
}

fn passes_miller_rabin(base: u64, d: u64, r: u32, n: u64) -> bool {
    let mut x = pow_mod(base, d, n);
    if x == 1 || x == n - 1 {
        return true;
    }
    for _ in 1..r {
        x = wide_mul_residue(x, x, n);
        if x == n - 1 {
            return true;
        }
    }
    false
}

fn pow_mod(base: u64, mut exponent: u64, n: u64) -> u64 {
    let mut acc = 1u64;
    let mut base = base % n;
    while exponent > 0 {
        if exponent & 1 == 1 {
            acc = wide_mul_residue(acc, base, n);
        }
        base = wide_mul_residue(base, base, n);
        exponent >>= 1;
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::{Felt, is_prime};

    #[test]
    fn new_reduces_into_the_field() {
        let p = 7;
        assert_eq!(Felt::new(0, p).value(), 0);
        assert_eq!(Felt::new(8, p).value(), 1);
        assert_eq!(Felt::new(-1, p).value(), 6);
        assert_eq!(Felt::new(-8, p).value(), 6);
    }

    #[test]
    fn arithmetic_stays_in_the_field() {
        let p = 7;
        let a = Felt::new(3, p);
        let b = Felt::new(6, p);
        assert_eq!(a.add(b, p).value(), 2);
        assert_eq!(a.sub(b, p).value(), 4);
        assert_eq!(b.sub(a, p).value(), 3);
        assert_eq!(a.mul(b, p).value(), 4);
        assert_eq!(a.neg(p).value(), 4);
        assert_eq!(Felt::new(0, p).neg(p).value(), 0);
        assert_eq!(a.mul(a.inv(p), p).value(), 1);
        assert_eq!(Felt::new(2, p).div(Felt::new(4, p), p).value(), 4);
    }

    #[test]
    fn primality_agrees_with_trial_division_below_ten_thousand() {
        for n in 0u64..10_000 {
            let trial = n >= 2 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0);
            assert_eq!(is_prime(n), trial, "disagreement at {n}");
        }
    }

    #[test]
    fn primality_decides_large_moduli() {
        assert!(is_prime(1_073_741_827));
        assert!(is_prime(4_611_686_018_427_387_847));
        assert!(!is_prime(4_611_686_018_427_387_849));
        assert!(!is_prime(1u64 << 62));
    }
}
