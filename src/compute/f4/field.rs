//! Field context for the F4 engine (design section 8).
//!
//! [`Small31`] implements [`FieldOps`] for a prime below 2^31. It
//! multiplies with Shoup's method and delays reduction until a sweep or a
//! read.
//!
//! A run builds one context and instantiates the kernel with it. A multiply
//! therefore carries no runtime method choice.

/// The exclusive upper bound on a modulus, 2^31.
///
/// `PolynomialRing::prime_field` accepts `2 <= p <= 2^31 - 1`, so every
/// prime the crate supports fits the kernel.
pub(crate) const MODULUS_LIMIT: u64 = 1 << 31;

/// The number of `axpy` applications one accumulator lane absorbs.
///
/// `capacity` is the largest value the lane holds, `start` the largest
/// value in a lane after a sweep, and `addend` the largest value one
/// application adds to one lane. `addend` must be at least 1. The count
/// is exact: `start + n * addend` is at most `capacity`, and one more
/// application passes `capacity`.
pub(crate) const fn applications_before_sweep(capacity: u64, start: u64, addend: u64) -> u64 {
    (capacity - start) / addend
}

/// Arithmetic in F_p for one run of the F4 engine.
///
/// A coefficient is a value below `p`. The reduction accumulator is one
/// `u64` lane per column. A lane holds a value below `p` after a scatter
/// or a sweep, and [`FieldOps::axpy`] adds to it without reducing. The
/// caller counts applications to the current accumulator and calls
/// [`FieldOps::sweep`] at [`FieldOps::applications_between_sweeps`]. That
/// count is what keeps a lane inside `u64`.
pub(crate) trait FieldOps {
    /// Every stored value is below `p`.
    type Coeff: Copy;

    /// Add two coefficients.
    ///
    /// `tests/f4_field.rs` includes this module by path. The library test
    /// target compiles the method without that caller.
    #[cfg(test)]
    #[allow(dead_code)]
    fn add(&self, a: Self::Coeff, b: Self::Coeff) -> Self::Coeff;

    /// Subtract `b` from `a`.
    fn sub(&self, a: Self::Coeff, b: Self::Coeff) -> Self::Coeff;

    /// Reduce one accumulator lane to its representative in `[0, p)`.
    fn reduce_acc(&self, acc: u64) -> Self::Coeff;

    /// The inverse of a nonzero coefficient.
    ///
    /// The caller passes a value in `1 ..= p - 1`. The engine calls this
    /// once per new pivot, never per column.
    fn inv(&self, a: Self::Coeff) -> Self::Coeff;

    /// Fill `out` with the multiply precomputation for `vals`.
    ///
    /// `out` is cleared first, so a workspace vector keeps its capacity
    /// and a caller can reserve it under the memory budget. The result
    /// belongs to one coefficient vector and stops being valid as soon as
    /// that vector changes. Rebuild it after normalization and after
    /// interreduction, before the row reduces anything.
    fn precompute(&self, vals: &[Self::Coeff], out: &mut Vec<u64>);

    /// Add `factor * vals[k]` into lane `cols[k]`, for every `k`.
    ///
    /// `cols`, `vals`, and `shoup` describe one row: `cols` and `vals`
    /// have the same length, and `shoup` is what
    /// [`FieldOps::precompute`] wrote for `vals`. Column indices are
    /// strictly increasing, so no lane takes two addends from one call,
    /// and each is below `acc.len()`. `factor` is below `p`. The lanes
    /// are not reduced.
    fn axpy(
        &self,
        acc: &mut [u64],
        cols: &[u32],
        vals: &[Self::Coeff],
        shoup: &[u64],
        factor: Self::Coeff,
    );

    /// Reduce every lane of `acc` to its representative in `[0, p)`.
    fn sweep(&self, acc: &mut [u64]);

    /// The number of [`FieldOps::axpy`] applications between two sweeps.
    ///
    /// The count holds for an accumulator whose lanes are below `p`, so
    /// it is counted from the scatter or from the last sweep. A debug
    /// build checks each lane in [`FieldOps::axpy`], so a missing sweep
    /// fails there instead of wrapping.
    fn applications_between_sweeps(&self) -> u64;
}

/// The high 64 bits of a 64-by-64 bit product.
///
/// The `u128` multiply is one widening machine multiply. Nothing else in
/// the kernel widens past `u64`.
#[inline]
const fn mulhi(a: u64, b: u64) -> u64 {
    (((a as u128) * (b as u128)) >> 64) as u64
}

/// Extended Euclid, the inverse of `a` modulo `p`.
fn inverse(a: u32, p: u32) -> u32 {
    debug_assert!(a != 0 && a < p, "the value must be in 1 ..= p - 1");
    let (mut r, mut new_r) = (p as i64, a as i64);
    let (mut t, mut new_t) = (0i64, 1i64);
    while new_r != 0 {
        let q = r / new_r;
        (r, new_r) = (new_r, r - q * new_r);
        (t, new_t) = (new_t, t - q * new_t);
    }
    debug_assert_eq!(r, 1, "the modulus must be prime");
    if t < 0 {
        (t + p as i64) as u32
    } else {
        t as u32
    }
}

/// The Shoup kernel for a prime below 2^31.
///
/// Shoup's method leaves a product below `2p` and the kernel does not
/// correct it. The accumulator is reduced when it is read, so the
/// correction buys nothing and costs a compare per nonzero. One
/// application therefore adds at most `2p - 1` to a lane, and a lane that
/// starts below `p` absorbs `(u64::MAX - (p - 1)) / (2p - 1)`
/// applications, which is above 2^32 for every prime the ring accepts.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Small31 {
    p: u32,
    p64: u64,
    /// Barrett's multiplier, `floor(2^64 / p)`.
    mu: u64,
    per_sweep: u64,
    /// The largest lane value that still takes one more application.
    lane_limit: u64,
}

impl Small31 {
    /// Build the context for the modulus `p`, with `2 <= p < 2^31`.
    pub(crate) fn new(p: u32) -> Self {
        assert!(p >= 2, "the modulus must be at least 2");
        assert!((p as u64) < MODULUS_LIMIT, "the modulus must be below 2^31");
        let p64 = p as u64;
        Self {
            p,
            p64,
            mu: ((1u128 << 64) / p64 as u128) as u64,
            per_sweep: applications_before_sweep(u64::MAX, p64 - 1, 2 * p64 - 1),
            lane_limit: u64::MAX - (2 * p64 - 1),
        }
    }

    /// Shoup's precomputation for one coefficient, `floor(w * 2^32 / p)`.
    ///
    /// `w < p < 2^31` keeps `w << 32` inside `u64`.
    #[inline]
    fn shoup_of(&self, w: u32) -> u64 {
        debug_assert!((w as u64) < self.p64, "a coefficient is below p");
        ((w as u64) << 32) / self.p64
    }

    fn debug_check_axpy(&self, cols: &[u32], vals: &[u32], shoup: &[u64], factor: u32) {
        debug_assert_eq!(cols.len(), vals.len(), "a row has one column per value");
        debug_assert_eq!(shoup.len(), vals.len(), "the Shoup array is stale");
        debug_assert!(factor < self.p, "the factor is below p");
        debug_assert!(
            cols.windows(2).all(|w| w[0] < w[1]),
            "column indices are strictly increasing"
        );
    }

    /// Return `value * factor` modulo `p`, in `[0, 2p)`.
    #[inline]
    fn shoup_product(&self, value: u32, precomputed: u64, factor: u64) -> u64 {
        let quotient = (precomputed * factor) >> 32;
        let product = u64::from(value) * factor - quotient * self.p64;
        debug_assert!(product < 2 * self.p64, "Shoup left a value at or above 2p");
        product
    }
}

impl FieldOps for Small31 {
    type Coeff = u32;

    #[cfg(test)]
    #[inline]
    fn add(&self, a: u32, b: u32) -> u32 {
        debug_assert!(a < self.p && b < self.p, "both terms are below p");
        let sum = a + b;
        if sum >= self.p { sum - self.p } else { sum }
    }

    #[inline]
    fn sub(&self, a: u32, b: u32) -> u32 {
        debug_assert!(a < self.p && b < self.p, "both terms are below p");
        // a + (p - b) is at most 2p - 1, which is inside u32.
        let d = a + (self.p - b);
        if d >= self.p { d - self.p } else { d }
    }

    #[inline]
    fn reduce_acc(&self, acc: u64) -> u32 {
        // mu = floor(2^64 / p) makes q either floor(acc / p) or one less,
        // so r is below 2p and one correction is enough. q * p is at most
        // acc, so it does not overflow.
        let q = mulhi(acc, self.mu);
        let r = acc - q * self.p64;
        let r = if r >= self.p64 { r - self.p64 } else { r };
        debug_assert!(r < self.p64, "Barrett left a value at or above p");
        r as u32
    }

    fn inv(&self, a: u32) -> u32 {
        inverse(a, self.p)
    }

    fn precompute(&self, vals: &[u32], out: &mut Vec<u64>) {
        out.clear();
        out.extend(vals.iter().map(|&w| self.shoup_of(w)));
    }

    #[inline]
    fn axpy(&self, acc: &mut [u64], cols: &[u32], vals: &[u32], shoup: &[u64], factor: u32) {
        self.debug_check_axpy(cols, vals, shoup, factor);
        let factor = u64::from(factor);
        for ((&col, &value), &precomputed) in cols.iter().zip(vals).zip(shoup) {
            let product = self.shoup_product(value, precomputed, factor);
            debug_assert!(
                acc[col as usize] <= self.lane_limit,
                "the lane needed a sweep"
            );
            acc[col as usize] += product;
        }
    }

    fn sweep(&self, acc: &mut [u64]) {
        for lane in acc.iter_mut() {
            *lane = self.reduce_acc(*lane) as u64;
        }
    }

    #[inline]
    fn applications_between_sweeps(&self) -> u64 {
        self.per_sweep
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BENCH_P: u32 = 1073741827;

    #[test]
    fn mulhi_matches_the_u128_product() {
        let cases = [
            (0u64, 0u64),
            (1, u64::MAX),
            (u64::MAX, u64::MAX),
            (1 << 63, 3),
            (0x0123_4567_89ab_cdef, 0xfedc_ba98_7654_3210),
        ];
        for (a, b) in cases {
            let want = ((a as u128 * b as u128) >> 64) as u64;
            assert_eq!(mulhi(a, b), want, "a = {a}, b = {b}");
        }
    }

    #[test]
    fn barrett_multiplier_is_the_floor_of_two_to_the_64() {
        for p in [2u32, 3, 7, 65521, 1 << 30, BENCH_P, (1 << 31) - 1] {
            let ctx = Small31::new(p);
            assert_eq!(ctx.mu, ((1u128 << 64) / p as u128) as u64, "p = {p}");
        }
    }

    #[test]
    fn shoup_precomputation_is_the_floor_of_the_scaled_quotient() {
        let ctx = Small31::new(BENCH_P);
        for w in [0u32, 1, 2, BENCH_P / 2, BENCH_P - 2, BENCH_P - 1] {
            let want = ((w as u128) << 32) / BENCH_P as u128;
            assert_eq!(ctx.shoup_of(w), want as u64, "w = {w}");
        }
    }

    #[test]
    fn reduce_acc_covers_the_boundary_values() {
        for p in [2u32, 3, 97, 65521, BENCH_P, (1 << 31) - 1] {
            let ctx = Small31::new(p);
            let p64 = p as u64;
            for x in [0, 1, p64 - 1, p64, p64 + 1, 2 * p64 - 1, 1 << 32, u64::MAX] {
                assert_eq!(ctx.reduce_acc(x) as u64, x % p64, "p = {p}, x = {x}");
            }
        }
    }

    #[test]
    fn the_sweep_bound_is_above_two_to_the_thirty_two() {
        assert!(Small31::new(BENCH_P).applications_between_sweeps() > 1 << 32);
    }

    #[test]
    fn the_shoup_product_stays_below_two_p() {
        let ctx = Small31::new(BENCH_P);
        let mut acc = [0u64; 1];
        for w in [1u32, 2, BENCH_P / 3, BENCH_P - 1] {
            for f in [1u32, 7, BENCH_P / 2, BENCH_P - 1] {
                let mut shoup = Vec::new();
                ctx.precompute(&[w], &mut shoup);
                acc[0] = 0;
                ctx.axpy(&mut acc, &[0], &[w], &shoup, f);
                assert!(acc[0] < 2 * BENCH_P as u64, "w = {w}, f = {f}");
                assert_eq!(
                    ctx.reduce_acc(acc[0]) as u64,
                    (w as u64 * f as u64) % BENCH_P as u64,
                    "w = {w}, f = {f}"
                );
            }
        }
    }

    #[test]
    fn the_sweep_bound_is_exact() {
        for (capacity, start, addend) in [(u16::MAX as u64, 6, 7), (255, 10, 1), (1000, 999, 1000)]
        {
            let n = applications_before_sweep(capacity, start, addend);
            assert!(start + n * addend <= capacity);
            assert!(start + (n + 1) * addend > capacity);
        }
    }
}
