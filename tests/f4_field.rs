//! Reference tests for the F4 field context (design section 8.3).
//!
//! Every operation is compared against a `u128` model that reduces with
//! `%`. The checks run through the `FieldOps` trait, which is the
//! interface the reduction kernel sees.

#[allow(dead_code)]
#[path = "../src/compute/f4/field.rs"]
mod field;

use field::{FieldOps, Small31, applications_before_sweep};
use proptest::prelude::*;

/// The modulus the release gate measures, 2^30 + 3.
const BENCH_P: u32 = 1073741827;

/// Primes used where a test wants one large modulus per size class.
const LARGE_PRIMES: [u32; 3] = [1073741789, BENCH_P, (1 << 31) - 1];

const SMALL_PRIMES: [u32; 6] = [2, 3, 5, 7, 13, 251];

fn add_model(a: u32, b: u32, p: u32) -> u32 {
    ((a as u128 + b as u128) % p as u128) as u32
}

fn sub_model(a: u32, b: u32, p: u32) -> u32 {
    ((a as i128 - b as i128).rem_euclid(p as i128)) as u32
}

fn mul_model(a: u32, b: u32, p: u32) -> u32 {
    ((a as u128 * b as u128) % p as u128) as u32
}

fn reduce_model(x: u64, p: u32) -> u32 {
    ((x as u128) % p as u128) as u32
}

/// The inverse by Fermat's theorem, which needs a prime `p`.
fn inv_model(a: u32, p: u32) -> u32 {
    let mut result: u128 = 1;
    let mut base = a as u128;
    let mut e = p - 2;
    let m = p as u128;
    while e > 0 {
        if e & 1 == 1 {
            result = result * base % m;
        }
        base = base * base % m;
        e >>= 1;
    }
    result as u32
}

/// One modular multiply through the kernel's `axpy` path.
fn mul_via_axpy<F: FieldOps<Coeff = u32>>(ctx: &F, a: u32, b: u32) -> u32 {
    let mut acc = [0u64; 1];
    let mut shoup = Vec::new();
    ctx.precompute(&[a], &mut shoup);
    ctx.axpy(&mut acc, &[0], &[a], &shoup, b);
    ctx.reduce_acc(acc[0])
}

/// A seeded generator, so a failing case is reproducible from its seed.
struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound
    }
}

fn check_all_pairs<F: FieldOps<Coeff = u32>>(ctx: &F, p: u32) {
    for a in 0..p {
        for b in 0..p {
            assert_eq!(ctx.add(a, b), add_model(a, b, p), "add, p = {p}");
            assert_eq!(ctx.sub(a, b), sub_model(a, b, p), "sub, p = {p}");
            assert_eq!(mul_via_axpy(ctx, a, b), mul_model(a, b, p), "mul, p = {p}");
        }
    }
}

#[test]
fn small_primes_agree_on_every_operand_pair() {
    for p in SMALL_PRIMES {
        check_all_pairs(&Small31::new(p), p);
    }
}

fn check_edge_values<F: FieldOps<Coeff = u32>>(ctx: &F, p: u32) {
    let edges = [0, 1, 2, p / 2, p - 2, p - 1];
    for a in edges {
        for b in edges {
            assert_eq!(ctx.add(a, b), add_model(a, b, p), "add {a} {b}, p = {p}");
            assert_eq!(ctx.sub(a, b), sub_model(a, b, p), "sub {a} {b}, p = {p}");
            assert_eq!(
                mul_via_axpy(ctx, a, b),
                mul_model(a, b, p),
                "mul {a} {b}, p = {p}"
            );
        }
    }
}

#[test]
fn large_primes_agree_on_the_edge_values() {
    for p in LARGE_PRIMES {
        check_edge_values(&Small31::new(p), p);
    }
}

#[test]
fn reduce_acc_agrees_on_the_boundary_values() {
    for p in LARGE_PRIMES.iter().chain(SMALL_PRIMES.iter()).copied() {
        let shoup = Small31::new(p);
        let p64 = p as u64;
        let widest = (p64 - 1) + shoup.applications_between_sweeps() * (p64 - 1);
        for x in [
            0,
            1,
            p64 - 1,
            p64,
            p64 + 1,
            2 * p64 - 1,
            2 * p64,
            1 << 32,
            widest,
            u64::MAX,
        ] {
            assert_eq!(shoup.reduce_acc(x), reduce_model(x, p), "p = {p}, x = {x}");
        }
    }
}

#[test]
fn inverse_is_correct_for_every_value_of_a_small_prime() {
    for p in [3u32, 5, 7, 13, 251] {
        let ctx = Small31::new(p);
        for a in 1..p {
            let a_inv = ctx.inv(a);
            assert_eq!(a_inv, inv_model(a, p), "p = {p}, a = {a}");
            assert_eq!(mul_model(a, a_inv, p), 1, "p = {p}, a = {a}");
        }
    }
}

#[test]
fn inverse_is_correct_at_the_benchmark_modulus() {
    let ctx = Small31::new(BENCH_P);
    let mut rng = XorShift(0x0051_71f4_0001);
    let mut values = vec![1, 2, BENCH_P - 1, BENCH_P / 2];
    for _ in 0..2000 {
        values.push(1 + rng.below(BENCH_P as u64 - 1) as u32);
    }
    for a in values {
        let a_inv = ctx.inv(a);
        assert_eq!(a_inv, inv_model(a, BENCH_P), "a = {a}");
        assert_eq!(mul_model(a, a_inv, BENCH_P), 1, "a = {a}");
    }
}

/// Apply `rows` through the kernel and through a per-element model.
///
/// The caller's sweep rule is the one section 3.7 states: count the
/// applications since the scatter or the last sweep, and sweep when the
/// count reaches `applications_between_sweeps`.
fn check_axpy<F: FieldOps<Coeff = u32>>(
    ctx: &F,
    p: u32,
    len: usize,
    rows: &[(Vec<u32>, Vec<u32>)],
    factors: &[u32],
) {
    let mut acc = vec![0u64; len];
    let mut model = vec![0u32; len];
    let mut shoup = Vec::new();
    let mut applied = 0u64;
    for ((cols, vals), &factor) in rows.iter().zip(factors) {
        if applied == ctx.applications_between_sweeps() {
            ctx.sweep(&mut acc);
            applied = 0;
        }
        ctx.precompute(vals, &mut shoup);
        ctx.axpy(&mut acc, cols, vals, &shoup, factor);
        applied += 1;
        for (&c, &w) in cols.iter().zip(vals) {
            let lane = model[c as usize];
            model[c as usize] = add_model(lane, mul_model(factor, w, p), p);
        }
    }
    for (i, &lane) in acc.iter().enumerate() {
        assert_eq!(ctx.reduce_acc(lane), model[i], "lane {i}, p = {p}");
    }
}

/// Random rows over `len` columns, with strictly increasing supports.
fn random_rows(rng: &mut XorShift, p: u32, len: usize, count: usize) -> Vec<(Vec<u32>, Vec<u32>)> {
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for c in 0..len {
            if rng.below(3) != 0 {
                continue;
            }
            cols.push(c as u32);
            vals.push(rng.below(p as u64) as u32);
        }
        rows.push((cols, vals));
    }
    rows
}

#[test]
fn axpy_matches_a_naive_modular_loop() {
    for p in [7u32, 251, 65521, BENCH_P, (1 << 31) - 1] {
        let mut rng = XorShift(0x9e37_79b9_7f4a_7c15 ^ p as u64);
        let len = 64;
        let rows = random_rows(&mut rng, p, len, 200);
        let factors: Vec<u32> = (0..rows.len())
            .map(|_| rng.below(p as u64) as u32)
            .collect();
        check_axpy(&Small31::new(p), p, len, &rows, &factors);
    }
}

#[test]
fn axpy_holds_for_the_extreme_coefficients() {
    let p = BENCH_P;
    let len = 8;
    let rows: Vec<(Vec<u32>, Vec<u32>)> = (0..64)
        .map(|_| ((0..len as u32).collect(), vec![p - 1; len]))
        .collect();
    let factors = vec![p - 1; rows.len()];
    check_axpy(&Small31::new(p), p, len, &rows, &factors);
}

/// A column may not repeat inside one row: the sweep bound counts one
/// addend per lane per application.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "strictly increasing")]
fn repeated_columns_are_rejected() {
    let ctx = Small31::new(BENCH_P);
    let vals = [1u32, 2];
    let mut shoup = Vec::new();
    ctx.precompute(&vals, &mut shoup);
    let mut acc = [0u64; 4];
    ctx.axpy(&mut acc, &[1, 1], &vals, &shoup, 3);
}

/// The sweep reduces every lane, whatever the lane holds.
#[test]
fn the_sweep_reduces_every_lane() {
    let ctx = Small31::new(BENCH_P);
    let p64 = BENCH_P as u64;
    let mut acc = vec![0, 1, p64 - 1, p64, p64 + 7, 2 * p64 - 1, 1 << 40, u64::MAX];
    let want: Vec<u32> = acc.iter().map(|&x| reduce_model(x, BENCH_P)).collect();
    ctx.sweep(&mut acc);
    for (lane, want) in acc.iter().zip(&want) {
        assert_eq!(*lane, *want as u64);
    }
}

/// A debug build catches the caller that does not sweep at the bound.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "needed a sweep")]
fn a_missing_sweep_fails_in_a_debug_build() {
    let ctx = Small31::new(BENCH_P);
    let vals = [BENCH_P - 1];
    let mut shoup = Vec::new();
    ctx.precompute(&vals, &mut shoup);
    // The lane is past the value one more application can hold, which is
    // the state the application count exists to prevent.
    let mut acc = [u64::MAX];
    ctx.axpy(&mut acc, &[0], &vals, &shoup, BENCH_P - 1);
}

#[test]
fn the_shoup_kernel_sweeps_above_two_to_the_thirty_two() {
    for p in LARGE_PRIMES {
        let ctx = Small31::new(p);
        let n = ctx.applications_between_sweeps();
        assert!(n > 1 << 32, "p = {p}, n = {n}");
        let p64 = p as u64;
        let addend = 2 * p64 - 1;
        assert_eq!(n, (u64::MAX - (p64 - 1)) / addend);
        assert!((p64 - 1).checked_add(n * addend).is_some());
        assert!((p64 - 1).checked_add((n + 1) * addend).is_none());
    }
}

/// The bound on a narrow accumulator, where the safe count is small
/// enough to drive to the last application.
#[test]
fn the_sweep_bound_is_exact_on_a_narrow_accumulator() {
    for capacity in [255u64, 4095, u16::MAX as u64] {
        for start in [0u64, 1, 7, 100] {
            for addend in [1u64, 2, 9, 64, 1000] {
                if start > capacity {
                    continue;
                }
                let n = applications_before_sweep(capacity, start, addend);
                let mut lane = start;
                for _ in 0..n {
                    lane += addend;
                    assert!(lane <= capacity, "{capacity} {start} {addend}");
                }
                assert!(lane + addend > capacity, "{capacity} {start} {addend}");
            }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    #[test]
    fn multiply_matches_the_model(
        pi in 0usize..LARGE_PRIMES.len(),
        a in any::<u32>(),
        b in any::<u32>(),
    ) {
        let p = LARGE_PRIMES[pi];
        let (a, b) = (a % p, b % p);
        prop_assert_eq!(mul_via_axpy(&Small31::new(p), a, b), mul_model(a, b, p));
    }

    #[test]
    fn add_and_subtract_match_the_model(
        pi in 0usize..LARGE_PRIMES.len(),
        a in any::<u32>(),
        b in any::<u32>(),
    ) {
        let p = LARGE_PRIMES[pi];
        let (a, b) = (a % p, b % p);
        let ctx = Small31::new(p);
        prop_assert_eq!(ctx.add(a, b), add_model(a, b, p));
        prop_assert_eq!(ctx.sub(a, b), sub_model(a, b, p));
        prop_assert_eq!(ctx.add(a, ctx.sub(0, b)), sub_model(a, b, p));
    }

    #[test]
    fn reduce_acc_matches_the_remainder(
        pi in 0usize..LARGE_PRIMES.len(),
        x in any::<u64>(),
    ) {
        let p = LARGE_PRIMES[pi];
        prop_assert_eq!(Small31::new(p).reduce_acc(x), reduce_model(x, p));
    }

    #[test]
    fn inverse_matches_the_model(pi in 0usize..LARGE_PRIMES.len(), a in any::<u32>()) {
        let p = LARGE_PRIMES[pi];
        let a = 1 + a % (p - 1);
        let ctx = Small31::new(p);
        prop_assert_eq!(ctx.inv(a), inv_model(a, p));
    }

    #[test]
    fn the_sweep_bound_is_exact(
        capacity in 1u64..u64::MAX,
        start in 0u64..1000,
        addend in 1u64..1_000_000,
    ) {
        prop_assume!(start <= capacity);
        let n = applications_before_sweep(capacity, start, addend) as u128;
        let (capacity, start, addend) = (capacity as u128, start as u128, addend as u128);
        prop_assert!(start + n * addend <= capacity);
        prop_assert!(start + (n + 1) * addend > capacity);
    }

    #[test]
    fn axpy_matches_the_model_on_random_rows(
        pi in 0usize..LARGE_PRIMES.len(),
        seed in any::<u64>(),
        rows in 1usize..40,
    ) {
        let p = LARGE_PRIMES[pi];
        let len = 32;
        let mut rng = XorShift(seed | 1);
        let generated = random_rows(&mut rng, p, len, rows);
        let factors: Vec<u32> = (0..rows).map(|_| rng.below(p as u64) as u32).collect();
        check_axpy(&Small31::new(p), p, len, &generated, &factors);
    }
}
