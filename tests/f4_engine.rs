//! The F4 backend against the classic backend.
//!
//! The reduced Gröbner basis of an ideal is unique, so two engines that
//! agree on it agree exactly, coefficients included. Every comparison here
//! is on the whole basis as a set of monic polynomials, never on the
//! leading monomials alone. `tests/differential.rs` runs the same
//! comparison against the independent Buchberger oracle.
//!
//! The families come from `benchmarks/gb-comparison/gen.py`, which states
//! the same definitions.

use std::collections::BTreeMap;
use std::time::{Duration, Instant};

use sylvester::{
    Backend, ComputeError, ComputeOptions, GroebnerBasis, Ideal, PolynomialRing, RingError,
};

/// The benchmark modulus, so the tests run the arithmetic the gate runs.
const P: u64 = 1073741827;

/// One polynomial as a coefficient and one exponent per variable.
type Terms = Vec<(i64, Vec<u16>)>;

fn ring(modulus: u64, nvars: usize) -> PolynomialRing {
    let names: Vec<String> = (0..nvars).map(|index| format!("x{index}")).collect();
    PolynomialRing::prime_field(modulus, names).expect("the modulus is prime")
}

fn ideal_of(ring: &PolynomialRing, system: &[Terms]) -> Ideal {
    let generators: Result<Vec<_>, RingError> = system
        .iter()
        .map(|terms| ring.polynomial(terms.iter().cloned()))
        .collect();
    ring.ideal(generators.expect("the exponents fit the ring"))
        .expect("the polynomials share a ring")
}

fn options(backend: Backend) -> ComputeOptions {
    ComputeOptions::new().backend(backend)
}

fn compute(ideal: &Ideal, backend: Backend) -> GroebnerBasis {
    ideal
        .groebner_basis(options(backend))
        .unwrap_or_else(|error| panic!("{backend:?}: {error}"))
}

/// The basis as one string per element, sorted, so two bases compare as
/// sets.
fn canonical(basis: &GroebnerBasis) -> Vec<String> {
    let mut out: Vec<String> = basis.iter().map(|poly| poly.to_string()).collect();
    out.sort();
    out
}

fn assert_monic(label: &str, basis: &GroebnerBasis) {
    for poly in basis {
        let (coeff, _) = poly.leading_term().expect("a basis element is nonzero");
        assert_eq!(coeff, 1, "{label}: every basis element is monic");
    }
}

/// Computes one system with both backends and compares the bases.
fn backends_agree(label: &str, ideal: &Ideal) -> GroebnerBasis {
    let f4 = compute(ideal, Backend::F4);
    assert_monic(label, &f4);
    let classic = compute(ideal, Backend::Classic);
    assert_eq!(
        canonical(&f4),
        canonical(&classic),
        "{label}: F4 against classic"
    );
    f4
}

/// One polynomial under construction, keyed by exponent vector.
type Build = BTreeMap<Vec<u16>, i64>;

fn add(poly: &mut Build, exps: Vec<u16>, coeff: i64) {
    *poly.entry(exps).or_insert(0) += coeff;
}

fn finish(poly: Build) -> Terms {
    poly.into_iter()
        .filter(|(_, coeff)| *coeff != 0)
        .map(|(exps, coeff)| (coeff, exps))
        .collect()
}

/// cyclic-n in n variables.
fn cyclic(n: usize) -> Vec<Terms> {
    let mut system = Vec::new();
    for d in 1..n {
        let mut poly = Build::new();
        for i in 0..n {
            let mut exps = vec![0u16; n];
            for j in 0..d {
                exps[(i + j) % n] += 1;
            }
            add(&mut poly, exps, 1);
        }
        system.push(finish(poly));
    }
    let mut last = Build::new();
    add(&mut last, vec![1u16; n], 1);
    add(&mut last, vec![0u16; n], -1);
    system.push(finish(last));
    system
}

/// katsura-n in n + 1 variables, the POSSO definition.
fn katsura(n: usize) -> Vec<Terms> {
    let nv = n + 1;
    let mut system = Vec::new();
    let mut first = Build::new();
    for i in 0..=n {
        let mut exps = vec![0u16; nv];
        exps[i] = 1;
        add(&mut first, exps, if i == 0 { 1 } else { 2 });
    }
    add(&mut first, vec![0u16; nv], -1);
    system.push(finish(first));

    let bound = n as i64;
    for m in 0..n as i64 {
        let mut poly = Build::new();
        for i in -bound..=bound {
            let j = m - i;
            if j.abs() > bound {
                continue;
            }
            let mut exps = vec![0u16; nv];
            exps[i.unsigned_abs() as usize] += 1;
            exps[j.unsigned_abs() as usize] += 1;
            add(&mut poly, exps, 1);
        }
        let mut exps = vec![0u16; nv];
        exps[m as usize] = 1;
        add(&mut poly, exps, -1);
        system.push(finish(poly));
    }
    system
}

/// eco-n in n variables.
fn eco(n: usize) -> Vec<Terms> {
    let mut system = Vec::new();
    for k in 1..n {
        let mut poly = Build::new();
        let mut exps = vec![0u16; n];
        exps[k - 1] += 1;
        exps[n - 1] += 1;
        add(&mut poly, exps, 1);
        for i in 1..n - k {
            let mut exps = vec![0u16; n];
            exps[i - 1] += 1;
            exps[i + k - 1] += 1;
            exps[n - 1] += 1;
            add(&mut poly, exps, 1);
        }
        add(&mut poly, vec![0u16; n], -(k as i64));
        system.push(finish(poly));
    }
    let mut last = Build::new();
    for i in 1..n {
        let mut exps = vec![0u16; n];
        exps[i - 1] = 1;
        add(&mut last, exps, 1);
    }
    add(&mut last, vec![0u16; n], 1);
    system.push(finish(last));
    system
}

/// noon-n in n variables.
fn noon(n: usize) -> Vec<Terms> {
    let mut system = Vec::new();
    for i in 0..n {
        let mut poly = Build::new();
        for j in 0..n {
            if j == i {
                continue;
            }
            let mut exps = vec![0u16; n];
            exps[i] += 1;
            exps[j] += 2;
            add(&mut poly, exps, 10);
        }
        let mut exps = vec![0u16; n];
        exps[i] = 1;
        add(&mut poly, exps, -11);
        add(&mut poly, vec![0u16; n], 10);
        system.push(finish(poly));
    }
    system
}

#[test]
fn the_backends_agree_on_the_counterexample_systems() {
    // The two systems of KNOWN_ISSUES.md, which the v0.1 engines got
    // wrong before the repair.
    let f2 = vec![
        vec![(1, vec![3, 0]), (1, vec![0, 3])],
        vec![(1, vec![1, 0]), (1, vec![2, 0]), (1, vec![2, 1])],
        vec![(1, vec![0, 1]), (1, vec![2, 0])],
    ];
    let f3 = vec![
        vec![(1, vec![2, 0, 0]), (1, vec![0, 2, 0])],
        vec![(1, vec![1, 0, 1]), (1, vec![1, 1, 0])],
        vec![(1, vec![0, 1, 0]), (1, vec![1, 1, 0])],
    ];
    let basis = backends_agree("F_2 counterexample", &ideal_of(&ring(2, 2), &f2));
    assert_eq!(canonical(&basis), ["x0", "x1"]);
    backends_agree("F_3 counterexample", &ideal_of(&ring(3, 3), &f3));
}

#[test]
fn the_backends_agree_on_cyclic() {
    for n in 4..=6 {
        backends_agree(&format!("cyclic-{n}"), &ideal_of(&ring(P, n), &cyclic(n)));
    }
}

#[test]
fn the_backends_agree_on_katsura() {
    for n in 4..=7 {
        let system = katsura(n);
        backends_agree(&format!("katsura-{n}"), &ideal_of(&ring(P, n + 1), &system));
    }
}

#[test]
fn the_backends_agree_on_noon() {
    for n in 3..=5 {
        backends_agree(&format!("noon-{n}"), &ideal_of(&ring(P, n), &noon(n)));
    }
}

#[test]
fn the_backends_agree_on_eco() {
    for n in 4..=6 {
        backends_agree(&format!("eco-{n}"), &ideal_of(&ring(P, n), &eco(n)));
    }
}

/// eco-8 costs seconds in a release build and minutes in a debug build.
/// `cargo test --release --test f4_engine -- --ignored` runs it. The
/// comparison is against the classic backend: the legacy batch backend did
/// not finish eco-8 inside 120 s.
#[test]
#[ignore = "slow; run in release"]
fn f4_agrees_with_classic_on_eco_8() {
    let ideal = ideal_of(&ring(P, 8), &eco(8));
    let f4 = compute(&ideal, Backend::F4);
    assert_monic("eco-8", &f4);
    let classic = compute(&ideal, Backend::Classic);
    assert_eq!(canonical(&f4), canonical(&classic), "eco-8");
}

struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

/// Four generators in three variables, degree at most 3.
fn random_system(rng: &mut Rng, p: u64) -> Vec<Terms> {
    (0..4)
        .map(|_| {
            let nterms = 1 + rng.below(4);
            let mut poly = Build::new();
            for _ in 0..nterms {
                let exps: Vec<u16> = (0..3).map(|_| rng.below(4) as u16).collect();
                let coeff = 1 + rng.below(p - 1);
                add(&mut poly, exps, coeff as i64);
            }
            let terms: Build = poly
                .into_iter()
                .map(|(exps, coeff)| (exps, coeff.rem_euclid(p as i64)))
                .collect();
            finish(terms)
        })
        .collect()
}

fn random_systems_agree(p: u64, seeds: u64) {
    let ring = ring(p, 3);
    for seed in 0..seeds {
        let mut rng = Rng(seed.wrapping_mul(0x0123_4567_89AB_CDEF).wrapping_add(p));
        let system = random_system(&mut rng, p);
        if system.iter().all(Vec::is_empty) {
            continue;
        }
        backends_agree(&format!("p={p}, seed={seed}"), &ideal_of(&ring, &system));
    }
}

#[test]
fn the_backends_agree_on_random_systems_over_f2() {
    random_systems_agree(2, 40);
}

#[test]
fn the_backends_agree_on_random_systems_over_f3() {
    random_systems_agree(3, 40);
}

#[test]
fn the_backends_agree_on_random_systems_over_f5() {
    random_systems_agree(5, 40);
}

#[test]
fn the_backends_agree_on_random_systems_over_f7() {
    random_systems_agree(7, 40);
}

#[test]
fn the_backends_agree_on_random_systems_over_f32003() {
    random_systems_agree(32003, 25);
}

#[test]
fn two_runs_write_the_same_bytes() {
    let ring = ring(P, 5);
    let ideal = ideal_of(&ring, &katsura(4));
    let first = compute(&ideal, Backend::F4);
    let second = compute(&ideal, Backend::F4);
    let render = |basis: &GroebnerBasis| {
        basis
            .iter()
            .map(|poly| poly.to_string())
            .collect::<Vec<_>>()
            .join("\n")
    };
    // The element order is part of the output, so this compares more than
    // `canonical` does.
    assert_eq!(render(&first), render(&second));
}

#[test]
fn the_basis_runs_strictly_descending_by_leading_monomial() {
    let ring = ring(P, 5);
    let basis = compute(&ideal_of(&ring, &katsura(4)), Backend::F4);
    let leads: Vec<Vec<u16>> = basis
        .iter()
        .map(|poly| poly.leading_term().expect("nonzero").1.to_vec())
        .collect();
    for pair in leads.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        let (da, db): (u32, u32) = (
            a.iter().map(|&e| u32::from(e)).sum(),
            b.iter().map(|&e| u32::from(e)).sum(),
        );
        let descending = da > db
            || (da == db
                && a.iter()
                    .zip(b)
                    .rev()
                    .find(|(x, y)| x != y)
                    .is_some_and(|(x, y)| x < y));
        assert!(descending, "{a:?} must come before {b:?}");
    }
}

#[test]
fn the_unit_ideal_gives_the_one_element_basis() {
    let ring = ring(P, 2);
    // x0 and x0 + 1 generate the whole ring.
    let system = vec![
        vec![(1, vec![1, 0])],
        vec![(1, vec![1, 0]), (1, vec![0, 0])],
    ];
    let basis = backends_agree("unit ideal", &ideal_of(&ring, &system));
    assert_eq!(canonical(&basis), ["1"]);
}

#[test]
fn a_constant_generator_gives_the_one_element_basis() {
    let ring = ring(P, 2);
    let system = vec![
        vec![(1, vec![2, 0]), (1, vec![0, 1])],
        vec![(7, vec![0, 0])],
    ];
    let basis = backends_agree("constant generator", &ideal_of(&ring, &system));
    assert_eq!(canonical(&basis), ["1"]);
}

#[test]
fn the_zero_ideal_gives_the_empty_basis() {
    let ring = ring(P, 3);
    let zero: Vec<Terms> = vec![Vec::new(), Vec::new()];
    let basis = backends_agree("zero generators", &ideal_of(&ring, &zero));
    assert!(basis.is_empty());
}

#[test]
fn no_generator_gives_the_empty_basis() {
    let ring = ring(P, 3);
    let basis = backends_agree("no generator", &ideal_of(&ring, &[]));
    assert!(basis.is_empty());
}

#[test]
fn one_generator_gives_that_generator_monic() {
    let ring = ring(P, 2);
    let system = vec![vec![(5, vec![2, 0]), (3, vec![0, 1])]];
    let basis = backends_agree("one generator", &ideal_of(&ring, &system));
    assert_eq!(basis.len(), 1);
    // 5*x0^2 + 3*x1 scaled by the inverse of 5.
    let expected = ring
        .polynomial([
            (1i64, vec![2u16, 0]),
            (3 * inverse(5, P) as i64, vec![0, 1]),
        ])
        .expect("the exponents fit the ring");
    assert_eq!(basis[0], expected);
}

#[test]
fn duplicate_generators_change_nothing() {
    let ring = ring(P, 3);
    let once = katsura(2);
    let mut twice = once.clone();
    twice.extend(once.iter().cloned());
    let single = backends_agree("katsura-2", &ideal_of(&ring, &once));
    let doubled = backends_agree("katsura-2 twice", &ideal_of(&ring, &twice));
    assert_eq!(canonical(&single), canonical(&doubled));
}

/// The inverse of `a` modulo the prime `p`, by Fermat.
fn inverse(a: u64, p: u64) -> u64 {
    let mut acc = 1u64;
    let mut base = a % p;
    let mut exp = p - 2;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = (acc as u128 * base as u128 % p as u128) as u64;
        }
        base = (base as u128 * base as u128 % p as u128) as u64;
        exp >>= 1;
    }
    acc
}

#[test]
fn an_exhausted_deadline_stops_the_run() {
    let ring = ring(P, 5);
    let ideal = ideal_of(&ring, &katsura(4));
    assert_eq!(
        ideal.groebner_basis(options(Backend::F4).timeout(Duration::ZERO)),
        Err(ComputeError::Timeout)
    );
    // A run with no pair still reports the budget rather than a basis.
    let single = ideal_of(&ring, &[vec![(1, vec![2, 0, 0, 0, 0])]]);
    assert_eq!(
        single.groebner_basis(options(Backend::F4).timeout(Duration::ZERO)),
        Err(ComputeError::Timeout)
    );
    let empty = ideal_of(&ring, &[]);
    assert_eq!(
        empty.groebner_basis(options(Backend::F4).timeout(Duration::ZERO)),
        Err(ComputeError::Timeout)
    );
}

#[test]
fn a_deadline_well_under_the_full_cost_stops_the_run() {
    // Self-calibrating: measure the free run, then ask for a twentieth of
    // it. The engine checks the clock between batches and inside the
    // symbolic, update, and row loops, so it must stop early.
    let ring = ring(P, 8);
    let ideal = ideal_of(&ring, &katsura(7));
    let start = Instant::now();
    ideal
        .groebner_basis(options(Backend::F4))
        .expect("the free run finishes");
    let full = start.elapsed();

    let start = Instant::now();
    let error = ideal
        .groebner_basis(options(Backend::F4).timeout(full / 20))
        .expect_err("the deadline stops the run");
    let elapsed = start.elapsed();
    assert_eq!(error, ComputeError::Timeout);
    assert!(
        elapsed < full / 2,
        "the run must stop at the deadline: {elapsed:?} of a full {full:?}"
    );
}

#[test]
fn a_timeout_too_large_to_represent_leaves_no_deadline() {
    let ring = ring(P, 4);
    let ideal = ideal_of(&ring, &cyclic(4));
    let basis = ideal
        .groebner_basis(options(Backend::F4).timeout(Duration::MAX))
        .expect("a timeout too large to add still lets the run finish");
    assert!(!basis.is_empty());
}

#[test]
fn an_exhausted_memory_budget_stops_the_run() {
    // Self-calibrating: the free run gives the basis, then the limit
    // shrinks until the run stops. The engine holds the basis, the pair
    // queue, the three monomial tables, and the batch, so a limit under
    // the size of one batch cannot hold it.
    let ring = ring(P, 8);
    let ideal = ideal_of(&ring, &katsura(7));
    let free = ideal
        .groebner_basis(options(Backend::F4))
        .expect("the free run finishes");
    let generous = ideal
        .groebner_basis(options(Backend::F4).memory_limit(1 << 30))
        .expect("a generous limit holds the run");
    assert_eq!(canonical(&free), canonical(&generous));

    for limit in [0, 1 << 8, 1 << 10] {
        assert_eq!(
            ideal.groebner_basis(options(Backend::F4).memory_limit(limit)),
            Err(ComputeError::MemoryLimitExceeded),
            "limit {limit}"
        );
    }
}

#[test]
fn an_exponent_past_the_width_is_a_typed_limit() {
    // The F4 engine bounds one exponent, not the degree of a pair's lcm.
    // Here the row `x1^10000 * f` needs the exponent 70000 on x1, which is
    // past the 65535 an exponent holds.
    let ring = ring(7, 2);
    let system = vec![
        vec![(1, vec![60000, 0]), (1, vec![0, 60000])],
        vec![(1, vec![30000, 10000]), (1, vec![0, 1])],
    ];
    assert_eq!(
        ideal_of(&ring, &system).groebner_basis(options(Backend::F4)),
        Err(ComputeError::ExponentLimit { limit: 65535 })
    );
}

#[test]
fn a_run_that_outgrows_the_narrow_lanes_gives_the_same_basis() {
    // Exponents near 100 fit the 127 bound of the eight-bit packing, and
    // their products do not. The run restarts at the wider packing and
    // must land on the same basis as the classic backend, which uses its
    // own owned-exponent arithmetic.
    let ring = ring(P, 2);
    let system = vec![
        vec![(1, vec![100, 30]), (1, vec![0, 1]), (3, vec![0, 0])],
        vec![(1, vec![30, 100]), (1, vec![1, 0]), (5, vec![0, 0])],
    ];
    backends_agree("wide exponents", &ideal_of(&ring, &system));
}

#[test]
fn the_report_counts_what_the_f4_run_did() {
    let ideal = ideal_of(&ring(P, 6), &cyclic(6));
    let (basis, report) = ideal
        .groebner_basis_with_report(options(Backend::F4))
        .expect("no limit");
    assert_eq!(canonical(&basis), canonical(&compute(&ideal, Backend::F4)));
    assert_eq!(report.backend, Backend::F4);
    let counters = report.counters.expect("the F4 backend counts");
    assert!(counters.batches > 0);
    assert!(counters.matrix_rows > 0);
    assert!(counters.matrix_nonzeros > 0);
    assert!(counters.new_pivots > 0);
    assert!(counters.pairs_generated > 0);
    assert_eq!(counters.lane_restarts, 0, "cyclic-6 fits 8 lanes");
    assert!(report.elapsed > Duration::ZERO);
}

#[test]
fn the_classic_backend_reports_no_counters() {
    let ideal = ideal_of(&ring(P, 4), &cyclic(4));
    let (basis, report) = ideal
        .groebner_basis_with_report(options(Backend::Classic))
        .expect("no limit");
    assert_eq!(
        canonical(&basis),
        canonical(&compute(&ideal, Backend::Classic))
    );
    assert_eq!(report.backend, Backend::Classic);
    assert_eq!(report.counters, None);
}

/// The output must not depend on the thread count (design section 5).
///
/// `ComputeOptions::threads` puts the whole run in a pool of that size,
/// and the kernel's parallel phase reduces the rows of a batch against
/// the frozen pivot set on that pool. katsura-7 and cyclic-6 both have
/// batches above the work threshold, so the parallel phase runs.
#[test]
fn the_thread_count_does_not_change_the_basis() {
    for (nvars, system) in [(8, katsura(7)), (6, cyclic(6))] {
        let ring = ring(P, nvars);
        let ideal = ideal_of(&ring, &system);
        let one = ideal
            .groebner_basis(options(Backend::F4).threads(1))
            .expect("no limit");
        for threads in [2, 8, 16] {
            let many = ideal
                .groebner_basis(options(Backend::F4).threads(threads))
                .expect("no limit");
            assert_eq!(
                canonical(&one),
                canonical(&many),
                "{threads} threads changed the basis"
            );
        }
    }
}

/// The report carries the threads the run had, not the threads it asked
/// for.
#[test]
fn the_report_carries_the_threads_the_run_had() {
    let ideal = ideal_of(&ring(P, 4), &cyclic(4));
    let asked = ideal
        .groebner_basis_with_report(options(Backend::F4).threads(3))
        .expect("no limit")
        .1;
    assert_eq!(asked.threads_used, 3);

    let one = ideal
        .groebner_basis_with_report(options(Backend::F4).threads(1))
        .expect("no limit")
        .1;
    assert_eq!(one.threads_used, 1);

    let classic = ideal
        .groebner_basis_with_report(options(Backend::Classic))
        .expect("no limit")
        .1;
    assert_eq!(classic.threads_used, 1, "the classic backend is sequential");

    let pool = ideal
        .groebner_basis_with_report(options(Backend::F4))
        .expect("no limit")
        .1;
    assert!(pool.threads_used >= 1, "the global pool holds a thread");
}

/// A run whose input alone passes the limit stops before its first batch.
///
/// The engine holds the interned generators and seeded basis before its
/// first batch, so one byte stops every nonempty input. A literally empty
/// input holds no engine data and succeeds under the hardened budget.
#[test]
fn a_limit_the_input_alone_passes_stops_the_run() {
    let ring = ring(P, 2);
    let unit = ideal_of(&ring, &[vec![(1, vec![0, 0])]]);
    assert_eq!(
        unit.groebner_basis(options(Backend::F4).memory_limit(1)),
        Err(ComputeError::MemoryLimitExceeded),
        "a constant generator"
    );
    assert_eq!(
        unit.groebner_basis(options(Backend::F4)),
        ideal_of(&ring, &[vec![(1, vec![0, 0])]]).groebner_basis(options(Backend::F4)),
        "the same input without a limit still finishes"
    );

    let empty = ideal_of(&ring, &[]);
    assert!(
        empty
            .groebner_basis(options(Backend::F4).memory_limit(0))
            .expect("the empty input holds no engine data")
            .is_empty()
    );

    let one_term = ideal_of(&ring, &[vec![(1, vec![3, 1])]]);
    assert_eq!(
        one_term.groebner_basis(options(Backend::F4).memory_limit(1)),
        Err(ComputeError::MemoryLimitExceeded),
        "one generator"
    );
}

/// The smallest memory limit `ideal` finishes under, by binary search.
///
/// The search assumes what the engine holds to: a run that finishes under
/// one limit finishes under every larger one.
fn smallest_memory_limit(ideal: &Ideal) -> usize {
    let succeeds = |limit: usize| {
        ideal
            .groebner_basis(options(Backend::F4).memory_limit(limit))
            .is_ok()
    };
    let mut high = 1 << 12;
    while !succeeds(high) {
        high *= 2;
        assert!(high <= 1 << 32, "no limit under 4 GiB finishes the run");
    }
    let mut low = 0;
    while high - low > 1 {
        let middle = low + (high - low) / 2;
        if succeeds(middle) {
            high = middle;
        } else {
            low = middle;
        }
    }
    high
}

/// The limit binds where the run peaks.
///
/// Self-calibrating: the search finds the smallest limit the run finishes
/// under. One byte less is a typed stop, and the limit itself gives the
/// basis the free run gives.
#[test]
fn a_limit_just_under_the_peak_stops_the_run_and_one_above_it_does_not() {
    let ring = ring(P, 5);
    let ideal = ideal_of(&ring, &katsura(4));
    let free = ideal
        .groebner_basis(options(Backend::F4))
        .expect("the free run finishes");

    let peak = smallest_memory_limit(&ideal);
    let held = ideal
        .groebner_basis(options(Backend::F4).memory_limit(peak))
        .expect("the peak limit holds the run");
    assert_eq!(canonical(&free), canonical(&held));
    assert_eq!(
        ideal.groebner_basis(options(Backend::F4).memory_limit(peak - 1)),
        Err(ComputeError::MemoryLimitExceeded),
        "one byte under the peak"
    );
}

/// A batch too large for the limit succeeds after the retry halves it.
///
/// Self-calibrating: the search finds the smallest limit the run finishes
/// under. That limit is under the peak of the whole pair batch, so the run
/// only finishes because a retry built a smaller batch. The retry must
/// release what the failed attempt held, or the estimate cannot fall and
/// the loop runs down to one pair and stops.
#[test]
fn a_batch_over_the_limit_finishes_after_the_retry_halves_it() {
    let ring = ring(P, 7);
    let ideal = ideal_of(&ring, &katsura(6));
    let free = ideal
        .groebner_basis(options(Backend::F4))
        .expect("the free run finishes");

    let peak = smallest_memory_limit(&ideal);
    let (basis, report) = ideal
        .groebner_basis_with_report(options(Backend::F4).memory_limit(peak))
        .expect("the peak limit holds the run");
    assert_eq!(canonical(&free), canonical(&basis));
    let counters = report.counters.expect("the F4 backend counts");
    assert!(
        counters.batch_retries > 0,
        "the smallest limit the run finishes under splits a batch: {counters:?}"
    );
}

/// A deadline inside one batch stops there, not after it.
///
/// Self-calibrating: the free run gives the cost, and the timeout is a
/// twentieth of it, which lands inside katsura-8's largest batch. The
/// engine reads the clock inside the elimination walk and once per row of
/// the parallel phase, so the run stops within twice the timeout at either
/// thread count.
#[test]
fn a_deadline_inside_one_batch_stops_within_twice_the_timeout() {
    let ring = ring(P, 9);
    let ideal = ideal_of(&ring, &katsura(8));
    for threads in [1, 4] {
        let start = Instant::now();
        ideal
            .groebner_basis(options(Backend::F4).threads(threads))
            .expect("the free run finishes");
        let full = start.elapsed();
        let timeout = full / 20;

        let start = Instant::now();
        let error = ideal
            .groebner_basis(options(Backend::F4).threads(threads).timeout(timeout))
            .expect_err("the deadline stops the run");
        let elapsed = start.elapsed();
        assert_eq!(error, ComputeError::Timeout);
        assert!(
            elapsed < timeout * 2,
            "{threads} threads: stopped after {elapsed:?} on a {timeout:?} timeout"
        );
    }
}
