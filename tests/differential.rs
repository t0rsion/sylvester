//! Differential tests: every backend against an independent Buchberger
//! oracle on random systems over F_2, F_3, F_5, F_7 (a larger shape), and
//! F_32003, the known counterexample systems in all generator orders, a
//! mutation-derived regression fixture, and (`#[ignore]`d) an exhaustive
//! small-system sweep and a degree-reversal stress sweep.
//!
//! The oracle below is self-contained (own polynomial arithmetic, grevlex
//! from the definition, textbook Buchberger completion plus interreduction);
//! it shares no code with `src/**`, so it cannot inherit an engine bug. The
//! comparison is exact: the reduced Gröbner basis of an ideal is unique, so
//! each engine's output must equal the oracle's reduced basis as a set of
//! monic polynomials.
//!
//! Seeds are deterministic; failures print the system and both bases.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use sylvester::{Backend, ComputeOptions, GroebnerBasis, Ideal, PolynomialRing};

fn sub_mod(a: u64, b: u64, p: u64) -> u64 {
    (a + p - b % p) % p
}

fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn inv_mod(a: u64, p: u64) -> u64 {
    assert!(!a.is_multiple_of(p), "attempted to invert zero mod {p}");
    let mut base = a % p;
    let mut exp = p - 2;
    let mut acc = 1u64;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base, p);
        }
        base = mul_mod(base, base, p);
        exp >>= 1;
    }
    acc
}

type Exps = Vec<u16>;

fn grevlex_cmp(a: &[u16], b: &[u16]) -> Ordering {
    assert_eq!(a.len(), b.len(), "monomials from different rings");
    let da: u64 = a.iter().map(|&e| e as u64).sum();
    let db: u64 = b.iter().map(|&e| e as u64).sum();
    match da.cmp(&db) {
        Ordering::Equal => {}
        unequal => return unequal,
    }
    for (ea, eb) in a.iter().zip(b.iter()).rev() {
        match ea.cmp(eb) {
            Ordering::Equal => continue,
            Ordering::Less => return Ordering::Greater,
            Ordering::Greater => return Ordering::Less,
        }
    }
    Ordering::Equal
}

fn divides(divisor: &[u16], multiple: &[u16]) -> bool {
    divisor.iter().zip(multiple.iter()).all(|(d, m)| d <= m)
}

fn exps_sub(a: &[u16], b: &[u16]) -> Exps {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}

fn exps_add(a: &[u16], b: &[u16]) -> Exps {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

fn exps_lcm(a: &[u16], b: &[u16]) -> Exps {
    a.iter().zip(b.iter()).map(|(x, y)| *x.max(y)).collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Poly {
    terms: BTreeMap<Exps, u64>,
}

impl Poly {
    fn zero() -> Self {
        Poly {
            terms: BTreeMap::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn leading_term(&self) -> Option<(&Exps, u64)> {
        self.terms
            .iter()
            .max_by(|(a, _), (b, _)| grevlex_cmp(a, b))
            .map(|(m, &c)| (m, c))
    }

    /// self -= factor * x^shift * g (mod p).
    fn sub_scaled(&mut self, factor: u64, shift: &[u16], g: &Poly, p: u64) {
        for (mono, &coeff) in &g.terms {
            let target = exps_add(mono, shift);
            let delta = mul_mod(factor, coeff, p);
            let entry = self.terms.entry(target).or_insert(0);
            *entry = sub_mod(*entry, delta, p);
        }
        self.terms.retain(|_, c| *c != 0);
    }

    fn make_monic(&self, p: u64) -> Poly {
        let Some((_, lead)) = self.leading_term() else {
            return Poly::zero();
        };
        let scale = inv_mod(lead, p);
        let mut out = Poly::zero();
        for (exps, &coeff) in &self.terms {
            out.terms.insert(exps.clone(), mul_mod(coeff, scale, p));
        }
        out
    }
}

fn normal_form(f: &Poly, basis: &[Poly], p: u64) -> Poly {
    let mut work = f.clone();
    let mut remainder = Poly::zero();
    while let Some((lead_mono, lead_coeff)) = work.leading_term() {
        let lead_mono = lead_mono.clone();
        let reducer = basis.iter().find_map(|g| {
            let (g_mono, g_coeff) = g.leading_term()?;
            divides(g_mono, &lead_mono).then_some((g, g_mono.clone(), g_coeff))
        });
        match reducer {
            Some((g, g_mono, g_coeff)) => {
                let shift = exps_sub(&lead_mono, &g_mono);
                let factor = mul_mod(lead_coeff, inv_mod(g_coeff, p), p);
                work.sub_scaled(factor, &shift, g, p);
            }
            None => {
                remainder.terms.insert(lead_mono.clone(), lead_coeff);
                work.terms.remove(&lead_mono);
            }
        }
    }
    remainder
}

fn s_polynomial(f: &Poly, g: &Poly, p: u64) -> Poly {
    let (f_mono, f_coeff) = f.leading_term().expect("s_polynomial: f is zero");
    let (g_mono, g_coeff) = g.leading_term().expect("s_polynomial: g is zero");
    let lcm = exps_lcm(f_mono, g_mono);
    let f_shift = exps_sub(&lcm, f_mono);
    let g_shift = exps_sub(&lcm, g_mono);
    let mut s = Poly::zero();
    s.sub_scaled(sub_mod(0, inv_mod(f_coeff, p), p), &f_shift, f, p);
    s.sub_scaled(inv_mod(g_coeff, p), &g_shift, g, p);
    s
}

/// Textbook Buchberger completion: criterion-free worklist of all pairs,
/// each processed exactly once, new pairs enqueued on every insertion.
/// Pairs are processed lowest lcm degree first (the normal selection
/// strategy); selection order does not affect correctness, only size.
fn buchberger(generators: &[Poly], p: u64) -> Vec<Poly> {
    fn lcm_degree(basis: &[Poly], i: usize, j: usize) -> u64 {
        let (mi, _) = basis[i].leading_term().expect("nonzero");
        let (mj, _) = basis[j].leading_term().expect("nonzero");
        exps_lcm(mi, mj).iter().map(|&e| e as u64).sum()
    }

    let mut basis: Vec<Poly> = generators
        .iter()
        .filter(|f| !f.is_zero())
        .cloned()
        .collect();
    let mut pairs: BTreeMap<u64, Vec<(usize, usize)>> = BTreeMap::new();
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            pairs
                .entry(lcm_degree(&basis, i, j))
                .or_default()
                .push((i, j));
        }
    }
    loop {
        let Some((&deg, _)) = pairs.iter().next() else {
            return basis;
        };
        let bucket = pairs.get_mut(&deg).expect("bucket exists");
        let (i, j) = bucket.pop().expect("bucket non-empty");
        if bucket.is_empty() {
            pairs.remove(&deg);
        }
        let s = s_polynomial(&basis[i], &basis[j], p);
        let r = normal_form(&s, &basis, p);
        if !r.is_zero() {
            let k = basis.len();
            basis.push(r.make_monic(p));
            for i in 0..k {
                pairs
                    .entry(lcm_degree(&basis, i, k))
                    .or_default()
                    .push((i, k));
            }
        }
    }
}

/// Reduced Gröbner basis from a Gröbner basis: minimalize leading terms,
/// then interreduce until stable.
fn reduce_basis(basis: &[Poly], p: u64) -> Vec<Poly> {
    let mut current: Vec<Poly> = basis
        .iter()
        .filter(|f| !f.is_zero())
        .map(|f| f.make_monic(p))
        .collect();
    let mut minimal: Vec<Poly> = Vec::new();
    for (i, f) in current.iter().enumerate() {
        let (lm_f, _) = f.leading_term().expect("nonzero");
        let redundant = current.iter().enumerate().any(|(j, g)| {
            if i == j {
                return false;
            }
            let (lm_g, _) = g.leading_term().expect("nonzero");
            divides(lm_g, lm_f) && (lm_g != lm_f || j < i)
        });
        if !redundant {
            minimal.push(f.clone());
        }
    }
    current = minimal;
    loop {
        let mut next: Vec<Poly> = Vec::new();
        for i in 0..current.len() {
            let others: Vec<Poly> = current
                .iter()
                .enumerate()
                .filter(|&(j, _)| i != j)
                .map(|(_, g)| g.clone())
                .collect();
            let r = normal_form(&current[i], &others, p);
            if !r.is_zero() {
                next.push(r.make_monic(p));
            }
        }
        if next == current {
            return next;
        }
        current = next;
    }
}

/// Splitmix64.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

#[derive(Clone, Copy)]
struct Shape {
    nvars: usize,
    max_deg: u64,
    ngens: usize,
    min_terms: u64,
    term_spread: u64,
}

/// 3 variables, 4 generators, degree <= 3, 1..=4 term-insertion attempts.
const SMALL: Shape = Shape {
    nvars: 3,
    max_deg: 3,
    ngens: 4,
    min_terms: 1,
    term_spread: 4,
};

/// 4 variables, 5 generators, degree <= 3, 2..=5 term-insertion attempts.
const HARDER: Shape = Shape {
    nvars: 4,
    max_deg: 3,
    ngens: 5,
    min_terms: 2,
    term_spread: 4,
};

/// A random monomial in `shape.nvars` variables of total degree <= max_deg.
fn random_mono(rng: &mut Rng, shape: Shape) -> Exps {
    loop {
        let exps: Exps = (0..shape.nvars)
            .map(|_| rng.below(shape.max_deg + 1) as u16)
            .collect();
        let deg: u64 = exps.iter().map(|&e| e as u64).sum();
        if deg <= shape.max_deg {
            return exps;
        }
    }
}

/// A random system of `shape.ngens` polynomials; term counts are insertion
/// attempts, so colliding monomials can leave fewer (or zero) terms.
fn random_system(rng: &mut Rng, p: u64, shape: Shape) -> Vec<Poly> {
    (0..shape.ngens)
        .map(|_| {
            let nterms = shape.min_terms + rng.below(shape.term_spread);
            let mut poly = Poly::zero();
            for _ in 0..nterms {
                let mono = random_mono(rng, shape);
                let coeff = 1 + rng.below(p - 1);
                let entry = poly.terms.entry(mono).or_insert(0);
                *entry = (*entry + coeff) % p;
            }
            poly.terms.retain(|_, c| *c != 0);
            poly
        })
        .collect()
}

/// The variable names the engine ring uses, one per oracle variable.
fn variable_names(nvars: usize) -> Vec<String> {
    (0..nvars).map(|index| format!("x{index}")).collect()
}

fn to_engine_ideal(system: &[Poly], p: u64, nvars: usize) -> Ideal {
    let ring = PolynomialRing::prime_field(p, variable_names(nvars)).expect("the modulus is prime");
    let generators: Vec<_> = system
        .iter()
        .map(|poly| {
            let terms: Vec<(i64, Vec<u16>)> = poly
                .terms
                .iter()
                .map(|(exps, &coeff)| {
                    let mut exps = exps.clone();
                    exps.resize(nvars, 0);
                    (coeff as i64, exps)
                })
                .collect();
            ring.polynomial(terms).expect("the exponents fit the ring")
        })
        .collect();
    ring.ideal(generators)
        .expect("the polynomials share a ring")
}

fn engine_output_to_polys(output: &GroebnerBasis, p: u64, nvars: usize) -> Vec<Poly> {
    output
        .iter()
        .map(|poly| {
            let mut out = Poly::zero();
            for (coeff, exps) in poly.terms() {
                assert!(exps.len() <= nvars, "engine output has too many variables");
                let mut exps = exps.to_vec();
                exps.resize(nvars, 0);
                let entry = out.terms.entry(exps).or_insert(0);
                *entry = (*entry + coeff.value() % p) % p;
            }
            out.terms.retain(|_, c| *c != 0);
            out
        })
        .collect()
}

fn canonical_set(basis: &[Poly], p: u64) -> BTreeSet<Vec<(Exps, u64)>> {
    basis
        .iter()
        .filter(|f| !f.is_zero())
        .map(|f| make_canonical(&f.make_monic(p)))
        .collect()
}

fn make_canonical(f: &Poly) -> Vec<(Exps, u64)> {
    f.terms.iter().map(|(e, &c)| (e.clone(), c)).collect()
}

fn render_system(system: &[Poly]) -> String {
    let rendered: Vec<String> = system.iter().map(|f| format!("{:?}", f.terms)).collect();
    format!("[{}]", rendered.join(", "))
}

fn check_backend(
    label: &str,
    system: &[Poly],
    oracle_reduced: &BTreeSet<Vec<(Exps, u64)>>,
    output: Result<GroebnerBasis, sylvester::ComputeError>,
    p: u64,
    nvars: usize,
    seed: u64,
) {
    let output = output.unwrap_or_else(|e| {
        panic!(
            "{label}: engine error {e} on seed {seed}, p={p}, system {}",
            render_system(system)
        )
    });
    let candidate = engine_output_to_polys(&output, p, nvars);
    for f in &candidate {
        assert!(
            !f.is_zero(),
            "{label}: zero polynomial in output on seed {seed}, p={p}, system {}",
            render_system(system)
        );
    }
    let candidate_set = canonical_set(&candidate, p);
    assert_eq!(
        candidate_set.len(),
        candidate.len(),
        "{label}: duplicate polynomial in output on seed {seed}, p={p}, system {}",
        render_system(system)
    );
    assert_eq!(
        &candidate_set,
        oracle_reduced,
        "{label}: mismatch on seed {seed}, p={p}\nsystem: {}\nengine: {:?}\noracle: {:?}",
        render_system(system),
        candidate_set,
        oracle_reduced
    );
}

/// Run both backends on one system and compare each against the oracle's
/// reduced basis.
fn check_system(system: &[Poly], p: u64, nvars: usize, seed: u64) {
    if system.iter().all(Poly::is_zero) {
        return;
    }
    let oracle = buchberger(system, p);
    let oracle_reduced = canonical_set(&reduce_basis(&oracle, p), p);
    let ideal = to_engine_ideal(system, p, nvars);

    for (label, backend) in [("f4", Backend::F4), ("classic", Backend::Classic)] {
        check_backend(
            label,
            system,
            &oracle_reduced,
            ideal.groebner_basis(ComputeOptions::new().backend(backend)),
            p,
            nvars,
            seed,
        );
    }
}

fn run_differential(p: u64, seeds: std::ops::Range<u64>, shape: Shape) {
    for seed in seeds {
        let mut rng = Rng(seed.wrapping_mul(0x0123_4567_89AB_CDEF).wrapping_add(p));
        let system = random_system(&mut rng, p, shape);
        check_system(&system, p, shape.nvars, seed);
    }
}

#[test]
fn differential_f2() {
    run_differential(2, 0..100, SMALL);
}

#[test]
fn differential_f3() {
    run_differential(3, 0..100, SMALL);
}

#[test]
fn differential_f5() {
    run_differential(5, 0..100, SMALL);
}

#[test]
fn differential_f7_harder_shape() {
    run_differential(7, 0..25, HARDER);
}

/// A prime past 2^15, so the coefficient arithmetic runs on wide values.
#[test]
fn differential_f32003() {
    run_differential(32003, 0..25, SMALL);
}

fn poly_from_terms(terms: &[(u64, &[u16])], p: u64) -> Poly {
    let mut poly = Poly::zero();
    for &(coeff, exps) in terms {
        let entry = poly.terms.entry(exps.to_vec()).or_insert(0);
        *entry = (*entry + coeff) % p;
    }
    poly.terms.retain(|_, c| *c != 0);
    poly
}

fn check_all_generator_orders(system: &[Poly], p: u64, nvars: usize) {
    let n = system.len();
    let mut indices: Vec<usize> = (0..n).collect();
    let mut orders: Vec<Vec<usize>> = Vec::new();
    permutations(&mut indices, 0, &mut orders);
    for (k, order) in orders.iter().enumerate() {
        let permuted: Vec<Poly> = order.iter().map(|&i| system[i].clone()).collect();
        check_system(&permuted, p, nvars, k as u64);
    }
}

fn permutations(items: &mut Vec<usize>, start: usize, out: &mut Vec<Vec<usize>>) {
    if start == items.len() {
        out.push(items.clone());
        return;
    }
    for i in start..items.len() {
        items.swap(start, i);
        permutations(items, start + 1, out);
        items.swap(start, i);
    }
}

/// The two hand-verified counterexample systems from KNOWN_ISSUES.md, in
/// every generator order: the incremental outer loop and POT signatures make
/// engine behavior order-sensitive, the ideal (and thus the oracle) is not.
#[test]
fn counterexample_systems_in_all_generator_orders() {
    let f2_system = vec![
        poly_from_terms(&[(1, &[3, 0]), (1, &[0, 3])], 2),
        poly_from_terms(&[(1, &[1, 0]), (1, &[2, 0]), (1, &[2, 1])], 2),
        poly_from_terms(&[(1, &[0, 1]), (1, &[2, 0])], 2),
    ];
    check_all_generator_orders(&f2_system, 2, 2);

    let f3_system = vec![
        poly_from_terms(&[(1, &[2, 0, 0]), (1, &[0, 2, 0])], 3),
        poly_from_terms(&[(1, &[1, 0, 1]), (1, &[1, 1, 0])], 3),
        poly_from_terms(&[(1, &[0, 1, 0]), (1, &[1, 1, 0])], 3),
    ];
    check_all_generator_orders(&f3_system, 3, 3);
}

/// Regression fixture for the rewritten criterion: reverting canonical-
/// rewriter substitution to pair deletion makes the matrix backend return a
/// basis of the wrong ideal on exactly this F_3 system (found by mutation
/// testing; batch_reversal_f3 seed 559).
#[test]
fn rewrite_substitution_regression_f3() {
    let system = vec![
        poly_from_terms(&[(1, &[0, 1, 0]), (1, &[4, 0, 0])], 3),
        poly_from_terms(&[(2, &[1, 0, 0]), (2, &[2, 0, 2])], 3),
        poly_from_terms(&[(2, &[0, 0, 1]), (2, &[1, 0, 0]), (1, &[1, 2, 1])], 3),
        poly_from_terms(&[(2, &[0, 0, 0]), (2, &[1, 0, 0]), (2, &[2, 1, 1])], 3),
    ];
    check_system(&system, 3, 3, 559);
}

/// All monomials and binomials over the 10 monomials of degree <= 2 in three
/// variables: 55 polynomials, all 55^3 ordered triples, exhaustively, F_2.
/// Slow; run explicitly with --ignored (release recommended).
#[test]
#[ignore = "exhaustive sweep; minutes of runtime, run with --ignored"]
fn exhaustive_f2_monomial_binomial_triples() {
    let monos = degree_two_monomials();
    assert_eq!(monos.len(), 10);
    let polys = monomials_and_binomials(&monos);
    assert_eq!(polys.len(), 55);
    check_ordered_triples(&polys);
}

fn degree_two_monomials() -> Vec<Exps> {
    let mut monos: Vec<Exps> = Vec::new();
    for a in 0..=2u16 {
        for b in 0..=2u16 {
            for c in 0..=2u16 {
                if a + b + c <= 2 {
                    monos.push(vec![a, b, c]);
                }
            }
        }
    }
    monos
}

fn monomials_and_binomials(monos: &[Exps]) -> Vec<Poly> {
    let mut polys: Vec<Poly> = Vec::new();
    for m in monos {
        polys.push(poly_from_terms(&[(1, m)], 2));
    }
    for i in 0..monos.len() {
        for j in (i + 1)..monos.len() {
            polys.push(poly_from_terms(&[(1, &monos[i]), (1, &monos[j])], 2));
        }
    }
    polys
}

fn check_ordered_triples(polys: &[Poly]) {
    for i in 0..polys.len() {
        for j in 0..polys.len() {
            for k in 0..polys.len() {
                let system = vec![polys[i].clone(), polys[j].clone(), polys[k].clone()];
                let seed = ((i * 55 + j) * 55 + k) as u64;
                check_system(&system, 2, 3, seed);
            }
        }
    }
}

/// Degree-reversal stress over F_3: each generator has one degree-4 leading
/// term plus a few terms of degree <= 1, so high-degree S-pairs generate
/// lower-degree rows and pairs across matrix degree batches.
#[test]
#[ignore = "sweep; run with --ignored"]
fn batch_reversal_f3() {
    for seed in 0..1000u64 {
        let mut rng = Rng(seed.wrapping_mul(0xD1B5_4A32_D192_ED03).wrapping_add(3));
        let system: Vec<Poly> = (0..4)
            .map(|_| {
                let lead = loop {
                    let exps: Exps = (0..3).map(|_| rng.below(5) as u16).collect();
                    if exps.iter().map(|&e| e as u64).sum::<u64>() == 4 {
                        break exps;
                    }
                };
                let mut terms: Vec<(u64, Exps)> = vec![(1 + rng.below(2), lead)];
                for _ in 0..(1 + rng.below(3)) {
                    let mono: Exps = loop {
                        let exps: Exps = (0..3).map(|_| rng.below(2) as u16).collect();
                        if exps.iter().map(|&e| e as u64).sum::<u64>() <= 1 {
                            break exps;
                        }
                    };
                    terms.push((1 + rng.below(2), mono));
                }
                let refs: Vec<(u64, &[u16])> =
                    terms.iter().map(|(c, e)| (*c, e.as_slice())).collect();
                poly_from_terms(&refs, 3)
            })
            .collect();
        check_system(&system, 3, 3, seed);
    }
}
