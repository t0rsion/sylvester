//! A self-contained Buchberger oracle for the differential suites.
//!
//! The oracle has its own polynomial arithmetic, grevlex from the
//! definition, textbook Buchberger completion, and interreduction. It shares
//! no code with `src/**`. The reduced Gröbner basis of an ideal is unique, so
//! each backend's output must equal the oracle's reduced basis as a set of
//! monic polynomials.
//!
//! Each test crate that includes this module uses a different subset, so the
//! whole module is exempt from `dead_code`.
#![allow(dead_code)]

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use sylvester::{Backend, ComputeError, ComputeOptions, GroebnerBasis, Ideal, PolynomialRing};

pub(crate) fn sub_mod(a: u64, b: u64, p: u64) -> u64 {
    (a + p - b % p) % p
}

pub(crate) fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

pub(crate) fn inv_mod(a: u64, p: u64) -> u64 {
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

pub(crate) type Exps = Vec<u16>;

pub(crate) fn grevlex_cmp(a: &[u16], b: &[u16]) -> Ordering {
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

pub(crate) fn divides(divisor: &[u16], multiple: &[u16]) -> bool {
    divisor.iter().zip(multiple.iter()).all(|(d, m)| d <= m)
}

pub(crate) fn exps_sub(a: &[u16], b: &[u16]) -> Exps {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}

pub(crate) fn exps_add(a: &[u16], b: &[u16]) -> Exps {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

pub(crate) fn exps_lcm(a: &[u16], b: &[u16]) -> Exps {
    a.iter().zip(b.iter()).map(|(x, y)| *x.max(y)).collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Poly {
    pub(crate) terms: BTreeMap<Exps, u64>,
}

impl Poly {
    pub(crate) fn zero() -> Self {
        Poly {
            terms: BTreeMap::new(),
        }
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub(crate) fn leading_term(&self) -> Option<(&Exps, u64)> {
        self.terms
            .iter()
            .max_by(|(a, _), (b, _)| grevlex_cmp(a, b))
            .map(|(m, &c)| (m, c))
    }

    /// Subtract `factor * x^shift * g` from `self`, modulo `p`.
    pub(crate) fn sub_scaled(&mut self, factor: u64, shift: &[u16], g: &Poly, p: u64) {
        for (mono, &coeff) in &g.terms {
            let target = exps_add(mono, shift);
            let delta = mul_mod(factor, coeff, p);
            let entry = self.terms.entry(target).or_insert(0);
            *entry = sub_mod(*entry, delta, p);
        }
        self.terms.retain(|_, c| *c != 0);
    }

    pub(crate) fn make_monic(&self, p: u64) -> Poly {
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

pub(crate) fn poly_from_terms(terms: &[(u64, &[u16])], p: u64) -> Poly {
    let mut poly = Poly::zero();
    for &(coeff, exps) in terms {
        let entry = poly.terms.entry(exps.to_vec()).or_insert(0);
        *entry = (*entry + coeff) % p;
    }
    poly.terms.retain(|_, c| *c != 0);
    poly
}

pub(crate) fn normal_form(f: &Poly, basis: &[Poly], p: u64) -> Poly {
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

pub(crate) fn s_polynomial(f: &Poly, g: &Poly, p: u64) -> Poly {
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
pub(crate) fn buchberger(generators: &[Poly], p: u64) -> Vec<Poly> {
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
pub(crate) fn reduce_basis(basis: &[Poly], p: u64) -> Vec<Poly> {
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

/// Deterministic splitmix64.
pub(crate) struct Rng(pub u64);

impl Rng {
    pub(crate) fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    pub(crate) fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

#[derive(Clone, Copy)]
pub(crate) struct Shape {
    pub(crate) nvars: usize,
    pub(crate) max_deg: u64,
    pub(crate) ngens: usize,
    pub(crate) min_terms: u64,
    pub(crate) term_spread: u64,
}

pub(crate) const fn shape(
    nvars: usize,
    max_deg: u64,
    ngens: usize,
    min_terms: u64,
    term_spread: u64,
) -> Shape {
    Shape {
        nvars,
        max_deg,
        ngens,
        min_terms,
        term_spread,
    }
}

/// 3 variables, 4 generators, degree <= 3, 1..=4 term-insertion attempts.
pub(crate) const SMALL: Shape = shape(3, 3, 4, 1, 4);

/// 4 variables, 5 generators, degree <= 3, 2..=5 term-insertion attempts.
pub(crate) const HARDER: Shape = shape(4, 3, 5, 2, 4);

/// A random monomial in `shape.nvars` variables of total degree <= max_deg.
pub(crate) fn random_mono(rng: &mut Rng, shape: Shape) -> Exps {
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
pub(crate) fn random_system(rng: &mut Rng, p: u64, shape: Shape) -> Vec<Poly> {
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
pub(crate) fn variable_names(nvars: usize) -> Vec<String> {
    (0..nvars).map(|index| format!("x{index}")).collect()
}

pub(crate) fn to_engine_ideal(system: &[Poly], p: u64, nvars: usize) -> Ideal {
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

pub(crate) fn engine_output_to_polys(output: &GroebnerBasis, p: u64, nvars: usize) -> Vec<Poly> {
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

pub(crate) fn canonical_set(basis: &[Poly], p: u64) -> BTreeSet<Vec<(Exps, u64)>> {
    basis
        .iter()
        .filter(|f| !f.is_zero())
        .map(|f| make_canonical(&f.make_monic(p)))
        .collect()
}

pub(crate) fn make_canonical(f: &Poly) -> Vec<(Exps, u64)> {
    f.terms.iter().map(|(e, &c)| (e.clone(), c)).collect()
}

pub(crate) fn render_system(system: &[Poly]) -> String {
    let rendered: Vec<String> = system.iter().map(|f| format!("{:?}", f.terms)).collect();
    format!("[{}]", rendered.join(", "))
}

fn check_backend(
    label: &str,
    system: &[Poly],
    oracle_reduced: &BTreeSet<Vec<(Exps, u64)>>,
    output: Result<GroebnerBasis, ComputeError>,
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

/// Runs both backends on one system and compares each against the oracle's
/// reduced basis.
pub(crate) fn check_system(system: &[Poly], p: u64, nvars: usize, seed: u64) {
    if system.iter().all(Poly::is_zero) {
        return;
    }
    let oracle = buchberger(system, p);
    let oracle_reduced = canonical_set(&reduce_basis(&oracle, p), p);
    let ideal = to_engine_ideal(system, p, nvars);

    for (label, backend) in [("classic", Backend::Classic), ("f4", Backend::F4)] {
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
