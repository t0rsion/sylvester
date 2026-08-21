//! Known-defect tests: verified counterexamples against both Gröbner backends.
//!
//! The checker below is independent of the engine code: it has
//! its own polynomial representation, its own mod-p arithmetic, a grevlex
//! comparison written from the definition, its own multivariate division,
//! and its own S-polynomials. It shares nothing with `src/**`, so it cannot
//! inherit an engine bug.
//!
//! Some tests validate the checker itself against hand-verified data. The
//! backend tests assert the correct behavior: they fail on the archived
//! snapshot, where they are `#[ignore]`d, and pass on the repaired backends
//! in this tree.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use sylvester::{Backend, ComputeOptions, GroebnerBasis, Ideal, PolynomialRing};

fn add_mod(a: u64, b: u64, p: u64) -> u64 {
    (a + b) % p
}

fn sub_mod(a: u64, b: u64, p: u64) -> u64 {
    (a + p - b % p) % p
}

fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

/// Modular inverse via Fermat's little theorem; `a` must be nonzero mod prime `p`.
fn inv_mod(a: u64, p: u64) -> u64 {
    assert!(a % p != 0, "attempted to invert zero mod {p}");
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

/// Graded reverse lexicographic order: first compare total degree; on a tie,
/// scan the exponents from the LAST variable towards the first, and at the
/// first position where they differ, the monomial with the SMALLER exponent
/// there is the LARGER monomial.
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

    fn from_terms(terms: &[(u64, &[u16])], p: u64) -> Self {
        let mut poly = Poly::zero();
        for (coeff, exps) in terms {
            let entry = poly.terms.entry(exps.to_vec()).or_insert(0);
            *entry = add_mod(*entry, coeff % p, p);
        }
        poly.terms.retain(|_, c| *c != 0);
        poly
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Leading term under grevlex, or None for the zero polynomial.
    fn leading_term(&self) -> Option<(&Exps, u64)> {
        self.terms
            .iter()
            .max_by(|(a, _), (b, _)| grevlex_cmp(a, b))
            .map(|(m, &c)| (m, c))
    }

    /// Subtract `factor * x^shift * g` from `self`, modulo `p`.
    fn sub_scaled(&mut self, factor: u64, shift: &[u16], g: &Poly, p: u64) {
        for (mono, &coeff) in &g.terms {
            let target = exps_add(mono, shift);
            let delta = mul_mod(factor, coeff, p);
            let entry = self.terms.entry(target).or_insert(0);
            *entry = sub_mod(*entry, delta, p);
        }
        self.terms.retain(|_, c| *c != 0);
    }

    /// Render for diagnostics, e.g. "x^2 + 2*y*z".
    fn render(&self, p: u64) -> String {
        if self.is_zero() {
            return "0".to_string();
        }
        let mut monomials: Vec<&Exps> = self.terms.keys().collect();
        monomials.sort_by(|a, b| grevlex_cmp(b, a));
        let var_name = |i: usize, nvars: usize| -> String {
            if nvars <= 3 {
                ["x", "y", "z"][i].to_string()
            } else {
                format!("x{i}")
            }
        };
        let pieces: Vec<String> = monomials
            .iter()
            .map(|mono| {
                let coeff = self.terms[*mono] % p;
                let mut factors: Vec<String> = Vec::new();
                if coeff != 1 || mono.iter().all(|&e| e == 0) {
                    factors.push(coeff.to_string());
                }
                for (i, &e) in mono.iter().enumerate() {
                    match e {
                        0 => {}
                        1 => factors.push(var_name(i, mono.len())),
                        _ => factors.push(format!("{}^{}", var_name(i, mono.len()), e)),
                    }
                }
                factors.join("*")
            })
            .collect();
        pieces.join(" + ")
    }
}

/// Normal form of `f` on division by `basis` (zero elements of `basis` are
/// skipped). Standard multivariate division: reduce the grevlex-leading term
/// while some basis leading monomial divides it, otherwise move it to the
/// remainder.
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

/// S-polynomial of `f` and `g` (both nonzero).
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

#[derive(Clone, Debug, PartialEq, Eq)]
enum Witness {
    /// S(candidate[i], candidate[j]) has this nonzero normal form.
    SPairDoesNotReduce {
        i: usize,
        j: usize,
        normal_form: Poly,
    },
    /// ideal_gens[generator] has this nonzero normal form modulo the candidate.
    GeneratorDoesNotReduce { generator: usize, normal_form: Poly },
}

/// Checks that `candidate` is a Gröbner basis whose ideal contains the ideal
/// generated by `ideal_gens`: (a) every S-pair of `candidate` reduces to zero
/// modulo `candidate` (Buchberger's criterion), and (b) every element of
/// `ideal_gens` reduces to zero modulo `candidate`.
///
/// This does not check that `candidate` lies in the ideal `ideal_gens`
/// generates, so it is necessary and not sufficient for `candidate` to be a
/// Gröbner basis of the input ideal. A candidate of {1} passes. The
/// `sylv-gb-cert-v1` certificate closes that gap by checking ideal equality
/// with the input; this checker does not.
fn is_groebner_basis_of(ideal_gens: &[Poly], candidate: &[Poly], p: u64) -> Result<(), Witness> {
    let nonzero: Vec<Poly> = candidate.iter().filter(|f| !f.is_zero()).cloned().collect();
    for i in 0..nonzero.len() {
        for j in (i + 1)..nonzero.len() {
            let s = s_polynomial(&nonzero[i], &nonzero[j], p);
            let r = normal_form(&s, &nonzero, p);
            if !r.is_zero() {
                return Err(Witness::SPairDoesNotReduce {
                    i,
                    j,
                    normal_form: r,
                });
            }
        }
    }
    for (generator, f) in ideal_gens.iter().enumerate() {
        let r = normal_form(f, &nonzero, p);
        if !r.is_zero() {
            return Err(Witness::GeneratorDoesNotReduce {
                generator,
                normal_form: r,
            });
        }
    }
    Ok(())
}

/// Counterexample 1 over F_2[x, y], grevlex: f1 = x^3 + y^3,
/// f2 = x + x^2 + x^2*y, f3 = y + x^2.
/// True reduced Gröbner basis: {x, y}.
fn f2_system() -> Vec<Poly> {
    vec![
        Poly::from_terms(&[(1, &[3, 0]), (1, &[0, 3])], 2),
        Poly::from_terms(&[(1, &[1, 0]), (1, &[2, 0]), (1, &[2, 1])], 2),
        Poly::from_terms(&[(1, &[0, 1]), (1, &[2, 0])], 2),
    ]
}

fn f2_system_ideal() -> Ideal {
    engine_ideal(2, &["x^3 + y^3", "x + x^2 + x^2*y", "y + x^2"])
}

fn f2_true_basis() -> Vec<Poly> {
    vec![
        Poly::from_terms(&[(1, &[1, 0])], 2), // x
        Poly::from_terms(&[(1, &[0, 1])], 2), // y
    ]
}

/// The classic backend's actual (wrong) output on the F_2 system: {y^2, x + y}.
fn f2_wrong_engine_output() -> Vec<Poly> {
    vec![
        Poly::from_terms(&[(1, &[0, 2])], 2),               // y^2
        Poly::from_terms(&[(1, &[1, 0]), (1, &[0, 1])], 2), // x + y
    ]
}

/// Counterexample 2 over F_3[x, y, z], grevlex: f1 = x^2 + y^2, f2 = xz + xy,
/// f3 = y + xy.
/// True reduced Gröbner basis: {yz^2 + y, x^2 + 2yz, xy + y, y^2 + yz, xz + 2y}.
fn f3_system() -> Vec<Poly> {
    vec![
        Poly::from_terms(&[(1, &[2, 0, 0]), (1, &[0, 2, 0])], 3),
        Poly::from_terms(&[(1, &[1, 0, 1]), (1, &[1, 1, 0])], 3),
        Poly::from_terms(&[(1, &[0, 1, 0]), (1, &[1, 1, 0])], 3),
    ]
}

fn f3_system_ideal() -> Ideal {
    engine_ideal(3, &["x^2 + y^2", "x*z + x*y", "y + x*y"])
}

/// The matrix backend's actual (wrong) 4-element output on the F_3 system.
fn f3_wrong_engine_output() -> Vec<Poly> {
    vec![
        Poly::from_terms(&[(1, &[2, 0, 0]), (2, &[0, 1, 1])], 3), // x^2 + 2yz
        Poly::from_terms(&[(1, &[1, 1, 0]), (1, &[0, 1, 0])], 3), // xy + y
        Poly::from_terms(&[(1, &[0, 2, 0]), (1, &[0, 1, 1])], 3), // y^2 + yz
        Poly::from_terms(&[(1, &[1, 0, 1]), (2, &[0, 1, 0])], 3), // xz + 2y
    ]
}

fn f3_true_basis() -> Vec<Poly> {
    let mut basis = vec![Poly::from_terms(&[(1, &[0, 1, 2]), (1, &[0, 1, 0])], 3)]; // yz^2 + y
    basis.extend(f3_wrong_engine_output());
    basis
}

/// Build the ideal of a system over `F_p`, from the engine's own syntax.
fn engine_ideal(p: u64, generators: &[&str]) -> Ideal {
    let names = ["x", "y", "z"];
    let nvars = if p == 2 { 2 } else { 3 };
    let ring = PolynomialRing::prime_field(p, names[..nvars].iter().copied())
        .expect("the modulus is prime");
    let polys: Vec<_> = generators
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the text parses"))
        .collect();
    ring.ideal(polys).expect("the polynomials share a ring")
}

fn engine_output_to_polys(output: &GroebnerBasis, nvars: usize, p: u64) -> Vec<Poly> {
    output
        .iter()
        .map(|poly| {
            let terms: Vec<(u64, Exps)> = poly
                .terms()
                .map(|(coeff, exps)| {
                    assert!(exps.len() <= nvars, "engine output has too many variables");
                    let mut exps = exps.to_vec();
                    exps.resize(nvars, 0);
                    (coeff % p, exps)
                })
                .collect();
            let borrowed: Vec<(u64, &[u16])> =
                terms.iter().map(|(c, e)| (*c, e.as_slice())).collect();
            Poly::from_terms(&borrowed, p)
        })
        .collect()
}

fn render_basis(basis: &[Poly], p: u64) -> String {
    let rendered: Vec<String> = basis.iter().map(|f| f.render(p)).collect();
    format!("{{{}}}", rendered.join(", "))
}

fn render_witness(witness: &Witness, p: u64) -> String {
    match witness {
        Witness::SPairDoesNotReduce { i, j, normal_form } => format!(
            "S-pair ({i}, {j}) has nonzero normal form {}",
            normal_form.render(p)
        ),
        Witness::GeneratorDoesNotReduce {
            generator,
            normal_form,
        } => format!(
            "input generator {generator} has nonzero normal form {}",
            normal_form.render(p)
        ),
    }
}

fn make_monic(f: &Poly, p: u64) -> Poly {
    let Some((_, lead)) = f.leading_term() else {
        return Poly::zero();
    };
    let scale = inv_mod(lead, p);
    let mut out = Poly::zero();
    for (exps, &coeff) in &f.terms {
        out.terms.insert(exps.clone(), mul_mod(coeff, scale, p));
    }
    out
}

/// Canonical form of a basis for exact comparison: nonzero elements, made
/// monic, rendered (the rendering is canonical because term maps are ordered).
fn normalized_basis_set(basis: &[Poly], p: u64) -> BTreeSet<String> {
    basis
        .iter()
        .filter(|f| !f.is_zero())
        .map(|f| make_monic(f, p).render(p))
        .collect()
}

fn assert_engine_output_is_basis(
    label: &str,
    gens: &[Poly],
    true_basis: &[Poly],
    output: &GroebnerBasis,
    nvars: usize,
    p: u64,
) {
    let candidate = engine_output_to_polys(output, nvars, p);
    if let Err(witness) = is_groebner_basis_of(gens, &candidate, p) {
        panic!(
            "{label}: engine returned {}, which fails the independent checker: {}",
            render_basis(&candidate, p),
            render_witness(&witness, p)
        );
    }
    // The containment checker alone cannot reject a proper superset ideal such
    // as {1}; for these fixed systems the reduced basis is known, so demand it
    // exactly (as a set of monic polynomials).
    assert_eq!(
        normalized_basis_set(&candidate, p),
        normalized_basis_set(true_basis, p),
        "{label}: engine returned {} instead of the known reduced basis {}",
        render_basis(&candidate, p),
        render_basis(true_basis, p)
    );
}

#[test]
fn checker_accepts_true_basis_for_f2_system() {
    assert_eq!(
        is_groebner_basis_of(&f2_system(), &f2_true_basis(), 2),
        Ok(())
    );
}

#[test]
fn checker_accepts_true_basis_for_f3_system() {
    assert_eq!(
        is_groebner_basis_of(&f3_system(), &f3_true_basis(), 3),
        Ok(())
    );
}

#[test]
fn checker_rejects_classic_f2_output_with_f2_generator_witness() {
    // {y^2, x + y} is a reduced Gröbner basis, but of the wrong ideal: every
    // S-pair reduces to zero, but the input generator f2 = x + x^2 + x^2*y
    // has normal form y, so the input ideal is not contained.
    let expected_normal_form = Poly::from_terms(&[(1, &[0, 1])], 2); // y
    let result = is_groebner_basis_of(&f2_system(), &f2_wrong_engine_output(), 2);
    assert_eq!(
        result,
        Err(Witness::GeneratorDoesNotReduce {
            generator: 1,
            normal_form: expected_normal_form,
        })
    );
}

#[test]
fn checker_rejects_matrix_f3_output_with_s_pair_witness() {
    // The 4-element matrix output is not a Gröbner basis at all:
    // S(x^2 + 2yz, xy + y) has nonzero normal form y + yz^2.
    let expected_normal_form = Poly::from_terms(&[(1, &[0, 1, 0]), (1, &[0, 1, 2])], 3); // y + yz^2
    let result = is_groebner_basis_of(&f3_system(), &f3_wrong_engine_output(), 3);
    match result {
        Err(Witness::SPairDoesNotReduce {
            i: 0,
            j: 1,
            normal_form,
        }) => assert_eq!(normal_form, expected_normal_form),
        other => panic!("expected the S(g0, g1) witness with normal form y + yz^2, got {other:?}"),
    }
}

#[test]
fn classic_backend_solves_the_f2_system() {
    let output = f2_system_ideal()
        .groebner_basis(ComputeOptions::new().backend(Backend::Classic))
        .expect("classic backend returned an error");
    assert_engine_output_is_basis(
        "classic backend, F_2 system",
        &f2_system(),
        &f2_true_basis(),
        &output,
        2,
        2,
    );
}

#[test]
fn matrix_backend_solves_the_f3_system() {
    let output = f3_system_ideal()
        .groebner_basis(ComputeOptions::new().backend(Backend::Matrix))
        .expect("matrix backend returned an error");
    assert_engine_output_is_basis(
        "matrix backend, F_3 system",
        &f3_system(),
        &f3_true_basis(),
        &output,
        3,
        3,
    );
}

#[test]
fn the_certified_path_solves_both_systems() {
    let certified = f2_system_ideal()
        .groebner_basis_certified(ComputeOptions::new())
        .expect("certification of the F_2 system returned an error");
    assert_engine_output_is_basis(
        "certified classic backend, F_2 system",
        &f2_system(),
        &f2_true_basis(),
        certified.basis(),
        2,
        2,
    );

    let certified = f3_system_ideal()
        .groebner_basis_certified(ComputeOptions::new())
        .expect("certification of the F_3 system returned an error");
    assert_engine_output_is_basis(
        "certified classic backend, F_3 system",
        &f3_system(),
        &f3_true_basis(),
        certified.basis(),
        3,
        3,
    );
}

#[test]
fn the_default_options_solve_both_systems() {
    let output = f2_system_ideal()
        .groebner_basis(ComputeOptions::new())
        .expect("the default options returned an error on the F_2 system");
    assert_engine_output_is_basis(
        "default options, F_2 system",
        &f2_system(),
        &f2_true_basis(),
        &output,
        2,
        2,
    );

    let output = f3_system_ideal()
        .groebner_basis(ComputeOptions::new())
        .expect("the default options returned an error on the F_3 system");
    assert_engine_output_is_basis(
        "default options, F_3 system",
        &f3_system(),
        &f3_true_basis(),
        &output,
        3,
        3,
    );
}

#[test]
fn exact_output_assertion_rejects_the_trivial_basis() {
    // {1} generates the whole ring, so it passes the containment checker;
    // the exact-set comparison is what rejects it.
    let one = Poly::from_terms(&[(1, &[0, 0])], 2);
    assert_eq!(
        is_groebner_basis_of(&f2_system(), std::slice::from_ref(&one), 2),
        Ok(())
    );
    assert_ne!(
        normalized_basis_set(std::slice::from_ref(&one), 2),
        normalized_basis_set(&f2_true_basis(), 2)
    );
}
