//! End-to-end tests for `Ideal::groebner_basis_certified`.
//!
//! Every round trip runs the certified entry point, hands the bytes back to
//! the independent verifier, and compares the basis against the raw
//! classic-backend output on the same input. The certificate holds by
//! construction, so a failure here is an engine defect, not a disagreement
//! about the contract.

use std::collections::{BTreeMap, BTreeSet};
use std::time::{Duration, Instant};

use sylvester::verify::{Poly, verify};
use sylvester::{
    Backend, CertifyError, ComputeError, ComputeOptions, GroebnerBasis, Ideal, Polynomial,
    PolynomialRing, RingError,
};

fn classic() -> ComputeOptions {
    ComputeOptions::new().backend(Backend::Classic)
}

fn ring(modulus: u64, nvars: usize) -> PolynomialRing {
    let names: Vec<String> = (0..nvars).map(|index| format!("x{index}")).collect();
    PolynomialRing::prime_field(modulus, names).expect("the modulus is prime")
}

/// Build an ideal from terms given as (coefficient, exponents).
fn ideal_of(ring: &PolynomialRing, system: &[Vec<(i64, Vec<u16>)>]) -> Ideal {
    let generators: Vec<Polynomial> = system
        .iter()
        .map(|terms| {
            ring.polynomial(terms.iter().cloned())
                .expect("the exponents fit the ring")
        })
        .collect();
    ring.ideal(generators)
        .expect("the polynomials share a ring")
}

/// A basis as a set of monic term maps, for comparison without an order.
fn canonical_set(basis: &[Polynomial], p: u64) -> BTreeSet<Vec<(Vec<u16>, u64)>> {
    basis
        .iter()
        .map(|poly| {
            let lead = poly
                .leading_term()
                .map(|(coeff, _)| coeff.value())
                .unwrap_or(1);
            let scale = mod_inverse(lead, p);
            poly.terms()
                .map(|(coeff, exps)| (exps.to_vec(), mul_mod(coeff.value(), scale, p)))
                .collect::<BTreeMap<Vec<u16>, u64>>()
                .into_iter()
                .collect()
        })
        .collect()
}

/// The certificate's own basis, read back as ring polynomials.
fn certificate_basis(ring: &PolynomialRing, verified: &[Poly]) -> Vec<Polynomial> {
    verified
        .iter()
        .map(|poly| {
            let terms: Vec<(i64, Vec<u16>)> = poly
                .terms()
                .iter()
                .map(|term| {
                    let exps = term.mono().exps().iter().map(|&e| e as u16).collect();
                    (term.coeff() as i64, exps)
                })
                .collect();
            ring.polynomial(terms).expect("the exponents fit the ring")
        })
        .collect()
}

fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn mod_inverse(a: u64, p: u64) -> u64 {
    let mut acc = 1u64;
    let mut base = a % p;
    let mut exp = p - 2;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base, p);
        }
        base = mul_mod(base, base, p);
        exp >>= 1;
    }
    acc
}

/// Certify one system and check every claim the API makes about it.
fn round_trip(label: &str, ideal: &Ideal) -> Vec<u8> {
    let ring = ideal.ring();
    let certified = ideal
        .groebner_basis_certified(classic())
        .unwrap_or_else(|error| panic!("{label}: certification failed: {error}"));

    let verified = verify(certified.certificate())
        .unwrap_or_else(|error| panic!("{label}: the verifier rejected the certificate: {error}"));
    assert_eq!(verified.modulus(), ring.modulus(), "{label}: modulus");
    assert_eq!(verified.nvars(), ring.nvars(), "{label}: nvars");

    let raw = ideal
        .groebner_basis(classic())
        .unwrap_or_else(|error| panic!("{label}: the raw entry point failed: {error}"));
    assert_eq!(
        certified.basis().to_vec(),
        raw.to_vec(),
        "{label}: the certified basis must match the raw basis"
    );
    assert_eq!(
        certificate_basis(ring, verified.basis()),
        certified.basis().to_vec(),
        "{label}: the certificate must carry the basis the API returns"
    );

    certified.into_parts().1
}

/// F_2[x0, x1]: f1 = x0^3 + x1^3, f2 = x0 + x0^2 + x0^2*x1, f3 = x1 + x0^2.
fn f2_system() -> Vec<Vec<(i64, Vec<u16>)>> {
    vec![
        vec![(1, vec![3, 0]), (1, vec![0, 3])],
        vec![(1, vec![1, 0]), (1, vec![2, 0]), (1, vec![2, 1])],
        vec![(1, vec![0, 1]), (1, vec![2, 0])],
    ]
}

/// The hand-verified reduced basis of the F_2 system: {x0, x1}.
fn f2_true_basis() -> Vec<Vec<(i64, Vec<u16>)>> {
    vec![vec![(1, vec![1, 0])], vec![(1, vec![0, 1])]]
}

/// F_3[x0, x1, x2]: f1 = x0^2 + x1^2, f2 = x0*x2 + x0*x1, f3 = x1 + x0*x1.
fn f3_system() -> Vec<Vec<(i64, Vec<u16>)>> {
    vec![
        vec![(1, vec![2, 0, 0]), (1, vec![0, 2, 0])],
        vec![(1, vec![1, 0, 1]), (1, vec![1, 1, 0])],
        vec![(1, vec![0, 1, 0]), (1, vec![1, 1, 0])],
    ]
}

/// The hand-verified reduced basis of the F_3 system.
fn f3_true_basis() -> Vec<Vec<(i64, Vec<u16>)>> {
    vec![
        vec![(1, vec![0, 1, 2]), (1, vec![0, 1, 0])],
        vec![(1, vec![2, 0, 0]), (2, vec![0, 1, 1])],
        vec![(1, vec![1, 1, 0]), (1, vec![0, 1, 0])],
        vec![(1, vec![0, 2, 0]), (1, vec![0, 1, 1])],
        vec![(1, vec![1, 0, 1]), (2, vec![0, 1, 0])],
    ]
}

/// Cyclic-4 over F_32003.
fn cyclic_4() -> Vec<Vec<(i64, Vec<u16>)>> {
    vec![
        vec![
            (1, vec![1, 0, 0, 0]),
            (1, vec![0, 1, 0, 0]),
            (1, vec![0, 0, 1, 0]),
            (1, vec![0, 0, 0, 1]),
        ],
        vec![
            (1, vec![1, 1, 0, 0]),
            (1, vec![0, 1, 1, 0]),
            (1, vec![0, 0, 1, 1]),
            (1, vec![1, 0, 0, 1]),
        ],
        vec![
            (1, vec![1, 1, 1, 0]),
            (1, vec![0, 1, 1, 1]),
            (1, vec![1, 0, 1, 1]),
            (1, vec![1, 1, 0, 1]),
        ],
        vec![(1, vec![1, 1, 1, 1]), (-1, vec![0, 0, 0, 0])],
    ]
}

/// Cyclic-5 over F_32003, from the text the ring parses.
fn cyclic_5(ring: &PolynomialRing) -> Ideal {
    let texts = [
        "x0 + x1 + x2 + x3 + x4",
        "x0*x1 + x1*x2 + x2*x3 + x3*x4 + x4*x0",
        "x0*x1*x2 + x1*x2*x3 + x2*x3*x4 + x3*x4*x0 + x4*x0*x1",
        "x0*x1*x2*x3 + x1*x2*x3*x4 + x2*x3*x4*x0 + x3*x4*x0*x1 + x4*x0*x1*x2",
        "x0*x1*x2*x3*x4 - 1",
    ];
    let generators: Vec<Polynomial> = texts
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the text parses"))
        .collect();
    ring.ideal(generators)
        .expect("the polynomials share a ring")
}

#[test]
fn certifies_the_f2_counterexample_system() {
    let ring = ring(2, 2);
    let bytes = round_trip("F_2 counterexample", &ideal_of(&ring, &f2_system()));
    let verified = verify(&bytes).expect("the certificate holds");
    assert_eq!(
        canonical_set(&certificate_basis(&ring, verified.basis()), 2),
        canonical_set(ideal_of(&ring, &f2_true_basis()).generators(), 2),
        "the certified basis must be the hand-verified reduced basis"
    );
}

#[test]
fn certifies_the_f3_counterexample_system() {
    let ring = ring(3, 3);
    let bytes = round_trip("F_3 counterexample", &ideal_of(&ring, &f3_system()));
    let verified = verify(&bytes).expect("the certificate holds");
    assert_eq!(
        canonical_set(&certificate_basis(&ring, verified.basis()), 3),
        canonical_set(ideal_of(&ring, &f3_true_basis()).generators(), 3),
        "the certified basis must be the hand-verified reduced basis"
    );
}

#[test]
fn certifies_cyclic_4() {
    round_trip("cyclic-4", &ideal_of(&ring(32003, 4), &cyclic_4()));
}

#[test]
fn certifies_a_single_generator() {
    let ring = ring(32003, 2);
    let system = vec![vec![(1i64, vec![2u16, 0]), (3, vec![0, 1])]];
    let bytes = round_trip("single generator", &ideal_of(&ring, &system));
    let verified = verify(&bytes).expect("the certificate holds");
    assert_eq!(verified.basis().len(), 1);
    assert_eq!(verified.input().len(), 1);
}

#[test]
fn certifies_the_unit_ideal() {
    let ring = ring(7, 2);
    let system = vec![
        vec![(1i64, vec![1u16, 0])],
        vec![(1, vec![1, 0]), (1, vec![0, 0])],
    ];
    let bytes = round_trip("unit ideal", &ideal_of(&ring, &system));
    let verified = verify(&bytes).expect("the certificate holds");
    assert_eq!(verified.basis().len(), 1);
    assert_eq!(verified.basis()[0].terms().len(), 1);
    assert_eq!(verified.basis()[0].terms()[0].coeff(), 1);
    assert_eq!(verified.basis()[0].terms()[0].mono().degree(), 0);
}

#[test]
fn certifies_the_zero_ideal() {
    // The coefficient 7 reduces to zero, so the generator vanishes and the
    // ideal is zero. The ring still fixes the variable count.
    let ring = ring(7, 2);
    let system = vec![vec![(7i64, vec![0u16, 0])]];
    let bytes = round_trip("zero ideal", &ideal_of(&ring, &system));
    let verified = verify(&bytes).expect("the certificate holds");
    assert!(verified.basis().is_empty());
    assert_eq!(verified.nvars(), 2);

    // A ring with no variables names none in the certificate either.
    let field = PolynomialRing::prime_field(7, Vec::<String>::new()).expect("7 is prime");
    let bytes = round_trip("no variables", &ideal_of(&field, &[Vec::new()]));
    let verified = verify(&bytes).expect("the certificate holds");
    assert!(verified.basis().is_empty());
    assert_eq!(verified.nvars(), 0);
}

#[test]
fn certifies_an_input_with_no_generators() {
    let bytes = round_trip("no generators", &ideal_of(&ring(7, 2), &[]));
    let verified = verify(&bytes).expect("the certificate holds");
    assert!(verified.basis().is_empty());
    assert!(verified.input().is_empty());
}

#[test]
fn the_certified_basis_derefs_to_a_slice() {
    let ideal = ideal_of(&ring(32003, 4), &cyclic_4());
    let certified = ideal
        .groebner_basis_certified(classic())
        .expect("the certificate holds");
    let basis: &GroebnerBasis = certified.basis();
    assert!(!basis.is_empty());
    assert_eq!(basis.iter().count(), basis.len());
    assert!(basis.iter().all(|poly| !poly.is_zero()));
}

/// Splitmix64, matching `tests/differential.rs`.
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

/// A random monomial in `nvars` variables of total degree at most 3.
fn random_mono(rng: &mut Rng, nvars: usize) -> Vec<u16> {
    loop {
        let exps: Vec<u16> = (0..nvars).map(|_| rng.below(4) as u16).collect();
        if exps.iter().map(|&e| e as u64).sum::<u64>() <= 3 {
            return exps;
        }
    }
}

/// Four generators over three variables; term counts are insertion
/// attempts, so colliding monomials can leave fewer terms.
fn random_system(rng: &mut Rng, p: u64) -> Vec<Vec<(i64, Vec<u16>)>> {
    (0..4)
        .map(|_| {
            let nterms = 1 + rng.below(4);
            let mut terms: BTreeMap<Vec<u16>, u64> = BTreeMap::new();
            for _ in 0..nterms {
                let mono = random_mono(rng, 3);
                let coeff = 1 + rng.below(p - 1);
                let entry = terms.entry(mono).or_insert(0);
                *entry = (*entry + coeff) % p;
            }
            terms
                .into_iter()
                .filter(|(_, coeff)| *coeff != 0)
                .map(|(exps, coeff)| (coeff as i64, exps))
                .collect()
        })
        .collect()
}

fn certify_random_systems(p: u64, seeds: std::ops::Range<u64>) {
    let ring = ring(p, 3);
    for seed in seeds {
        let mut rng = Rng(seed);
        let system = random_system(&mut rng, p);
        round_trip(
            &format!("random p={p} seed={seed}"),
            &ideal_of(&ring, &system),
        );
    }
}

#[test]
fn certifies_random_systems_over_f2() {
    certify_random_systems(2, 0..12);
}

#[test]
fn certifies_random_systems_over_f3() {
    certify_random_systems(3, 0..12);
}

#[test]
fn certifies_random_systems_over_f5() {
    certify_random_systems(5, 0..12);
}

#[test]
fn two_runs_give_the_same_certificate_bytes() {
    for (label, ideal) in [
        ("F_2 counterexample", ideal_of(&ring(2, 2), &f2_system())),
        ("F_3 counterexample", ideal_of(&ring(3, 3), &f3_system())),
        ("cyclic-4", ideal_of(&ring(32003, 4), &cyclic_4())),
    ] {
        let first = ideal
            .groebner_basis_certified(classic())
            .unwrap_or_else(|error| panic!("{label}: {error}"));
        let second = ideal
            .groebner_basis_certified(classic())
            .unwrap_or_else(|error| panic!("{label}: {error}"));
        assert_eq!(
            first.certificate(),
            second.certificate(),
            "{label}: the bytes must be a function of the input"
        );
        assert_eq!(first.basis(), second.basis(), "{label}");
    }
}

#[test]
fn a_timeout_too_large_to_represent_does_not_panic() {
    // `Instant::now() + Duration::MAX` overflows. The run must still
    // complete, with no deadline rather than a panic.
    let ideal = ideal_of(&ring(32003, 4), &cyclic_4());
    let certified = ideal
        .groebner_basis_certified(classic().timeout(Duration::MAX))
        .expect("a timeout too large to add still lets the run finish");
    assert!(!certified.basis().is_empty());
}

#[test]
fn an_exhausted_timeout_stops_certification() {
    let ideal = ideal_of(&ring(32003, 4), &cyclic_4());
    assert_eq!(
        ideal.groebner_basis_certified(classic().timeout(Duration::ZERO)),
        Err(CertifyError::Engine(ComputeError::Timeout))
    );
}

#[test]
fn a_pair_past_the_degree_limit_stops_certification() {
    let ring = ring(7, 2);
    let system = vec![
        vec![(1i64, vec![65535u16, 0]), (1, vec![0, 65535])],
        vec![(1, vec![0, 65535])],
    ];
    let ideal = ideal_of(&ring, &system);
    assert_eq!(
        ideal.groebner_basis_certified(classic()),
        Err(CertifyError::Engine(ComputeError::DegreeLimit {
            limit: 65535
        }))
    );
    assert_eq!(
        ideal.groebner_basis(classic()),
        Err(ComputeError::DegreeLimit { limit: 65535 })
    );
}

#[test]
fn an_exhausted_timeout_stops_a_run_with_no_pairs() {
    // A single generator makes no critical pair, and no generator makes no
    // basis. Neither run may report success past an exhausted budget.
    let ring = ring(32003, 2);
    let single = vec![vec![(1i64, vec![2u16, 0]), (3, vec![0, 1])]];
    assert_eq!(
        ideal_of(&ring, &single).groebner_basis(classic().timeout(Duration::ZERO)),
        Err(ComputeError::Timeout)
    );
    assert_eq!(
        ideal_of(&ring, &[]).groebner_basis(classic().timeout(Duration::ZERO)),
        Err(ComputeError::Timeout)
    );
}

#[test]
fn the_certificate_format_follows_the_backend() {
    // The classic backend writes `sylv-gb-cert-v1`, which is JSON, and the
    // F4 backend writes `sylv-gb-cert-v2`, which starts with the magic.
    let ideal = ideal_of(&ring(3, 3), &f3_system());
    let classic_run = ideal
        .groebner_basis_certified(classic())
        .expect("the certificate holds");
    assert_eq!(classic_run.certificate()[0], b'{');
    let f4_run = ideal
        .groebner_basis_certified(ComputeOptions::new())
        .expect("the certificate holds");
    assert_eq!(&f4_run.certificate()[..8], b"SYLVGB\x02\x00");
    assert_eq!(
        canonical_set(classic_run.basis(), 3),
        canonical_set(f4_run.basis(), 3)
    );
}

#[test]
fn a_non_prime_modulus_has_no_ring() {
    assert_eq!(
        PolynomialRing::prime_field(4, ["x", "y"]),
        Err(RingError::ModulusNotPrime { modulus: 4 })
    );
}

#[test]
fn the_uncertified_path_still_solves_the_counterexample_systems() {
    let f2_ring = ring(2, 2);
    let output = ideal_of(&f2_ring, &f2_system())
        .groebner_basis(classic())
        .expect("the classic backend runs");
    assert_eq!(
        canonical_set(&output, 2),
        canonical_set(ideal_of(&f2_ring, &f2_true_basis()).generators(), 2),
        "the raw classic output on the F_2 system must stay the known reduced basis"
    );

    let f3_ring = ring(3, 3);
    let output = ideal_of(&f3_ring, &f3_system())
        .groebner_basis(classic())
        .expect("the classic backend runs");
    assert_eq!(
        canonical_set(&output, 3),
        canonical_set(ideal_of(&f3_ring, &f3_true_basis()).generators(), 3),
        "the raw classic output on the F_3 system must stay the known reduced basis"
    );
}

#[test]
fn the_certified_basis_is_the_basis_the_certificate_carries() {
    // The value the caller reads comes from the accepted bytes. Decoding
    // the certificate with the verifier must give back the same basis, in
    // the same order.
    let ring = ring(32003, 4);
    let certified = ideal_of(&ring, &cyclic_4())
        .groebner_basis_certified(classic())
        .expect("the certificate holds");
    let verified = verify(certified.certificate()).expect("the certificate holds");

    assert_eq!(
        certificate_basis(&ring, verified.basis()),
        certified.basis().to_vec(),
        "the basis must be the one decoded from the certificate"
    );
    assert_eq!(certified.basis().ring(), &ring);
    let (basis, bytes) = certified.into_parts();
    assert_eq!(
        certificate_basis(&ring, verified.basis()),
        basis.into_polynomials(),
        "into_parts must hand out the same basis"
    );
    verify(&bytes).expect("the bytes the caller owns still hold");
}

#[test]
fn a_memory_limit_changes_no_certificate_byte() {
    for (label, ideal) in [
        ("F_3 counterexample", ideal_of(&ring(3, 3), &f3_system())),
        ("cyclic-4", ideal_of(&ring(32003, 4), &cyclic_4())),
    ] {
        let free = ideal
            .groebner_basis_certified(classic())
            .unwrap_or_else(|error| panic!("{label}: {error}"));
        let bounded = ideal
            .groebner_basis_certified(classic().memory_limit(64 << 20))
            .unwrap_or_else(|error| panic!("{label}: {error}"));
        assert_eq!(
            free.certificate(),
            bounded.certificate(),
            "{label}: a memory limit that holds the run changes no byte"
        );
        assert_eq!(free.basis(), bounded.basis(), "{label}");
    }
}

#[test]
fn a_tight_timeout_stops_the_certified_path_early() {
    // The certified path costs many times the raw run. It tracks the
    // cofactors, divides every input and every S-polynomial by the basis,
    // and verifies the bytes it wrote. A deadline well under the full cost
    // must stop it, and the value must report the budget, not the basis.
    let ring = ring(32003, 5);
    let ideal = cyclic_5(&ring);
    let start = Instant::now();
    ideal
        .groebner_basis_certified(classic())
        .expect("the certificate holds");
    let full = start.elapsed();

    let start = Instant::now();
    let error = ideal
        .groebner_basis_certified(classic().timeout(full / 20))
        .expect_err("the deadline stops the run");
    let elapsed = start.elapsed();

    assert_eq!(error, CertifyError::Engine(ComputeError::Timeout));
    assert!(
        elapsed < full / 2,
        "the run must stop at the deadline: {elapsed:?} of a full {full:?}"
    );
}

#[test]
fn a_tight_memory_limit_stops_the_certified_path() {
    // The raw engine fits this limit. The tracked cofactors, the
    // representations over the basis, and the certificate do not.
    const LIMIT: usize = 64 << 10;
    let ring = ring(32003, 5);
    let ideal = cyclic_5(&ring);
    ideal
        .groebner_basis(classic().memory_limit(LIMIT))
        .expect("the raw run fits the limit");

    assert_eq!(
        ideal.groebner_basis_certified(classic().memory_limit(LIMIT)),
        Err(CertifyError::Engine(ComputeError::MemoryLimitExceeded))
    );
}

#[test]
fn a_zero_memory_limit_stops_certification_of_the_empty_ideal() {
    // The engine holds nothing for an empty ideal, but the certificate
    // bytes are not empty: the schema, the order, and the empty arrays
    // still cost bytes. A limit of zero must not let those bytes through
    // uncounted.
    let ideal = ideal_of(&ring(7, 2), &[]);
    assert_eq!(
        ideal.groebner_basis_certified(classic().memory_limit(0)),
        Err(CertifyError::WriterExhausted(
            ComputeError::MemoryLimitExceeded
        ))
    );
}

#[test]
fn no_memory_limit_turns_certification_into_a_rejection() {
    // Every limit either holds the whole run or stops it with a typed
    // exhaustion. A rejection would say the certificate does not hold, and
    // no budget can prove that.
    let ideal = ideal_of(&ring(32003, 4), &cyclic_4());
    let mut held = 0;
    for step in 1..32 {
        let limit = step * 2048;
        match ideal.groebner_basis_certified(classic().memory_limit(limit)) {
            Ok(_) => held += 1,
            Err(
                CertifyError::Engine(_)
                | CertifyError::WriterExhausted(_)
                | CertifyError::VerifierExhausted(_),
            ) => {}
            Err(other) => panic!("a limit of {limit} bytes reported {other}"),
        }
    }
    assert!(held > 0, "a generous limit must hold the whole run");
}

#[test]
fn the_emitter_shares_no_code_with_the_verifier() {
    // `src/cert` is the untrusted emitter: it writes certificate bytes and
    // never reads them back. `groebner_basis_certified` is the one place
    // that hands those bytes to the independent verifier. `tests.rs`
    // checks its own output against `verify` as an end-to-end sanity
    // check, which does not weaken that boundary in the shipped code, so
    // it is the one file this scan skips.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/cert");
    let mut checked = 0;
    let mut pending = vec![root];
    while let Some(directory) = pending.pop() {
        for entry in std::fs::read_dir(&directory).expect("the emitter directory is readable") {
            let entry = entry.expect("the directory entry is readable").path();
            if entry.is_dir() {
                pending.push(entry);
                continue;
            }
            if entry.file_name().and_then(|name| name.to_str()) == Some("tests.rs") {
                continue;
            }
            let text = std::fs::read_to_string(&entry).expect("the file is readable");
            assert!(
                !text.contains("use crate::verify") && !text.contains("use super::verify"),
                "{} imports the verifier",
                entry.display()
            );
            checked += 1;
        }
    }
    assert!(
        checked >= 4,
        "expected the emitter module tree, found {checked} files"
    );
}
