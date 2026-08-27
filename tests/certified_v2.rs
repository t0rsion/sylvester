//! End-to-end tests for the `sylv-gb-cert-v2` certified path.
//!
//! Every run here goes through `Ideal::groebner_basis_certified` with the
//! F4 backend, which writes a v2 certificate and hands it to the
//! independent verifier. The certificate holds by construction, so a
//! failure is an engine or a writer defect, not a disagreement about the
//! contract.
//!
//! The five fixtures under `tests/fixtures/cert-v2` are the test vectors
//! of `docs/certificate-v2.md` section 11.

use std::collections::BTreeMap;
use std::time::{Duration, Instant};

use sylvester::verify::verify;
use sylvester::{
    Backend, CertifyError, ComputeError, ComputeOptions, GroebnerBasis, Ideal, PolynomialRing,
    RingError,
};

#[path = "support/cert_v2.rs"]
mod fixtures;

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

/// The basis as one sorted string per element, so two bases compare as
/// sets.
fn canonical(basis: &GroebnerBasis) -> Vec<String> {
    let mut out: Vec<String> = basis.iter().map(|poly| poly.to_string()).collect();
    out.sort();
    out
}

/// The bytes of one fixture.
fn fixture(name: &str) -> Vec<u8> {
    match name {
        "tiny" => fixtures::TINY,
        "square" => fixtures::SQUARE,
        "unit" => fixtures::UNIT,
        "zero" => fixtures::ZERO,
        "empty" => fixtures::EMPTY,
        other => panic!("unknown v2 fixture {other}"),
    }
    .to_vec()
}

/// The certificate of a system over `F_7[x, y]`, written from its terms.
fn small(texts: &[&str]) -> Vec<u8> {
    let ring = PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime");
    let generators: Vec<_> = texts
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the syntax holds"))
        .collect();
    let ideal = ring
        .ideal(generators)
        .expect("the polynomials share a ring");
    ideal
        .groebner_basis_certified(ComputeOptions::new())
        .expect("the certificate holds")
        .into_parts()
        .1
}

fn hex(bytes: &[u8]) -> String {
    bytes
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Certifies one system and compares the basis against the raw F4 run.
fn certify(label: &str, ideal: &Ideal) -> Vec<u8> {
    let raw = ideal
        .groebner_basis(ComputeOptions::new())
        .unwrap_or_else(|error| panic!("{label}: {error}"));
    let certified = ideal
        .groebner_basis_certified(ComputeOptions::new())
        .unwrap_or_else(|error| panic!("{label}: {error}"));
    assert_eq!(
        canonical(certified.basis()),
        canonical(&raw),
        "{label}: the certified basis is the raw basis"
    );
    // The verifier already accepted these bytes inside the call. Reading
    // them back from the outside is the same check a caller can run on a
    // certificate it received.
    let verified = verify(certified.certificate())
        .unwrap_or_else(|error| panic!("{label}: the verifier rejected the bytes: {error}"));
    assert_eq!(verified.basis().len(), raw.len(), "{label}: basis length");
    certified.into_parts().1
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

/// The two systems of KNOWN_ISSUES.md, which the v0.1 engines got wrong.
fn counterexamples() -> [(u64, usize, Vec<Terms>); 2] {
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
    [(2, 2, f2), (3, 3, f3)]
}

#[test]
fn the_verifier_accepts_the_five_contract_fixtures() {
    for name in ["tiny", "square", "unit", "zero", "empty"] {
        let bytes = fixture(name);
        verify(&bytes).unwrap_or_else(|error| panic!("{name}: {error}"));
    }
}

#[test]
fn the_writer_reproduces_the_four_pinned_fixtures() {
    // Section 11 pins these four: each trace holds at most two nodes, and
    // section 9 fixes every one of them.
    assert_eq!(hex(&small(&["x^2*y", "x*y^2"])), hex(&fixture("square")));
    assert_eq!(hex(&small(&["3"])), hex(&fixture("unit")));
    assert_eq!(hex(&small(&["0"])), hex(&fixture("zero")));
    assert_eq!(hex(&small(&[])), hex(&fixture("empty")));
}

#[test]
fn the_tiny_system_matches_the_listed_fields() {
    // The listing of section 10 is a hand-built derivation, so the writer's
    // bytes are not pinned to it. The header and the values are.
    let listing = fixture("tiny");
    let listed = verify(&listing).expect("the listing holds");
    let bytes = small(&["x^2 + y", "x"]);
    let written = verify(&bytes).expect("the certificate holds");
    assert_eq!(bytes[..37], listing[..37]);
    assert_eq!(written.modulus(), listed.modulus());
    assert_eq!(written.nvars(), listed.nvars());
    assert_eq!(written.input(), listed.input());
    assert_eq!(written.basis(), listed.basis());
}

#[test]
fn two_runs_write_the_same_bytes() {
    let ideal = ideal_of(&ring(P, 5), &cyclic(5));
    let first = certify("cyclic-5", &ideal);
    let second = certify("cyclic-5", &ideal);
    assert_eq!(first, second);
}

#[test]
fn the_thread_count_changes_no_byte() {
    // The recorder reduces every row on the calling thread, so the trace
    // is the same whatever the pool holds.
    let ideal = ideal_of(&ring(P, 5), &katsura(4));
    let one = ideal
        .groebner_basis_certified(ComputeOptions::new().threads(1))
        .expect("the certificate holds");
    let eight = ideal
        .groebner_basis_certified(ComputeOptions::new().threads(8))
        .expect("the certificate holds");
    assert_eq!(one.certificate(), eight.certificate());
}

#[test]
fn the_classic_backend_still_writes_a_v1_certificate() {
    let ideal = ideal_of(&ring(P, 4), &cyclic(4));
    let classic = ideal
        .groebner_basis_certified(ComputeOptions::new().backend(Backend::Classic))
        .expect("the certificate holds");
    assert_eq!(classic.certificate()[0], b'{');
    let f4 = ideal
        .groebner_basis_certified(ComputeOptions::new())
        .expect("the certificate holds");
    assert_eq!(&f4.certificate()[..6], b"SYLVGB");
    assert_eq!(canonical(classic.basis()), canonical(f4.basis()));
}

#[test]
fn certifies_the_counterexample_systems() {
    for (modulus, nvars, system) in counterexamples() {
        let ring = ring(modulus, nvars);
        certify(
            &format!("F_{modulus} counterexample"),
            &ideal_of(&ring, &system),
        );
    }
}

#[test]
fn certifies_cyclic() {
    for n in 4..=6 {
        certify(&format!("cyclic-{n}"), &ideal_of(&ring(P, n), &cyclic(n)));
    }
}

#[test]
fn certifies_katsura() {
    for n in 4..=7 {
        let system = katsura(n);
        certify(&format!("katsura-{n}"), &ideal_of(&ring(P, n + 1), &system));
    }
}

#[test]
fn certifies_noon() {
    for n in 3..=5 {
        certify(&format!("noon-{n}"), &ideal_of(&ring(P, n), &noon(n)));
    }
}

#[test]
fn certifies_eco_8() {
    certify("eco-8", &ideal_of(&ring(P, 8), &eco(8)));
}

#[test]
fn certifies_a_run_that_restarts_at_a_wider_lane_width() {
    // The product of the two leads passes the 127 bound of the narrowest
    // packing, so the run starts again wider. The writer drops the nodes
    // of the attempt that overflowed and records the new one from the
    // start.
    let ring = ring(P, 2);
    let system = vec![
        vec![(1i64, vec![100u16, 30]), (-1, vec![0, 0])],
        vec![(1i64, vec![30u16, 100]), (-1, vec![0, 0])],
    ];
    let ideal = ideal_of(&ring, &system);
    let free = certify("wide exponents", &ideal);

    // The attempt that overflowed wrote nodes the restart drops. Their
    // bytes go back to the budget, so a limit the finished attempt fits
    // under certifies the run.
    let limit = smallest_limit(|limit| {
        ideal
            .groebner_basis_certified(ComputeOptions::new().memory_limit(limit))
            .is_ok()
    });
    let bounded = ideal
        .groebner_basis_certified(ComputeOptions::new().memory_limit(limit))
        .expect("the restart releases what the dropped attempt held");
    assert_eq!(free, bounded.certificate(), "the limit changes no byte");
}

/// One memory limit covers the engine, the writer, and the verifier.
///
/// Self-calibrating: the search finds the smallest limit the certified run
/// finishes under. The run finishes there and writes the bytes the free
/// run writes.
#[test]
fn one_limit_covers_the_engine_the_writer_and_the_verifier() {
    let ideal = ideal_of(&ring(P, 5), &cyclic(5));
    let free = certify("cyclic-5", &ideal);
    let limit = smallest_limit(|limit| {
        ideal
            .groebner_basis_certified(ComputeOptions::new().memory_limit(limit))
            .is_ok()
    });
    let bounded = ideal
        .groebner_basis_certified(ComputeOptions::new().memory_limit(limit))
        .expect("the smallest limit that fits holds the whole run");
    assert_eq!(free, bounded.certificate());
    // One byte under what the run needs stops it typed. Which part runs
    // out first is the budget's business, so the check names the class and
    // not the part.
    let tight = ideal.groebner_basis_certified(ComputeOptions::new().memory_limit(limit - 1));
    assert!(
        matches!(
            tight,
            Err(CertifyError::Engine(ComputeError::MemoryLimitExceeded))
                | Err(CertifyError::VerifierExhausted(_))
        ),
        "one byte under what the run needs: {tight:?}"
    );
}

#[test]
fn certifies_a_unit_basis_a_batch_found() {
    // Seeding finds no constant here. The pair of the two generators
    // reduces to one, so the constant is a pivot of a batch and the
    // trace names it by its pivot slot.
    let ring = ring(P, 2);
    let system = vec![
        vec![(1i64, vec![1u16, 0])],
        vec![(1i64, vec![1u16, 0]), (1, vec![0, 0])],
    ];
    let ideal = ideal_of(&ring, &system);
    let certified = certify("unit from a batch", &ideal);
    let verified = verify(&certified).expect("the certificate holds");
    assert_eq!(verified.basis().len(), 1);
    assert_eq!(verified.basis()[0].terms().len(), 1);
    assert_eq!(verified.basis()[0].terms()[0].coeff(), 1);
    assert_eq!(verified.basis()[0].terms()[0].mono().degree(), 0);
}

#[test]
fn certifies_an_input_with_zero_duplicate_and_scaled_generators() {
    // The trace carries the normalization: the scale that makes a
    // generator monic is a node, and dropping a zero or a duplicate only
    // changes which node the later references name.
    let ring = PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime");
    let texts = ["2*x^2 + 2*y", "0", "x^2 + y", "x"];
    let generators: Vec<_> = texts
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the syntax holds"))
        .collect();
    let ideal = ring
        .ideal(generators)
        .expect("the polynomials share a ring");
    let certified = ideal
        .groebner_basis_certified(ComputeOptions::new())
        .expect("the certificate holds");
    assert_eq!(canonical(certified.basis()), ["x", "y"]);
    let verified = verify(certified.certificate()).expect("the certificate holds");
    assert_eq!(verified.input().len(), 4, "the input keeps every generator");
}

#[test]
fn an_exhausted_timeout_stops_the_certified_path() {
    let ideal = ideal_of(&ring(P, 4), &cyclic(4));
    assert_eq!(
        ideal.groebner_basis_certified(ComputeOptions::new().timeout(Duration::ZERO)),
        Err(CertifyError::Engine(ComputeError::Timeout))
    );
}

#[test]
fn a_tight_memory_limit_stops_the_certified_path() {
    // One node costs more than one byte, so the recorder stops on the
    // first node it writes and the writer writes no certificate.
    let ideal = ideal_of(&ring(P, 4), &cyclic(4));
    assert_eq!(
        ideal.groebner_basis_certified(ComputeOptions::new().memory_limit(1)),
        Err(CertifyError::Engine(ComputeError::MemoryLimitExceeded))
    );
}

#[test]
fn a_memory_limit_that_fits_changes_no_byte() {
    let ideal = ideal_of(&ring(P, 5), &cyclic(5));
    let free = certify("cyclic-5", &ideal);
    let bounded = ideal
        .groebner_basis_certified(ComputeOptions::new().memory_limit(1 << 30))
        .expect("the certificate holds");
    assert_eq!(free, bounded.certificate());
}

#[test]
fn a_tampered_fixture_is_rejected_or_still_true() {
    // The tamper corpus of section 11, from the writer's side. The
    // verifier's own suite carries the structural cases; this one covers
    // what the writer can state without a second decoder.
    //
    // A flipped byte is not always a rejection. Byte 36 is the variable
    // count, and every monomial here is a sparse support, so the same
    // derivation holds over a ring of more variables. Byte 52 of `square`
    // is an exponent, and `[x^3*y, x*y^2]` is the reduced basis of its own
    // ideal as much as `[x^2*y, x*y^2]` is. So the sweep requires no
    // panic, and the cases below require rejection.
    for name in ["tiny", "square", "unit", "zero", "empty"] {
        let bytes = fixture(name);
        tamper_each_byte(&bytes);
        reject_prefixes(name, &bytes);
        reject_trailing_byte(name, &bytes);
        reject_fixed_header(name, &bytes);
        reject_section_lengths(name, &bytes);
        reject_composite_modulus(name, &bytes);
    }
}

fn tamper_each_byte(bytes: &[u8]) {
    for offset in 0..bytes.len() {
        for mask in [0x01u8, 0x80, 0xff] {
            let mut broken = bytes.to_vec();
            broken[offset] ^= mask;
            std::hint::black_box(verify(&broken).is_err());
        }
    }
}

fn reject_prefixes(name: &str, bytes: &[u8]) {
    for offset in 0..bytes.len() {
        assert!(
            verify(&bytes[..offset]).is_err(),
            "{name}: the first {offset} bytes were accepted"
        );
    }
}

fn reject_trailing_byte(name: &str, bytes: &[u8]) {
    let mut longer = bytes.to_vec();
    longer.push(0);
    assert!(
        verify(&longer).is_err(),
        "{name}: a trailing byte was accepted"
    );
}

fn reject_fixed_header(name: &str, bytes: &[u8]) {
    for offset in 0..37 {
        if offset == 35 || offset == 36 {
            continue;
        }
        let mut broken = bytes.to_vec();
        broken[offset] ^= 0x01;
        assert!(
            verify(&broken).is_err(),
            "{name}: the header byte {offset} was accepted"
        );
    }
}

fn reject_section_lengths(name: &str, bytes: &[u8]) {
    for offset in 37..43 {
        let mut broken = bytes.to_vec();
        broken[offset] = broken[offset].wrapping_add(1);
        assert!(
            verify(&broken).is_err(),
            "{name}: the section length at {offset} was accepted"
        );
    }
}

fn reject_composite_modulus(name: &str, bytes: &[u8]) {
    let mut composite = bytes.to_vec();
    composite[35] = 9;
    assert!(
        verify(&composite).is_err(),
        "{name}: a composite modulus was accepted"
    );
}

/// The cost of certification, as a table.
///
/// `cargo test --release --test certified_v2 -- --ignored --nocapture`
/// prints the table. The numbers are the minimum of three runs, so they
/// belong to one machine and one build, and the release note quotes the
/// record they come from.
#[test]
#[ignore = "a measurement, not a check"]
fn the_cost_of_certification() {
    let cells: Vec<(String, Ideal)> = vec![
        ("cyclic-6".to_string(), ideal_of(&ring(P, 6), &cyclic(6))),
        ("katsura-7".to_string(), ideal_of(&ring(P, 8), &katsura(7))),
        ("katsura-8".to_string(), ideal_of(&ring(P, 9), &katsura(8))),
        ("cyclic-7".to_string(), ideal_of(&ring(P, 7), &cyclic(7))),
        ("eco-9".to_string(), ideal_of(&ring(P, 9), &eco(9))),
        ("noon-6".to_string(), ideal_of(&ring(P, 6), &noon(6))),
    ];
    println!("| cell | raw s | certified s | verify s | bytes | certified/raw | verify/raw |");
    println!("|---|---|---|---|---|---|---|");
    for (name, ideal) in &cells {
        let mut raw = f64::MAX;
        let mut certified = f64::MAX;
        let mut checked = f64::MAX;
        let mut bytes = 0;
        for _ in 0..3 {
            let start = Instant::now();
            let basis = ideal
                .groebner_basis(ComputeOptions::new())
                .expect("a basis");
            raw = raw.min(start.elapsed().as_secs_f64());

            let start = Instant::now();
            let value = ideal
                .groebner_basis_certified(ComputeOptions::new())
                .expect("a certificate");
            certified = certified.min(start.elapsed().as_secs_f64());
            assert_eq!(canonical(value.basis()), canonical(&basis));

            let start = Instant::now();
            verify(value.certificate()).expect("the verifier accepts");
            checked = checked.min(start.elapsed().as_secs_f64());
            bytes = value.certificate().len();
        }
        println!(
            "| {name} | {raw:.3} | {certified:.3} | {checked:.3} | {bytes} | {:.1} | {:.1} |",
            certified / raw,
            checked / raw
        );
    }
}

/// The smallest memory limit `run` finishes under, by binary search.
///
/// The search assumes what the budget holds to: a run that finishes under
/// one limit finishes under every larger one.
fn smallest_limit(run: impl Fn(usize) -> bool) -> usize {
    let mut high = 1 << 12;
    while !run(high) {
        high *= 2;
        assert!(high <= 1 << 32, "no limit under 4 GiB finishes the run");
    }
    let mut low = 0;
    while high - low > 1 {
        let middle = low + (high - low) / 2;
        if run(middle) {
            high = middle;
        } else {
            low = middle;
        }
    }
    high
}

/// The writer stops typed on a limit the engine alone fits under.
///
/// Self-calibrating: the search finds the smallest limit the raw run
/// finishes under. The certificate is more data than the run that wrote
/// it, so the same limit stops the certified path, and it stops with a
/// typed budget error and no bytes.
#[test]
fn a_limit_the_engine_fits_but_the_certificate_does_not_stops_typed() {
    let ideal = ideal_of(&ring(P, 5), &cyclic(5));
    let raw_limit = smallest_limit(|limit| {
        ideal
            .groebner_basis(ComputeOptions::new().memory_limit(limit))
            .is_ok()
    });
    assert!(
        ideal
            .groebner_basis(ComputeOptions::new().memory_limit(raw_limit))
            .is_ok()
    );
    // The engine's own estimate counts the recorder's nodes, so which
    // part of the run notices first is the budget's business. The check
    // names the class.
    let stopped = ideal.groebner_basis_certified(ComputeOptions::new().memory_limit(raw_limit));
    assert!(
        matches!(
            stopped,
            Err(CertifyError::Engine(ComputeError::MemoryLimitExceeded))
        ),
        "the writer needs more than the engine: {stopped:?}"
    );
}
