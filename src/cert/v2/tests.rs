use super::*;
use crate::compute::ComputeOptions;
use crate::compute::f4::{Limits, solve_recorded};

/// The four certificates that section 11 pins byte for byte.
#[path = "../../../tests/support/cert_v2.rs"]
mod fixtures;

use fixtures::{EMPTY, SQUARE, TINY, UNIT, ZERO};

fn ring() -> PolynomialRing {
    PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime")
}

fn polys(ring: &PolynomialRing, texts: &[&str]) -> Vec<Polynomial> {
    texts
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the syntax holds"))
        .collect()
}

/// Run F4 with the recorder and write the certificate.
fn certificate(ring: &PolynomialRing, input: &[Polynomial]) -> Vec<u8> {
    let (bytes, _) = write_certificate(ring, input);
    bytes
}

fn write_certificate(ring: &PolynomialRing, input: &[Polynomial]) -> (Vec<u8>, Vec<Polynomial>) {
    let mut recorder = Recorder::new(ring, input, None, None);
    let options = ComputeOptions::new();
    let (basis, _) = solve_recorded(ring, input, &options, &Limits::of(&options), &mut recorder)
        .expect("the run finishes");
    let (bytes, _) = assemble(ring, input, &basis, recorder).expect("the writer holds");
    (bytes, basis)
}

fn hex(bytes: &[u8]) -> String {
    bytes
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<Vec<_>>()
        .join(" ")
}

#[test]
fn the_square_fixture_is_written_byte_for_byte() {
    let ring = ring();
    let input = polys(&ring, &["x^2*y", "x*y^2"]);
    let bytes = certificate(&ring, &input);
    assert_eq!(hex(&bytes), hex(SQUARE));
}

#[test]
fn the_unit_fixture_is_written_byte_for_byte() {
    let ring = ring();
    let input = polys(&ring, &["3"]);
    let bytes = certificate(&ring, &input);
    assert_eq!(hex(&bytes), hex(UNIT));
}

#[test]
fn the_zero_fixture_is_written_byte_for_byte() {
    let ring = ring();
    let input = polys(&ring, &["0"]);
    let bytes = certificate(&ring, &input);
    assert_eq!(hex(&bytes), hex(ZERO));
}

#[test]
fn the_empty_fixture_is_written_byte_for_byte() {
    let ring = ring();
    let bytes = certificate(&ring, &[]);
    assert_eq!(hex(&bytes), hex(EMPTY));
}

#[test]
fn the_tiny_system_writes_the_header_of_the_listing() {
    let ring = ring();
    let input = polys(&ring, &["x^2 + y", "x"]);
    let (bytes, basis) = write_certificate(&ring, &input);
    // The trace depends on how F4 batches the run, so section 11 pins the
    // header and the values, not the bytes.
    let listing: &[u8] = TINY;
    assert_eq!(bytes[..37], listing[..37]);
    assert_eq!(basis.len(), 2);
    assert_eq!(format!("{}", basis[0]), "x");
    assert_eq!(format!("{}", basis[1]), "y");
}

/// The benchmark modulus, so the measurement runs the arithmetic the gate
/// runs.
const P: u64 = 1073741827;

fn wide_ring(nvars: usize) -> PolynomialRing {
    let names: Vec<String> = (0..nvars).map(|index| format!("x{index}")).collect();
    PolynomialRing::prime_field(P, names).expect("the modulus is prime")
}

/// cyclic-n in n variables, as `benchmarks/gb-comparison/gen.py` defines
/// it.
fn cyclic(ring: &PolynomialRing, n: usize) -> Vec<Polynomial> {
    let mut system = Vec::new();
    for d in 1..n {
        let mut terms = Vec::new();
        for i in 0..n {
            let mut exps = vec![0u16; n];
            for j in 0..d {
                exps[(i + j) % n] += 1;
            }
            terms.push((1i64, exps));
        }
        system.push(ring.polynomial(terms).expect("the exponents fit"));
    }
    let last = vec![(1i64, vec![1u16; n]), (-1, vec![0u16; n])];
    system.push(ring.polynomial(last).expect("the exponents fit"));
    system
}

/// katsura-n in n + 1 variables, the POSSO definition.
fn katsura(ring: &PolynomialRing, n: usize) -> Vec<Polynomial> {
    let nv = n + 1;
    let mut system = Vec::new();
    let mut first = Vec::new();
    for i in 0..=n {
        let mut exps = vec![0u16; nv];
        exps[i] = 1;
        first.push((if i == 0 { 1i64 } else { 2 }, exps));
    }
    first.push((-1, vec![0u16; nv]));
    system.push(ring.polynomial(first).expect("the exponents fit"));

    let bound = n as i64;
    for m in 0..n as i64 {
        let mut terms: Vec<(i64, Vec<u16>)> = Vec::new();
        for i in -bound..=bound {
            let j = m - i;
            if j.abs() > bound {
                continue;
            }
            let mut exps = vec![0u16; nv];
            exps[i.unsigned_abs() as usize] += 1;
            exps[j.unsigned_abs() as usize] += 1;
            terms.push((1, exps));
        }
        let mut exps = vec![0u16; nv];
        exps[m as usize] = 1;
        terms.push((-1, exps));
        system.push(ring.polynomial(terms).expect("the exponents fit"));
    }
    system
}

/// The node cap and the budget read the marking pass, not a pruned trace.
///
/// `reachable` gives the kept node count and the kept combination steps
/// before `prune` allocates anything, so a trace above a cap costs the
/// marking pass alone.
#[test]
fn the_kept_node_count_comes_before_the_pruned_trace() {
    let recording = Recording {
        nodes: vec![
            Node::Input { index: 0 },
            Node::Input { index: 1 },
            // Nothing reaches this one.
            Node::Scale { src: 1, scalar: 3 },
            Node::Comb {
                steps: vec![(0, 1), (1, 2)],
            },
        ],
        basis: vec![3],
    };
    let keep = reachable(&recording);
    assert_eq!(keep, vec![true, true, false, true]);
    let kept = keep.iter().filter(|&&mark| mark).count();
    assert_eq!(kept, 3);
    assert!(pruned_bytes(kept, 2, 1) > 0);

    let (nodes, roots, uses) = prune(recording, &keep, kept);
    assert_eq!(
        nodes.len(),
        kept,
        "prune keeps what the marking pass counted"
    );
    assert_eq!(roots, vec![2], "the dropped node renumbers the rest");
    assert_eq!(uses, vec![1, 1, 1]);
}

/// A restart gives the budget back what the dropped attempt held.
///
/// Contract section 9.2: a lane-width restart drops every node of the
/// failed attempt. Their bytes go with them, or the next attempt runs
/// against a budget that still holds them and stops for no reason.
#[test]
fn a_restart_releases_what_the_dropped_nodes_held() {
    use crate::compute::f4::trace::{BatchRow, BatchRows, RowId, Trace};

    let ring = ring();
    let input = polys(&ring, &["x^2*y", "x*y^2"]);
    let mut recorder = Recorder::new(&ring, &input, None, None);
    let seeded = recorder.held_bytes();
    assert!(seeded > 0, "the seed writes one node per generator");

    recorder.rows(&BatchRows {
        rows: &[BatchRow {
            place: RowId::Pivot(0),
            source: 0,
        }],
        mults: &[1, 0],
        nvars: 2,
    });
    assert!(
        recorder.held_bytes() > seeded,
        "the multiplier node costs the budget"
    );

    recorder.restart();
    assert_eq!(
        recorder.held_bytes(),
        seeded,
        "the restart holds the seed and nothing of the dropped attempt"
    );
}

/// What recording costs the run, as a ratio against the raw run.
///
/// Run it with
/// `cargo test --release --lib cert::v2::tests::the_cost_of_recording -- --ignored --nocapture`.
/// The numbers are the minimum of three runs on one machine and one build.
#[test]
#[ignore = "a measurement, not a check"]
fn the_cost_of_recording() {
    use crate::compute::f4::solve;
    use std::time::Instant;

    let cells: Vec<(String, PolynomialRing, Vec<Polynomial>)> = vec![
        {
            let ring = wide_ring(6);
            let system = cyclic(&ring, 6);
            ("cyclic-6".to_string(), ring, system)
        },
        {
            let ring = wide_ring(8);
            let system = katsura(&ring, 7);
            ("katsura-7".to_string(), ring, system)
        },
        {
            let ring = wide_ring(9);
            let system = katsura(&ring, 8);
            ("katsura-8".to_string(), ring, system)
        },
    ];
    println!("| cell | raw s | recorded s | recorded/raw |");
    println!("|---|---|---|---|");
    for (name, ring, system) in &cells {
        let options = ComputeOptions::new();
        let mut raw = f64::MAX;
        let mut recorded = f64::MAX;
        for _ in 0..3 {
            let start = Instant::now();
            solve(ring, system, &options).expect("a basis");
            raw = raw.min(start.elapsed().as_secs_f64());

            let mut recorder = Recorder::new(ring, system, None, None);
            let start = Instant::now();
            solve_recorded(ring, system, &options, &Limits::of(&options), &mut recorder)
                .expect("a basis");
            recorded = recorded.min(start.elapsed().as_secs_f64());
        }
        println!(
            "| {name} | {raw:.4} | {recorded:.4} | {:.2} |",
            recorded / raw
        );
    }
}
