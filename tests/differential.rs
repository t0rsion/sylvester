//! Both backends against the independent Buchberger oracle in
//! `tests/oracle`.
//!
//! The default suite runs random systems over F_2, F_3, F_5, F_7 (a larger
//! shape), and F_32003, then the two known counterexample systems in every
//! generator order, then a mutation-derived regression fixture. Two sweeps
//! are `#[ignore]`d: an exhaustive small-system sweep and a degree-reversal
//! stress sweep.
//!
//! Seeds are deterministic; failures print the system and both bases.

mod oracle;

use oracle::*;

fn run_differential(p: u64, seeds: std::ops::Range<u64>, shape: Shape) {
    for seed in seeds {
        let mut rng = Rng(seed.wrapping_mul(0x0123_4567_89AB_CDEF).wrapping_add(p));
        let system = random_system(&mut rng, p, shape);
        check_system(&system, p, shape.nvars, seed);
    }
}

#[test]
fn both_backends_match_the_oracle_over_f2() {
    run_differential(2, 0..100, SMALL);
}

#[test]
fn both_backends_match_the_oracle_over_f3() {
    run_differential(3, 0..100, SMALL);
}

#[test]
fn both_backends_match_the_oracle_over_f5() {
    run_differential(5, 0..100, SMALL);
}

#[test]
fn both_backends_match_the_oracle_over_f7_on_the_larger_shape() {
    run_differential(7, 0..25, HARDER);
}

/// The largest test prime, so a coefficient needs 15 bits rather than 3.
#[test]
fn both_backends_match_the_oracle_over_f32003() {
    run_differential(32003, 0..25, SMALL);
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

/// The two hand-verified counterexample systems, in every generator order.
/// The incremental outer loop and POT signatures make backend behavior
/// order-sensitive. The ideal, and thus the oracle, is not.
#[test]
fn both_backends_match_the_oracle_in_every_generator_order() {
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
/// rewriter substitution to pair deletion made the legacy batch engine return a
/// basis of the wrong ideal on exactly this F_3 system (found by mutation
/// testing; the degree-batch sweep below, seed 559).
#[test]
fn pair_deletion_instead_of_rewriter_substitution_returns_the_wrong_ideal() {
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
/// The sweep takes minutes. `cargo test --release --test differential -- --ignored` runs it.
#[test]
#[ignore = "exhaustive sweep; minutes of runtime, run with --ignored"]
fn both_backends_match_the_oracle_on_every_f2_triple() {
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
fn both_backends_match_the_oracle_across_degree_batch_boundaries() {
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
