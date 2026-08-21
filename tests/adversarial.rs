//! Adversarial differential and boundary tests against the Buchberger
//! oracle in `tests/oracle`.
//!
//! Five tests run in the default suite. Two sweeps, `adversarial_largest_modulus`
//! and `adversarial_huge_exponents`, are `#[ignore]`d for their runtime; run
//! them with `cargo +1.92 test --release --test adversarial -- --ignored`.

mod oracle;

use oracle::*;
use sylvester::{Backend, ComputeOptions, PolynomialRing};

const P_BIG: u64 = 2_147_483_647; // 2^31 - 1, the largest modulus a ring accepts

/// Shapes the shipped differential sweep never reaches: one and two
/// variables, five variables, higher degrees, and generator counts of two,
/// three, and eight.
#[test]
fn adversarial_shapes() {
    let shapes = [
        shape(1, 6, 3, 1, 4),
        shape(1, 9, 2, 2, 3),
        shape(2, 5, 3, 1, 5),
        shape(2, 6, 2, 2, 4),
        shape(3, 5, 2, 2, 4),
        shape(3, 3, 8, 1, 3),
        shape(4, 2, 8, 1, 3),
        shape(5, 2, 4, 1, 4),
    ];
    for (si, &s) in shapes.iter().enumerate() {
        for p in [2u64, 3, 7, 32003] {
            for seed in 0..40u64 {
                let mut rng = Rng((seed + 1_000 * si as u64)
                    .wrapping_mul(0x0123_4567_89AB_CDEF)
                    .wrapping_add(p));
                let system = random_system(&mut rng, p, s);
                check_system(&system, p, s.nvars, seed + 1_000 * si as u64);
            }
        }
    }
}

/// The largest modulus a ring accepts, where every field product runs at
/// the top of the `u64` range the crate promises fits.
#[test]
#[ignore = "sweep; minutes of runtime, run with --ignored"]
fn adversarial_largest_modulus() {
    for p in [P_BIG, 1_073_741_827] {
        for seed in 0..40u64 {
            let mut rng = Rng(seed.wrapping_mul(0x0123_4567_89AB_CDEF).wrapping_add(p));
            let system = random_system(&mut rng, p, SMALL);
            check_system(&system, p, SMALL.nvars, seed);
        }
        for seed in 0..15u64 {
            let mut rng = Rng(seed.wrapping_mul(0xDEAD_BEEF_1234_5678).wrapping_add(p));
            let system = random_system(&mut rng, p, HARDER);
            check_system(&system, p, HARDER.nvars, seed);
        }
    }
}

/// Degenerate generator lists: zero polynomials mixed in, repeated
/// generators, and a unit that makes the ideal the whole ring.
#[test]
fn adversarial_degenerate_inputs() {
    for p in [2u64, 3, 7, 32003] {
        for seed in 0..60u64 {
            let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(p));
            let base = random_system(&mut rng, p, SMALL);
            if base.iter().all(Poly::is_zero) {
                continue;
            }

            // A zero generator in front, in the middle, and at the end.
            let mut with_zeros = vec![Poly::zero()];
            with_zeros.extend(base.iter().cloned());
            with_zeros.push(Poly::zero());
            with_zeros.insert(2, Poly::zero());
            check_system(&with_zeros, p, SMALL.nvars, seed);

            // Every generator repeated.
            let mut doubled: Vec<Poly> = Vec::new();
            for f in &base {
                doubled.push(f.clone());
                doubled.push(f.clone());
            }
            check_system(&doubled, p, SMALL.nvars, seed);

            // A unit somewhere in the list.
            let mut one = Poly::zero();
            one.terms.insert(vec![0u16; SMALL.nvars], 1);
            let mut with_unit = base.clone();
            with_unit.insert(base.len() / 2, one);
            check_system(&with_unit, p, SMALL.nvars, seed);
        }
    }
}

/// A single generator, and a ring with no variables at all.
#[test]
fn adversarial_tiny_inputs() {
    for p in [2u64, 5, 32003] {
        for seed in 0..40u64 {
            let mut rng = Rng(seed.wrapping_mul(0xA24B_AED4_963E_E407).wrapping_add(p));
            let s = shape(2, 4, 1, 1, 4);
            let system = random_system(&mut rng, p, s);
            check_system(&system, p, s.nvars, seed);
        }

        // No variables: the ring is the field itself.
        let ring = PolynomialRing::prime_field(p, Vec::<String>::new()).expect("prime");
        let zero = ring
            .polynomial(Vec::<(i64, Vec<u16>)>::new())
            .expect("fits");
        let one = ring.polynomial([(1i64, Vec::<u16>::new())]).expect("fits");
        for gens in [
            vec![],
            vec![zero.clone()],
            vec![one.clone()],
            vec![zero, one],
        ] {
            let ideal = ring.ideal(gens.clone()).expect("same ring");
            for backend in [Backend::Classic, Backend::Matrix] {
                let out = ideal
                    .groebner_basis(ComputeOptions::new().backend(backend))
                    .expect("no budget");
                let expect_unit = gens.iter().any(|f| !f.is_zero());
                assert_eq!(
                    out.len(),
                    usize::from(expect_unit),
                    "p={p} backend={backend:?} gens={gens:?}"
                );
            }
        }
    }
}

/// Both backends against each other on shapes the oracle is too slow for.
#[test]
fn adversarial_backend_agreement_on_larger_shapes() {
    let shapes = [
        shape(4, 3, 4, 2, 4),
        shape(5, 2, 5, 2, 4),
        shape(6, 2, 4, 1, 3),
    ];
    for (si, &s) in shapes.iter().enumerate() {
        for p in [2u64, 32003] {
            for seed in 0..12u64 {
                let mut rng = Rng((seed + 7_000 * si as u64)
                    .wrapping_mul(0x0123_4567_89AB_CDEF)
                    .wrapping_add(p));
                let system = random_system(&mut rng, p, s);
                if system.iter().all(Poly::is_zero) {
                    continue;
                }
                let ideal = to_engine_ideal(&system, p, s.nvars);
                let classic = ideal
                    .groebner_basis(ComputeOptions::new().backend(Backend::Classic))
                    .expect("no budget");
                let matrix = ideal
                    .groebner_basis(ComputeOptions::new().backend(Backend::Matrix))
                    .expect("no budget");
                let a = canonical_set(&engine_output_to_polys(&classic, p, s.nvars), p);
                let b = canonical_set(&engine_output_to_polys(&matrix, p, s.nvars), p);
                assert_eq!(
                    a,
                    b,
                    "backends disagree, shape {si}, p={p}, seed={seed}, system {}",
                    render_system(&system)
                );
            }
        }
    }
}

/// Exponents near `u16::MAX`, where a product can leave the width one
/// exponent holds. Every run must end in a value or in a typed error, never
/// in a panic.
#[test]
#[ignore = "sweep; minutes of runtime, run with --ignored"]
fn adversarial_huge_exponents() {
    const SKIP_OVERSIZED: bool = true;
    let pool: [u16; 8] = [0, 1, 2, 3, 21845, 32767, 40000, 65535];
    for nvars in 1..=3usize {
        for p in [2u64, 32003] {
            for seed in 0..120u64 {
                let mut rng = Rng(seed.wrapping_mul(0x2545_F491_4F6C_DD1D).wrapping_add(p * 7));
                let ngens = 2 + rng.below(2) as usize;
                let ring = PolynomialRing::prime_field(p, variable_names(nvars))
                    .expect("the modulus is prime");
                let generators: Vec<_> = (0..ngens)
                    .map(|_| {
                        let nterms = 1 + rng.below(3);
                        let terms: Vec<(i64, Vec<u16>)> = (0..nterms)
                            .map(|_| {
                                let exps: Vec<u16> = (0..nvars)
                                    .map(|_| pool[rng.below(pool.len() as u64) as usize])
                                    .collect();
                                (1 + rng.below(p - 1) as i64, exps)
                            })
                            .collect();
                        ring.polynomial(terms).expect("the exponents fit the ring")
                    })
                    .collect();
                // Skip any input whose own total degree is past the width
                // one exponent holds: that is the case under test elsewhere.
                if SKIP_OVERSIZED
                    && generators
                        .iter()
                        .any(|f| f.degree().is_some_and(|d| d > 65535))
                {
                    continue;
                }
                let ideal = ring.ideal(generators).expect("one ring");
                // A few near-limit configurations run a genuinely large
                // computation, so each carries a timeout. The claim under
                // test is that every run ends in a value or a typed error,
                // never a panic; a stopped run reports the typed error.
                let options = ComputeOptions::new().timeout(std::time::Duration::from_millis(750));
                for backend in [Backend::Classic, Backend::Matrix] {
                    match ideal.groebner_basis(options.clone().backend(backend)) {
                        Ok(_)
                        | Err(sylvester::ComputeError::DegreeLimit { .. })
                        | Err(sylvester::ComputeError::Timeout) => {}
                        Err(other) => panic!(
                            "unexpected error {other} for nvars={nvars} p={p} seed={seed} backend={backend:?}"
                        ),
                    }
                }
            }
        }
    }
}

/// The certified path must return the same basis the oracle does.
#[test]
fn adversarial_certified_matches_the_oracle() {
    for p in [2u64, 3, 7, 32003] {
        for seed in 0..40u64 {
            let mut rng = Rng(seed.wrapping_mul(0x0123_4567_89AB_CDEF).wrapping_add(p));
            let system = random_system(&mut rng, p, SMALL);
            if system.iter().all(Poly::is_zero) {
                continue;
            }
            let oracle = buchberger(&system, p);
            let oracle_reduced = canonical_set(&reduce_basis(&oracle, p), p);
            let ideal = to_engine_ideal(&system, p, SMALL.nvars);
            let certified = ideal
                .groebner_basis_certified(ComputeOptions::new())
                .unwrap_or_else(|e| {
                    panic!(
                        "certified error {e} on seed {seed}, p={p}, {}",
                        render_system(&system)
                    )
                });
            let got = canonical_set(
                &engine_output_to_polys(certified.basis(), p, SMALL.nvars),
                p,
            );
            assert_eq!(
                got,
                oracle_reduced,
                "certified mismatch on seed {seed}, p={p}, {}",
                render_system(&system)
            );
        }
    }
}
