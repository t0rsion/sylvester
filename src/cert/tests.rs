//! Emission tests: every certificate is checked byte for byte.
//!
//! The expected strings are the ones `tests/verifier.rs` accepts by hand, so
//! the writer and the hand-written corpus stay one format. The verifier runs
//! over them there, not here: this module never imports `crate::verify`, and
//! `tests/certified.rs` covers the emitter and the verifier end to end.

use std::time::{Duration, Instant};

use super::{Budget, assemble, divide::divide, origin::Origin, represent};
use crate::certificate::{CertifyError, EmitterFault, Place};
use crate::compute::ComputeError;
use crate::poly::Polynomial;
use crate::ring::PolynomialRing;

const P: u64 = 7;

/// A budget that stops nothing, for the two-variable ring.
fn budget() -> Budget {
    Budget::unlimited(2)
}

fn ring() -> PolynomialRing {
    PolynomialRing::prime_field(P, ["x", "y"]).expect("7 is prime")
}

/// Assemble over the two-variable ring with an unlimited budget.
fn emit(
    input: &[Polynomial],
    basis: &[Polynomial],
    origins: &[Origin],
) -> Result<Vec<u8>, CertifyError> {
    assemble(
        &ring(),
        input,
        basis.to_vec(),
        origins.to_vec(),
        &mut budget(),
    )
}

/// Build polynomials of the two-variable ring from text.
fn list(texts: &[&str]) -> Vec<Polynomial> {
    let ring = ring();
    texts
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the text parses"))
        .collect()
}

fn one(text: &str) -> Polynomial {
    list(&[text]).pop().expect("the list holds one polynomial")
}

/// The product of two polynomials, over the engine arithmetic.
fn product(a: &Polynomial, b: &Polynomial) -> Polynomial {
    let mut out = a.zero_like();
    for term in &a.terms {
        out = out.add(
            &b.scale_monomial(term.coeff, &term.mono, P)
                .expect("the product fits"),
            P,
        );
    }
    out
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8(bytes.to_vec()).expect("the certificate is UTF-8")
}

/// F = {x}, G = {x}, with x = 1 * x.
#[test]
fn emits_the_single_generator_certificate() {
    let input = list(&["x"]);
    let basis = input.clone();
    let origins = vec![list(&["1"])];

    let bytes = emit(&input, &basis, &origins).expect("the basis holds");
    assert_eq!(
        text(&bytes),
        concat!(
            r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
            r#""input":[[[1,[1,0]]]],"basis":[[[1,[1,0]]]],"#,
            r#""origin":[[[[1,[0,0]]]]],"membership":[[[[1,[0,0]]]]],"spairs":[]}"#
        )
    );
}

/// F = {x, x + 1}, G = {1}, with 1 = -x + (x + 1).
#[test]
fn emits_the_unit_ideal_certificate() {
    let input = list(&["x", "x + 1"]);
    let basis = list(&["1"]);
    let origins = vec![list(&["-1", "1"])];

    let bytes = emit(&input, &basis, &origins).expect("the basis holds");
    assert_eq!(
        text(&bytes),
        concat!(
            r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
            r#""input":[[[1,[1,0]]],[[1,[1,0]],[1,[0,0]]]],"basis":[[[1,[0,0]]]],"#,
            r#""origin":[[[[6,[0,0]]],[[1,[0,0]]]]],"#,
            r#""membership":[[[[1,[1,0]]]],[[[1,[1,0]],[1,[0,0]]]]],"spairs":[]}"#
        )
    );
}

/// F = {x^2 - 1, xy - 1}, G = [y^2 - 1, x - y]. The only pair is coprime,
/// so the certificate carries no S-pair entry.
fn pair_system() -> (Vec<Polynomial>, Vec<Polynomial>, Vec<Origin>) {
    let input = list(&["x^2 - 1", "x*y - 1"]);
    let basis = list(&["y^2 - 1", "x - y"]);
    let origins = vec![list(&["-y^2", "x*y + 1"]), list(&["y", "-x"])];
    (input, basis, origins)
}

const PAIR: &str = concat!(
    r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
    r#""input":[[[1,[2,0]],[6,[0,0]]],[[1,[1,1]],[6,[0,0]]]],"#,
    r#""basis":[[[1,[0,2]],[6,[0,0]]],[[1,[1,0]],[6,[0,1]]]],"#,
    r#""origin":[[[[6,[0,2]]],[[1,[1,1]],[1,[0,0]]]],[[[1,[0,1]]],[[6,[1,0]]]]],"#,
    r#""membership":[[[[1,[0,0]]],[[1,[1,0]],[1,[0,1]]]],[[[1,[0,0]]],[[1,[0,1]]]]],"#,
    r#""spairs":[]}"#
);

#[test]
fn emits_a_two_element_basis_with_hand_derived_origins() {
    let (input, basis, origins) = pair_system();
    let bytes = emit(&input, &basis, &origins).expect("the basis holds");
    assert_eq!(text(&bytes), PAIR);
}

#[test]
fn sorts_the_basis_and_carries_the_origins_with_it() {
    let (input, basis, origins) = pair_system();
    let swapped: Vec<Polynomial> = basis.iter().rev().cloned().collect();
    let swapped_origins: Vec<Origin> = origins.iter().rev().cloned().collect();

    let bytes = emit(&input, &swapped, &swapped_origins).expect("the basis holds");
    assert_eq!(text(&bytes), PAIR);
}

/// F = {x^2 - y, xy - 1}, G = [x^2 - y, xy - 1, y^2 - x]. The pairs (0,1)
/// and (1,2) share a variable; the pair (0,2) is coprime.
fn base_system() -> (Vec<Polynomial>, Vec<Polynomial>, Vec<Origin>) {
    let input = list(&["x^2 - y", "x*y - 1"]);
    let basis = list(&["x^2 - y", "x*y - 1", "y^2 - x"]);
    let origins = vec![list(&["1", "0"]), list(&["0", "1"]), list(&["-y", "x"])];
    (input, basis, origins)
}

const BASE: &str = concat!(
    r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
    r#""input":[[[1,[2,0]],[6,[0,1]]],[[1,[1,1]],[6,[0,0]]]],"#,
    r#""basis":[[[1,[2,0]],[6,[0,1]]],[[1,[1,1]],[6,[0,0]]],[[1,[0,2]],[6,[1,0]]]],"#,
    r#""origin":[[[[1,[0,0]]],[]],[[],[[1,[0,0]]]],[[[6,[0,1]]],[[1,[1,0]]]]],"#,
    r#""membership":[[[[1,[0,0]]],[],[]],[[],[[1,[0,0]]],[]]],"#,
    r#""spairs":[[0,1,[[],[],[[6,[0,0]]]]],[1,2,[[[1,[0,0]]],[],[]]]]}"#
);

#[test]
fn emits_the_spair_entries_of_a_three_element_basis() {
    let (input, basis, origins) = base_system();
    let bytes = emit(&input, &basis, &origins).expect("the basis holds");
    assert_eq!(text(&bytes), BASE);
}

#[test]
fn leaves_out_the_pair_with_coprime_leading_monomials() {
    let (_, basis, _) = base_system();
    let entries = represent::spair_representations(&basis, P, &mut budget())
        .expect("the basis is a Gröbner basis");
    let pairs: Vec<(usize, usize)> = entries.iter().map(|entry| (entry.i, entry.j)).collect();
    assert_eq!(pairs, vec![(0, 1), (1, 2)]);
    assert_eq!(entries[0].cofactors.len(), basis.len());
}

/// F = {0}, G = {}: the zero ideal.
#[test]
fn emits_the_zero_ideal_certificate() {
    let ring = ring();
    // The coefficient 7 reduces to zero, so the polynomial vanishes and
    // keeps its place.
    let input = vec![ring.polynomial([(7, [0, 0])]).expect("fits")];
    assert!(input[0].is_zero());

    let bytes = assemble(&ring, &input, Vec::new(), Vec::new(), &mut budget())
        .expect("the empty basis holds");
    assert_eq!(
        text(&bytes),
        concat!(
            r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
            r#""input":[[]],"basis":[],"origin":[],"membership":[[]],"spairs":[]}"#
        )
    );
}

/// F = {1}, G = {1} with no variables.
#[test]
fn emits_a_certificate_without_variables() {
    let ring = PolynomialRing::prime_field(P, Vec::<String>::new()).expect("7 is prime");
    let input = vec![ring.one()];
    let basis = input.clone();
    let origins = vec![input.clone()];

    let bytes = assemble(&ring, &input, basis.clone(), origins.clone(), &mut budget())
        .expect("the basis holds");
    assert_eq!(
        text(&bytes),
        concat!(
            r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":0,"#,
            r#""input":[[[1,[]]]],"basis":[[[1,[]]]],"#,
            r#""origin":[[[[1,[]]]]],"membership":[[[[1,[]]]]],"spairs":[]}"#
        )
    );
}

#[test]
fn two_assemblies_of_one_system_give_the_same_bytes() {
    let (input, basis, origins) = base_system();
    let first = emit(&input, &basis, &origins).expect("the basis holds");
    let second = emit(&input, &basis, &origins).expect("the basis holds");
    assert_eq!(first, second);

    let (pair_input, pair_basis, pair_origins) = pair_system();
    let third = emit(&pair_input, &pair_basis, &pair_origins).expect("the basis holds");
    let fourth = emit(&pair_input, &pair_basis, &pair_origins).expect("the basis holds");
    assert_eq!(third, fourth);
    assert_ne!(first, third);
}

#[test]
fn a_wrong_origin_cofactor_goes_into_the_bytes_unrepaired() {
    let (input, basis, mut origins) = pair_system();
    // The cofactor of x - y over xy - 1 is -x. The sign flip keeps the
    // shape and breaks the identity. The writer must copy it out as it
    // stands; `tests/certified.rs` checks that such bytes are rejected.
    origins[1] = list(&["y", "x"]);

    let bytes = emit(&input, &basis, &origins).expect("the emitter writes what it gets");
    assert_eq!(
        text(&bytes),
        PAIR.replace(
            r#"[[[1,[0,1]]],[[6,[1,0]]]]"#,
            r#"[[[1,[0,1]]],[[1,[1,0]]]]"#
        )
    );
}

#[test]
fn a_basis_that_is_not_groebner_is_a_typed_defect() {
    // G = [x^2 - y, xy - 1] generates the right ideal and is monic, sorted,
    // and interreduced, but the S-polynomial of the pair (0,1) does not
    // reduce to zero over G.
    let input = list(&["x^2 - y", "x*y - 1"]);
    let basis = input.clone();
    let origins = vec![list(&["1", "0"]), list(&["0", "1"])];

    assert_eq!(
        emit(&input, &basis, &origins),
        Err(CertifyError::Emitter(EmitterFault::NotAGroebnerBasis {
            i: 0,
            j: 1
        }))
    );
    assert_eq!(
        represent::spair_representations(&basis, P, &mut budget()),
        Err(CertifyError::Emitter(EmitterFault::NotAGroebnerBasis {
            i: 0,
            j: 1
        }))
    );
}

#[test]
fn an_input_that_does_not_reduce_is_a_typed_defect() {
    let input = list(&["y"]);
    let basis = list(&["x"]);
    let origins = vec![list(&["0"])];

    assert_eq!(
        emit(&input, &basis, &origins),
        Err(CertifyError::Emitter(EmitterFault::InputHasRemainder {
            input: 0
        }))
    );
    assert_eq!(
        represent::membership_representations(&input, &basis, P, &mut budget()),
        Err(CertifyError::Emitter(EmitterFault::InputHasRemainder {
            input: 0
        }))
    );
}

#[test]
fn a_count_that_does_not_match_is_a_typed_defect() {
    let (input, basis, origins) = pair_system();
    assert_eq!(
        emit(&input, &basis, &origins[..1]),
        Err(CertifyError::Emitter(EmitterFault::OriginCount {
            found: 1,
            expected: 2
        }))
    );

    let short = vec![origins[0][..1].to_vec(), origins[1].clone()];
    assert_eq!(
        emit(&input, &basis, &short),
        Err(CertifyError::Emitter(EmitterFault::OriginEntryCount {
            basis: 0,
            found: 1,
            expected: 2
        }))
    );

    let zero_basis = vec![ring().zero()];
    assert_eq!(
        emit(&input, &zero_basis, &[Vec::new()]),
        Err(CertifyError::Emitter(EmitterFault::BasisElementZero {
            index: 0
        }))
    );

    // A width the certificate does not declare stops the emitter before it
    // compares monomials of different widths.
    let wide_ring = PolynomialRing::prime_field(P, ["x", "y", "z"]).expect("7 is prime");
    let wide = vec![wide_ring.polynomial([(1, [1, 0, 0])]).expect("fits")];
    assert_eq!(
        emit(&wide, &basis, &origins),
        Err(CertifyError::Emitter(EmitterFault::ExponentCount {
            at: Place::Input(0),
            found: 3,
            expected: 2
        }))
    );
}

#[test]
fn division_records_quotients_that_rebuild_the_dividend() {
    let divisors = list(&["x^2 - y", "x*y - 1"]);
    for text in ["x^3 + 2*y", "3*x*y^2 + x*y + 5", "x^2*y^2", "0"] {
        let f = one(text);
        let (quotients, remainder) =
            divide(&f, &divisors, P, &budget()).expect("the budget stops nothing");
        assert_eq!(quotients.len(), divisors.len());

        let mut sum = remainder.clone();
        for (quotient, divisor) in quotients.iter().zip(&divisors) {
            sum = sum.add(&product(quotient, divisor), P);
        }
        assert_eq!(sum, f, "f = sum q_j g_j + r");

        for (quotient, divisor) in quotients.iter().zip(&divisors) {
            if quotient.is_zero() {
                continue;
            }
            let lead = quotient
                .lm()
                .expect("a nonzero quotient has a leading monomial")
                .checked_mul(divisor.lm().expect("a divisor is nonzero"))
                .expect("the product fits");
            assert!(
                Some(&lead) <= f.lm(),
                "lm(q_j g_j) must stay at or below lm(f)"
            );
        }
    }
}

/// An instant in the past, so every deadline check stops the work.
fn passed() -> Option<Instant> {
    Some(Instant::now() - Duration::from_secs(1))
}

#[test]
fn a_passed_deadline_stops_the_emitter() {
    let (input, basis, origins) = base_system();
    assert_eq!(
        assemble(
            &ring(),
            &input,
            basis.clone(),
            origins.clone(),
            &mut Budget::new(passed(), None, 2),
        ),
        Err(CertifyError::Engine(ComputeError::Timeout))
    );

    let divisors = list(&["x^2 - y", "x*y - 1"]);
    assert_eq!(
        divide(
            &one("x^3 + 2*y"),
            &divisors,
            P,
            &Budget::new(passed(), None, 2)
        ),
        Err(ComputeError::Timeout)
    );
    assert_eq!(
        represent::spair_representations(&basis, P, &mut Budget::new(passed(), None, 2)),
        Err(CertifyError::Engine(ComputeError::Timeout))
    );
}

#[test]
fn a_tight_memory_limit_stops_the_emitter() {
    let (input, basis, origins) = base_system();
    for limit in [1, 64, 256] {
        assert_eq!(
            assemble(
                &ring(),
                &input,
                basis.clone(),
                origins.clone(),
                &mut Budget::new(None, Some(limit), 2),
            ),
            Err(CertifyError::Engine(ComputeError::MemoryLimitExceeded)),
            "a limit of {limit} bytes cannot hold the emission"
        );
    }

    let divisors = list(&["x^2 - y", "x*y - 1"]);
    assert_eq!(
        divide(
            &one("x^3 + 2*y"),
            &divisors,
            P,
            &Budget::new(None, Some(1), 2)
        ),
        Err(ComputeError::MemoryLimitExceeded)
    );
}

#[test]
fn a_limit_that_holds_the_emission_changes_no_byte() {
    let (input, basis, origins) = base_system();
    let free = emit(&input, &basis, &origins).expect("the basis holds");
    let bounded = assemble(
        &ring(),
        &input,
        basis,
        origins,
        &mut Budget::new(
            Some(Instant::now() + Duration::from_secs(60)),
            Some(1 << 20),
            2,
        ),
    )
    .expect("the basis holds inside the limit");
    assert_eq!(free, bounded);
}
