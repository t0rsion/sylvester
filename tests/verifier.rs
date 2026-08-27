//! End-to-end tests for the certificate verifier.
//!
//! Every certificate here is written by hand as canonical JSON. The valid
//! ones are accepted from raw bytes. The tamper corpus starts from one
//! valid certificate and changes one thing per test.
//!
//! The base system is F = {x^2 - y, xy - 1} over F_7 in two variables. Its
//! reduced grevlex basis is G = [x^2 - y, xy - 1, y^2 - x]. The pairs (0,1)
//! and (1,2) share a variable and carry entries. The pair (0,2) has coprime
//! leading monomials and is left out.

use std::time::{Duration, Instant};

use sylvester::verify::{
    BasisFault, Cap, Limits, Location, PolyFault, Syntax, VerifyError, verify, verify_with_limits,
};

const BASE: &str = concat!(
    r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
    r#""input":[[[1,[2,0]],[6,[0,1]]],[[1,[1,1]],[6,[0,0]]]],"#,
    r#""basis":[[[1,[2,0]],[6,[0,1]]],[[1,[1,1]],[6,[0,0]]],[[1,[0,2]],[6,[1,0]]]],"#,
    r#""origin":[[[[1,[0,0]]],[]],[[],[[1,[0,0]]]],[[[6,[0,1]]],[[1,[1,0]]]]],"#,
    r#""membership":[[[[1,[0,0]]],[],[]],[[],[[1,[0,0]]],[]]],"#,
    r#""spairs":[[0,1,[[],[],[[6,[0,0]]]]],[1,2,[[[1,[0,0]]],[],[]]]]}"#
);

/// A certificate for F = {x}, G = {x} in two variables.
const SINGLE: &str = concat!(
    r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
    r#""input":[[[1,[1,0]]]],"basis":[[[1,[1,0]]]],"#,
    r#""origin":[[[[1,[0,0]]]]],"membership":[[[[1,[0,0]]]]],"spairs":[]}"#
);

/// A certificate for F = {x, x + 1}, G = {1}, with 1 = (x + 1) - x.
const UNIT: &str = concat!(
    r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
    r#""input":[[[1,[1,0]]],[[1,[1,0]],[1,[0,0]]]],"basis":[[[1,[0,0]]]],"#,
    r#""origin":[[[[6,[0,0]]],[[1,[0,0]]]]],"#,
    r#""membership":[[[[1,[1,0]]]],[[[1,[1,0]],[1,[0,0]]]]],"spairs":[]}"#
);

/// A certificate for F = {x^2 - 1, xy - 1}, G = [y^2 - 1, x - y] over F_7.
const PAIR: &str = concat!(
    r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
    r#""input":[[[1,[2,0]],[6,[0,0]]],[[1,[1,1]],[6,[0,0]]]],"#,
    r#""basis":[[[1,[0,2]],[6,[0,0]]],[[1,[1,0]],[6,[0,1]]]],"#,
    r#""origin":[[[[6,[0,2]]],[[1,[1,1]],[1,[0,0]]]],[[[1,[0,1]]],[[6,[1,0]]]]],"#,
    r#""membership":[[[[1,[0,0]]],[[1,[1,0]],[1,[0,1]]]],[[[1,[0,0]]],[[1,[0,1]]]]],"#,
    r#""spairs":[[0,1,[[[1,[0,1]]],[[6,[0,0]]]]]]}"#
);

/// A certificate for F = {0}, G = {}, the zero ideal.
const ZERO: &str = concat!(
    r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
    r#""input":[[]],"basis":[],"origin":[],"membership":[[]],"spairs":[]}"#
);

/// A certificate for F = {1}, G = {1} with no variables.
const NO_VARS: &str = concat!(
    r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":0,"#,
    r#""input":[[[1,[]]]],"basis":[[[1,[]]]],"#,
    r#""origin":[[[[1,[]]]]],"membership":[[[[1,[]]]]],"spairs":[]}"#
);

fn tampered(from: &str, to: &str) -> Vec<u8> {
    assert_eq!(
        BASE.matches(from).count(),
        1,
        "the tamper pattern must occur once: {from}"
    );
    BASE.replace(from, to).into_bytes()
}

fn reject(from: &str, to: &str) -> VerifyError {
    verify(&tampered(from, to)).expect_err("the verifier must reject the change")
}

fn malformed(from: &str, to: &str) -> Syntax {
    match reject(from, to) {
        VerifyError::Malformed { reason, .. } => reason,
        other => panic!("expected a malformed certificate, got {other}"),
    }
}

#[test]
fn accepts_a_single_generator() {
    let verified = verify(SINGLE.as_bytes()).expect("the certificate holds");
    assert_eq!(verified.modulus(), 7);
    assert_eq!(verified.nvars(), 2);
    assert_eq!(verified.input().len(), 1);
    assert_eq!(verified.basis().len(), 1);
    assert_eq!(verified.basis()[0].terms()[0].mono().exps(), &[1, 0]);
}

#[test]
fn accepts_the_unit_ideal() {
    let verified = verify(UNIT.as_bytes()).expect("the certificate holds");
    assert_eq!(verified.basis().len(), 1);
    assert_eq!(verified.basis()[0].terms().len(), 1);
    assert_eq!(verified.basis()[0].terms()[0].coeff(), 1);
    assert_eq!(verified.basis()[0].terms()[0].mono().degree(), 0);
}

#[test]
fn accepts_a_two_element_basis_with_hand_derived_cofactors() {
    let verified = verify(PAIR.as_bytes()).expect("the certificate holds");
    assert_eq!(verified.basis().len(), 2);
    assert_eq!(verified.basis()[0].terms()[0].mono().exps(), &[0, 2]);
    assert_eq!(verified.basis()[1].terms()[0].mono().exps(), &[1, 0]);
    assert_eq!(verified.input().len(), 2);
}

#[test]
fn accepts_the_zero_ideal() {
    let verified = verify(ZERO.as_bytes()).expect("the certificate holds");
    assert!(verified.basis().is_empty());
    assert_eq!(verified.input().len(), 1);
    assert!(verified.input()[0].is_zero());
}

#[test]
fn accepts_a_certificate_without_variables() {
    let verified = verify(NO_VARS.as_bytes()).expect("the certificate holds");
    assert_eq!(verified.nvars(), 0);
    assert_eq!(verified.basis().len(), 1);
}

#[test]
fn accepts_an_omitted_coprime_spair() {
    // The only pair of PAIR has coprime leading monomials, so the entry may
    // go. The base certificate leaves out its coprime pair (0,2) as well.
    let bytes = PAIR.replace(
        r#""spairs":[[0,1,[[[1,[0,1]]],[[6,[0,0]]]]]]"#,
        r#""spairs":[]"#,
    );
    verify(bytes.as_bytes()).expect("the certificate holds");
    verify(BASE.as_bytes()).expect("the certificate holds");
}

#[test]
fn accepts_the_base_certificate() {
    let verified = verify(BASE.as_bytes()).expect("the certificate holds");
    assert_eq!(verified.basis().len(), 3);
    assert_eq!(verified.input().len(), 2);
    let (input, basis) = verified.into_parts();
    assert_eq!(input.len(), 2);
    assert_eq!(basis[2].terms()[0].mono().exps(), &[0, 2]);
}

#[test]
fn rejects_a_unit_basis_for_a_proper_ideal() {
    // The claim 1 = 1 * x is the failure mode an output-only check accepts.
    let bytes = concat!(
        r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
        r#""input":[[[1,[1,0]]]],"basis":[[[1,[0,0]]]],"#,
        r#""origin":[[[[1,[0,0]]]]],"membership":[[[[1,[1,0]]]]],"spairs":[]}"#
    );
    assert_eq!(
        verify(bytes.as_bytes()),
        Err(VerifyError::OriginIdentity { basis: 0 })
    );
}

#[test]
fn rejects_a_basis_that_is_not_groebner() {
    // G = [x^2 - y, xy - 1] generates the right ideal and is monic, sorted,
    // and interreduced, but it is not a Gröbner basis: the pair (0,1) is not
    // coprime and its S-polynomial has no standard representation over G.
    let bytes = concat!(
        r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
        r#""input":[[[1,[2,0]],[6,[0,1]]],[[1,[1,1]],[6,[0,0]]]],"#,
        r#""basis":[[[1,[2,0]],[6,[0,1]]],[[1,[1,1]],[6,[0,0]]]],"#,
        r#""origin":[[[[1,[0,0]]],[]],[[],[[1,[0,0]]]]],"#,
        r#""membership":[[[[1,[0,0]]],[]],[[],[[1,[0,0]]]]],"spairs":[]}"#
    );
    assert_eq!(
        verify(bytes.as_bytes()),
        Err(VerifyError::SpairMissing { i: 0, j: 1 })
    );
}

#[test]
fn rejects_a_wrong_schema_string() {
    assert_eq!(
        reject("sylv-gb-cert-v1", "sylv-gb-cert-v2"),
        VerifyError::Schema {
            found: "sylv-gb-cert-v2".to_string()
        }
    );
}

#[test]
fn rejects_a_wrong_order_string() {
    assert_eq!(
        reject("grevlex-v1", "lex-v1"),
        VerifyError::Order {
            found: "lex-v1".to_string()
        }
    );
}

#[test]
fn rejects_a_composite_modulus() {
    assert_eq!(
        reject(r#""modulus":7"#, r#""modulus":9"#),
        VerifyError::Modulus {
            found: 9,
            composite: true
        }
    );
}

#[test]
fn rejects_a_modulus_out_of_range() {
    assert_eq!(
        reject(r#""modulus":7"#, r#""modulus":2147483649"#),
        VerifyError::Modulus {
            found: 2147483649,
            composite: false
        }
    );
}

#[test]
fn rejects_nvars_out_of_range() {
    assert_eq!(
        reject(r#""nvars":2"#, r#""nvars":257"#),
        VerifyError::Nvars {
            found: 257,
            max: 256
        }
    );
}

#[test]
fn rejects_a_zero_coefficient() {
    assert_eq!(
        reject(r#""input":[[[1,[2,0]]"#, r#""input":[[[0,[2,0]]"#),
        VerifyError::Polynomial {
            at: Location::Input(0),
            fault: PolyFault::CoefficientZero { term: 0 }
        }
    );
}

#[test]
fn rejects_a_coefficient_at_the_modulus() {
    assert_eq!(
        reject(r#""input":[[[1,[2,0]]"#, r#""input":[[[7,[2,0]]"#),
        VerifyError::Polynomial {
            at: Location::Input(0),
            fault: PolyFault::CoefficientOutOfRange { term: 0, coeff: 7 }
        }
    );
}

#[test]
fn rejects_unsorted_terms() {
    assert_eq!(
        reject(
            r#""input":[[[1,[2,0]],[6,[0,1]]]"#,
            r#""input":[[[6,[0,1]],[1,[2,0]]]"#
        ),
        VerifyError::Polynomial {
            at: Location::Input(0),
            fault: PolyFault::NotDescending { term: 1 }
        }
    );
}

#[test]
fn rejects_a_duplicate_monomial() {
    assert_eq!(
        reject(
            r#""input":[[[1,[2,0]],[6,[0,1]]]"#,
            r#""input":[[[1,[2,0]],[6,[2,0]]]"#
        ),
        VerifyError::Polynomial {
            at: Location::Input(0),
            fault: PolyFault::DuplicateMonomial { term: 1 }
        }
    );
}

#[test]
fn rejects_an_exponent_vector_of_the_wrong_length() {
    assert_eq!(
        reject(r#""input":[[[1,[2,0]]"#, r#""input":[[[1,[2,0,0]]"#),
        VerifyError::Polynomial {
            at: Location::Input(0),
            fault: PolyFault::ExponentCount {
                term: 0,
                found: 3,
                expected: 2
            }
        }
    );
}

#[test]
fn rejects_an_exponent_array_above_the_contract_width() {
    // The array holds far more than 256 entries, well under the default
    // byte cap, so the width check is what stops it, not `max_bytes`.
    let huge: String = std::iter::repeat_n("0,", 10_000).collect::<String>() + "0";
    let from = r#""input":[[[1,[2,0]]"#;
    let to = format!(r#""input":[[[1,[{huge}]]"#);
    assert!(tampered(from, &to).len() < Limits::default().max_bytes);
    assert_eq!(malformed(from, &to), Syntax::ExponentWidth { max: 256 });
}

#[test]
fn rejects_a_non_monic_basis_element() {
    assert_eq!(
        reject(r#""basis":[[[1,[2,0]]"#, r#""basis":[[[2,[2,0]]"#),
        VerifyError::Basis {
            index: 0,
            fault: BasisFault::NotMonic { lc: 2 }
        }
    );
}

#[test]
fn rejects_an_unsorted_basis() {
    assert_eq!(
        reject(
            r#""basis":[[[1,[2,0]],[6,[0,1]]],[[1,[1,1]],[6,[0,0]]],[[1,[0,2]],[6,[1,0]]]]"#,
            r#""basis":[[[1,[2,0]],[6,[0,1]]],[[1,[0,2]],[6,[1,0]]],[[1,[1,1]],[6,[0,0]]]]"#
        ),
        VerifyError::Basis {
            index: 2,
            fault: BasisFault::NotDescending
        }
    );
}

#[test]
fn rejects_a_basis_that_is_not_interreduced() {
    // The tail of element 1 becomes y^2, which the leading monomial of
    // element 2 divides.
    assert_eq!(
        reject(
            r#",[[1,[1,1]],[6,[0,0]]],[[1,[0,2]],[6,[1,0]]]]"#,
            r#",[[1,[1,1]],[6,[0,2]]],[[1,[0,2]],[6,[1,0]]]]"#
        ),
        VerifyError::Basis {
            index: 1,
            fault: BasisFault::TailReducible { term: 1, by: 2 }
        }
    );
}

#[test]
fn rejects_a_redundant_basis_element() {
    // Element 0 becomes x^2*y - y, and the leading monomial xy of element 1
    // divides x^2*y.
    assert_eq!(
        reject(r#""basis":[[[1,[2,0]]"#, r#""basis":[[[1,[2,1]]"#),
        VerifyError::Basis {
            index: 0,
            fault: BasisFault::LeadDivisible { by: 1 }
        }
    );
}

#[test]
fn rejects_a_wrong_origin_cofactor() {
    assert_eq!(
        reject(
            r#"[[[6,[0,1]]],[[1,[1,0]]]]"#,
            r#"[[[1,[0,1]]],[[1,[1,0]]]]"#
        ),
        VerifyError::OriginIdentity { basis: 2 }
    );
}

#[test]
fn rejects_a_wrong_membership_cofactor() {
    assert_eq!(
        reject(
            r#""membership":[[[[1,[0,0]]],[],[]],[[],[[1,[0,0]]],[]]]"#,
            r#""membership":[[[[1,[0,0]]],[],[]],[[],[[2,[0,0]]],[]]]"#
        ),
        VerifyError::MembershipIdentity { input: 1 }
    );
}

#[test]
fn rejects_a_dropped_spair_for_a_non_coprime_pair() {
    assert_eq!(
        reject(r#",[1,2,[[[1,[0,0]]],[],[]]]"#, ""),
        VerifyError::SpairMissing { i: 1, j: 2 }
    );
}

#[test]
fn rejects_a_wrong_spair_cofactor() {
    assert_eq!(
        reject(
            r#"[1,2,[[[1,[0,0]]],[],[]]]"#,
            r#"[1,2,[[[2,[0,0]]],[],[]]]"#
        ),
        VerifyError::SpairIdentity { i: 1, j: 2 }
    );
}

#[test]
fn rejects_a_summand_above_the_spair_bound() {
    // The sum stays correct: the added multiples of g_0 and g_1 cancel each
    // other, and each one leads with x^3*y, far above the S-polynomial.
    assert_eq!(
        reject(
            r#"[0,1,[[],[],[[6,[0,0]]]]]"#,
            r#"[0,1,[[[1,[1,1]],[6,[0,0]]],[[6,[2,0]],[1,[0,1]]],[[6,[0,0]]]]]"#
        ),
        VerifyError::SpairBound {
            i: 0,
            j: 1,
            summand: 0
        }
    );
}

#[test]
fn rejects_a_duplicate_spair_entry() {
    assert_eq!(
        reject(
            r#""spairs":[[0,1,"#,
            r#""spairs":[[0,1,[[],[],[[6,[0,0]]]]],[0,1,"#
        ),
        VerifyError::SpairDuplicate { i: 0, j: 1 }
    );
}

#[test]
fn rejects_unsorted_spair_entries() {
    assert_eq!(
        reject(
            r#""spairs":[[0,1,[[],[],[[6,[0,0]]]]],[1,2,[[[1,[0,0]]],[],[]]]]"#,
            r#""spairs":[[1,2,[[[1,[0,0]]],[],[]]],[0,1,[[],[],[[6,[0,0]]]]]]"#
        ),
        VerifyError::SpairUnsorted { i: 0, j: 1 }
    );
}

#[test]
fn rejects_a_spair_index_out_of_range() {
    assert_eq!(
        reject(
            r#"[1,2,[[[1,[0,0]]],[],[]]]"#,
            r#"[1,3,[[[1,[0,0]]],[],[]]]"#
        ),
        VerifyError::IndexOutOfRange {
            what: "the second S-pair index",
            index: 3,
            bound: 3
        }
    );
}

#[test]
fn rejects_spair_indices_that_do_not_increase() {
    assert_eq!(
        reject(
            r#"[1,2,[[[1,[0,0]]],[],[]]]"#,
            r#"[2,1,[[[1,[0,0]]],[],[]]]"#
        ),
        VerifyError::SpairIndexOrder { i: 2, j: 1 }
    );
}

#[test]
fn rejects_a_spair_entry_with_the_wrong_cofactor_count() {
    assert_eq!(
        reject(r#"[1,2,[[[1,[0,0]]],[],[]]]"#, r#"[1,2,[[[1,[0,0]]],[]]]"#),
        VerifyError::CountMismatch {
            what: "S-pair entry",
            index: Some(1),
            found: 2,
            expected: 3
        }
    );
}

#[test]
fn rejects_an_origin_count_below_the_basis_count() {
    assert_eq!(
        reject(r#",[[[6,[0,1]]],[[1,[1,0]]]]]"#, "]"),
        VerifyError::CountMismatch {
            what: "origin",
            index: None,
            found: 2,
            expected: 3
        }
    );
}

#[test]
fn rejects_an_origin_entry_that_does_not_match_the_input_count() {
    assert_eq!(
        reject(
            r#""origin":[[[[1,[0,0]]],[]],"#,
            r#""origin":[[[[1,[0,0]]],[],[]],"#
        ),
        VerifyError::CountMismatch {
            what: "origin entry",
            index: Some(0),
            found: 3,
            expected: 2
        }
    );
}

#[test]
fn rejects_a_membership_entry_that_does_not_match_the_basis_count() {
    assert_eq!(
        reject(
            r#""membership":[[[[1,[0,0]]],[],[]],"#,
            r#""membership":[[[[1,[0,0]]],[]],"#
        ),
        VerifyError::CountMismatch {
            what: "membership entry",
            index: Some(0),
            found: 2,
            expected: 3
        }
    );
}

#[test]
fn rejects_a_duplicate_key() {
    assert_eq!(
        malformed(
            r#"{"schema":"sylv-gb-cert-v1","#,
            r#"{"schema":"sylv-gb-cert-v1","schema":"sylv-gb-cert-v1","#
        ),
        Syntax::DuplicateKey("schema".to_string())
    );
}

#[test]
fn rejects_an_unknown_key() {
    assert_eq!(
        malformed(r#""input":"#, r#""extra":1,"input":"#),
        Syntax::UnknownKey("extra".to_string())
    );
}

#[test]
fn rejects_a_missing_key() {
    assert_eq!(
        malformed(
            r#","spairs":[[0,1,[[],[],[[6,[0,0]]]]],[1,2,[[[1,[0,0]]],[],[]]]]}"#,
            "}"
        ),
        Syntax::MissingKey("spairs")
    );
}

#[test]
fn rejects_keys_out_of_order() {
    assert_eq!(
        malformed(
            r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","#,
            r#"{"order":"grevlex-v1","schema":"sylv-gb-cert-v1","#
        ),
        Syntax::KeyOutOfOrder {
            found: "order".to_string(),
            expected: "schema"
        }
    );
}

#[test]
fn rejects_inserted_whitespace() {
    assert_eq!(
        malformed(r#"{"schema""#, r#"{ "schema""#),
        Syntax::Whitespace
    );
    assert_eq!(
        malformed(r#""nvars":2,"#, "\"nvars\":2,\n"),
        Syntax::Whitespace
    );
}

#[test]
fn rejects_a_leading_zero() {
    assert_eq!(
        malformed(r#""modulus":7"#, r#""modulus":07"#),
        Syntax::LeadingZero
    );
}

#[test]
fn rejects_a_float() {
    assert_eq!(
        malformed(r#""modulus":7"#, r#""modulus":7.0"#),
        Syntax::Float
    );
}

#[test]
fn rejects_an_exponent_part() {
    assert_eq!(
        malformed(r#""modulus":7"#, r#""modulus":7e1"#),
        Syntax::Exponent
    );
}

#[test]
fn rejects_a_signed_integer() {
    assert_eq!(malformed(r#""modulus":7"#, r#""modulus":-7"#), Syntax::Sign);
}

#[test]
fn rejects_trailing_bytes() {
    let mut bytes = BASE.as_bytes().to_vec();
    bytes.push(b' ');
    match verify(&bytes).expect_err("the verifier must reject trailing bytes") {
        VerifyError::Malformed { reason, offset } => {
            assert_eq!(reason, Syntax::TrailingBytes);
            assert_eq!(offset, BASE.len());
        }
        other => panic!("expected a malformed certificate, got {other}"),
    }
}

#[test]
fn rejects_truncated_bytes() {
    let bytes = &BASE.as_bytes()[..BASE.len() - 1];
    match verify(bytes).expect_err("the verifier must reject truncated bytes") {
        VerifyError::Malformed { reason, .. } => assert_eq!(reason, Syntax::Truncated),
        other => panic!("expected a malformed certificate, got {other}"),
    }
}

#[test]
fn rejects_bytes_above_the_byte_cap() {
    let limits = Limits {
        max_bytes: 32,
        ..Limits::default()
    };
    let error = verify_with_limits(BASE.as_bytes(), &limits).expect_err("the cap holds");
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::Bytes,
            limit: 32
        }
    );
    assert!(error.is_exhaustion());
}

#[test]
fn rejects_a_polynomial_count_above_the_cap() {
    let limits = Limits {
        max_polys: 3,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(BASE.as_bytes(), &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::Polynomials,
            limit: 3
        })
    );
}

#[test]
fn rejects_a_term_count_above_the_per_polynomial_cap() {
    let limits = Limits {
        max_terms_per_poly: 1,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(BASE.as_bytes(), &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::TermsPerPolynomial,
            limit: 1
        })
    );
}

#[test]
fn rejects_a_term_count_above_the_total_cap() {
    let limits = Limits {
        max_total_terms: 3,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(BASE.as_bytes(), &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::TotalTerms,
            limit: 3
        })
    );
}

#[test]
fn rejects_work_after_the_deadline() {
    let limits = Limits {
        deadline: Some(Instant::now() - Duration::from_secs(1)),
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(BASE.as_bytes(), &limits),
        Err(VerifyError::DeadlineExceeded)
    );
}

#[test]
fn stops_an_oversized_polynomial_count_before_it_allocates() {
    // The decoder stops at the cap, so it never reads the rest of the array.
    let oversized_polys = 200_000;
    let mut bytes = String::with_capacity(oversized_polys * 16);
    bytes.push_str(
        r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"input":["#,
    );
    for index in 0..oversized_polys {
        if index > 0 {
            bytes.push(',');
        }
        bytes.push_str("[[1,[1,0]]]");
    }
    bytes.push_str(r#"],"basis":[],"origin":[],"membership":[],"spairs":[]}"#);
    assert!(bytes.len() > 2_000_000);

    let limits = Limits {
        max_polys: 8,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(bytes.as_bytes(), &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::Polynomials,
            limit: 8
        })
    );

    let limits = Limits {
        max_bytes: 1024,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(bytes.as_bytes(), &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::Bytes,
            limit: 1024
        })
    );
}

#[test]
fn accepts_a_present_entry_for_a_coprime_pair() {
    // The pair (0,2) may be left out. An entry for it is still checked, and
    // this one holds: S(g_0, g_2) = x*g_0 - y*g_2.
    let bytes = BASE.replace(
        r#""spairs":[[0,1,[[],[],[[6,[0,0]]]]],"#,
        r#""spairs":[[0,1,[[],[],[[6,[0,0]]]]],[0,2,[[[1,[1,0]]],[],[[6,[0,1]]]]],"#,
    );
    let verified = verify(bytes.as_bytes()).expect("the certificate holds");
    assert_eq!(verified.basis().len(), 3);
}

#[test]
fn accepts_a_certificate_over_the_smallest_field() {
    let bytes = concat!(
        r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":2,"nvars":2,"#,
        r#""input":[[[1,[1,0]]]],"basis":[[[1,[1,0]]]],"#,
        r#""origin":[[[[1,[0,0]]]]],"membership":[[[[1,[0,0]]]]],"spairs":[]}"#
    );
    let verified = verify(bytes.as_bytes()).expect("the certificate holds");
    assert_eq!(verified.modulus(), 2);
}

#[test]
fn rejects_a_wrong_cofactor_in_an_entry_for_a_coprime_pair() {
    let bytes = BASE.replace(
        r#""spairs":[[0,1,[[],[],[[6,[0,0]]]]],"#,
        r#""spairs":[[0,1,[[],[],[[6,[0,0]]]]],[0,2,[[],[],[]]],"#,
    );
    assert_eq!(
        verify(bytes.as_bytes()),
        Err(VerifyError::SpairIdentity { i: 0, j: 2 })
    );
}

#[test]
fn rejects_a_zero_basis_element() {
    assert_eq!(
        reject(r#",[[1,[0,2]],[6,[1,0]]]],"origin""#, r#",[]],"origin""#),
        VerifyError::Basis {
            index: 2,
            fault: BasisFault::Zero
        }
    );
}

#[test]
fn rejects_a_repeated_basis_element() {
    // A repeated element shares its leading monomial with the element
    // before it, so the strict sort rejects it.
    assert_eq!(
        reject(
            r#",[[1,[0,2]],[6,[1,0]]]],"origin""#,
            r#",[[1,[1,1]],[6,[0,0]]]],"origin""#
        ),
        VerifyError::Basis {
            index: 2,
            fault: BasisFault::NotDescending
        }
    );
}

#[test]
fn rejects_a_membership_count_below_the_input_count() {
    assert_eq!(
        reject(
            r#""membership":[[[[1,[0,0]]],[],[]],[[],[[1,[0,0]]],[]]]"#,
            r#""membership":[[[[1,[0,0]]],[],[]]]"#
        ),
        VerifyError::CountMismatch {
            what: "membership",
            index: None,
            found: 1,
            expected: 2
        }
    );
}

#[test]
fn rejects_a_modulus_below_two() {
    assert_eq!(
        reject(r#""modulus":7"#, r#""modulus":0"#),
        VerifyError::Modulus {
            found: 0,
            composite: false
        }
    );
    assert_eq!(
        reject(r#""modulus":7"#, r#""modulus":1"#),
        VerifyError::Modulus {
            found: 1,
            composite: false
        }
    );
}

#[test]
fn rejects_a_composite_modulus_that_passes_the_first_two_witnesses() {
    // 1373653 is the smallest composite that passes Miller-Rabin for the
    // witnesses 2 and 3. The verifier uses twelve witnesses, so it names
    // it composite.
    assert_eq!(
        reject(r#""modulus":7"#, r#""modulus":1373653"#),
        VerifyError::Modulus {
            found: 1373653,
            composite: true
        }
    );
}

#[test]
fn rejects_an_exponent_above_the_contract_range() {
    assert_eq!(
        malformed(r#""input":[[[1,[2,0]]"#, r#""input":[[[1,[65536,0]]"#),
        Syntax::IntegerOutOfRange {
            what: "an exponent",
            value: 65536,
            max: 65535
        }
    );
}

#[test]
fn rejects_a_coefficient_above_the_field() {
    assert_eq!(
        reject(r#""input":[[[1,[2,0]]"#, r#""input":[[[4294967296,[2,0]]"#),
        VerifyError::Polynomial {
            at: Location::Input(0),
            fault: PolyFault::CoefficientOutOfRange {
                term: 0,
                coeff: 4294967296
            }
        }
    );
}

#[test]
fn rejects_a_negative_coefficient() {
    assert_eq!(
        malformed(r#""input":[[[1,[2,0]]"#, r#""input":[[[-1,[2,0]]"#),
        Syntax::Sign
    );
}

#[test]
fn rejects_a_negative_exponent() {
    assert_eq!(
        malformed(r#""input":[[[1,[2,0]]"#, r#""input":[[[1,[-2,0]]"#),
        Syntax::Sign
    );
}

#[test]
fn rejects_a_coefficient_above_64_bits() {
    assert_eq!(
        malformed(
            r#""input":[[[1,[2,0]]"#,
            r#""input":[[[99999999999999999999999,[2,0]]"#
        ),
        Syntax::IntegerOverflow {
            what: "a coefficient"
        }
    );
}

#[test]
fn rejects_a_byte_outside_the_canonical_encoding() {
    let inside = BASE.find("grevlex-v1").expect("the order string is there") + 2;
    let mut bytes = BASE.as_bytes().to_vec();
    bytes[inside] = 0xff;
    match verify(&bytes).expect_err("the verifier must reject the byte") {
        VerifyError::Malformed { reason, .. } => assert_eq!(reason, Syntax::BadString),
        other => panic!("expected a malformed certificate, got {other}"),
    }

    let outside = BASE.find(r#""modulus""#).expect("the modulus key is there");
    let mut bytes = BASE.as_bytes().to_vec();
    bytes[outside] = 0xff;
    match verify(&bytes).expect_err("the verifier must reject the byte") {
        VerifyError::Malformed { reason, .. } => assert_eq!(
            reason,
            Syntax::UnexpectedByte {
                found: 0xff,
                expected: "a string"
            }
        ),
        other => panic!("expected a malformed certificate, got {other}"),
    }
}

/// A polynomial of `count` terms of degree `count - 1` in two variables,
/// written as canonical JSON.
fn dense_poly(count: usize) -> String {
    let mut text = String::from("[");
    for index in (0..count).rev() {
        if index + 1 < count {
            text.push(',');
        }
        text.push_str(&format!("[1,[{index},{}]]", count - 1 - index));
    }
    text.push(']');
    text
}

#[test]
fn rejects_a_product_above_the_intermediate_byte_cap() {
    let limits = Limits {
        max_intermediate_bytes: 1,
        ..Limits::default()
    };
    let error = verify_with_limits(BASE.as_bytes(), &limits).expect_err("the cap holds");
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit: 1
        }
    );
    assert!(error.is_exhaustion());
}

#[test]
fn stops_an_oversized_product_under_the_default_caps() {
    // The cofactor and the input polynomial each hold 20000 terms, so the
    // product buffer would hold 20000^2 = 400 million terms. Every decoder
    // cap holds. The product cap stops the work before it allocates.
    let wide = dense_poly(20_000);
    let bytes = format!(
        concat!(
            r#"{{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
            r#""input":[{wide}],"basis":[[[1,[1,0]]]],"origin":[[{wide}]],"#,
            r#""membership":[[[]]],"spairs":[]}}"#
        ),
        wide = wide
    );
    let limits = Limits::default();
    assert!(bytes.len() < limits.max_bytes);
    let error = verify(bytes.as_bytes()).expect_err("the cap holds");
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit: limits.max_intermediate_bytes
        }
    );
    assert!(error.is_exhaustion());
}

/// A polynomial of `count` terms over 256 variables, written as canonical
/// JSON.
///
/// The terms are the powers of the first variable, in descending degree, so
/// the encoding is canonical and the widest the contract allows.
fn wide_poly(count: usize) -> String {
    let zeros = ",0".repeat(255);
    let mut text = String::from("[");
    for degree in (1..=count).rev() {
        if degree < count {
            text.push(',');
        }
        text.push_str(&format!("[1,[{degree}{zeros}]]"));
    }
    text.push(']');
    text
}

#[test]
fn stops_a_product_that_just_crosses_the_intermediate_byte_cap() {
    // One term over 256 variables costs 1064 bytes, and the cap counts
    // bytes. The cofactor and the input polynomial each hold 711 terms, so
    // the product buffer would hold 711^2 = 505,521 terms, just past what
    // the default cap allows at this width. Every other cap holds: the
    // certificate is under a megabyte and holds 1423 terms. Under a cap of
    // eight million terms this product ran to the end.
    let limits = Limits::default();
    let bytes = wide_origin(711);
    assert!(bytes.len() < limits.max_bytes);

    let start = Instant::now();
    let error = verify(bytes.as_bytes()).expect_err("the cap holds");
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit: limits.max_intermediate_bytes
        }
    );
    assert!(error.is_exhaustion());
    // The cap stops the work before the first term of the product is
    // written, so the rejection costs the decoded certificate and nothing
    // else.
    assert!(
        start.elapsed() < Duration::from_secs(10),
        "{:?}",
        start.elapsed()
    );
}

/// A certificate over 256 variables whose origin cofactor and input
/// polynomial each hold `terms` terms. The origin identity is false.
fn wide_origin(terms: usize) -> String {
    format!(
        concat!(
            r#"{{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":32003,"nvars":256,"#,
            r#""input":[{width}],"basis":[[[1,[1{zeros}]]]],"origin":[[{width}]],"#,
            r#""membership":[[[]]],"spairs":[]}}"#
        ),
        width = wide_poly(terms),
        zeros = ",0".repeat(255),
    )
}

#[test]
fn the_intermediate_byte_cap_admits_the_product_that_exactly_fills_it() {
    // The product buffer holds 20 * 20 = 400 terms over 256 variables. A
    // cap of exactly that many terms' worth of bytes admits it, and the
    // check runs to its answer: the identity is false. One byte less stops
    // the work before the buffer is reserved.
    let bytes = wide_origin(20);
    let term_bytes =
        size_of::<sylvester::verify::Term>() + 256 * size_of::<sylvester::verify::Exp>();
    let fits = Limits {
        max_intermediate_bytes: 400 * term_bytes,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(bytes.as_bytes(), &fits),
        Err(VerifyError::OriginIdentity { basis: 0 })
    );

    let tight = Limits {
        max_intermediate_bytes: 400 * term_bytes - 1,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(bytes.as_bytes(), &tight),
        Err(VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit: 400 * term_bytes - 1
        })
    );
}

#[test]
fn a_deadline_outranks_an_invalid_certificate() {
    // The origin identity of this certificate is false. A verifier that ran
    // past its deadline must say so, not hand back a finding the caller did
    // not pay for.
    let bad = BASE.replace(
        r#""origin":[[[[1,[0,0]]],[]]"#,
        r#""origin":[[[[2,[0,0]]],[]]"#,
    );
    assert_eq!(
        verify(bad.as_bytes()),
        Err(VerifyError::OriginIdentity { basis: 0 })
    );
    let limits = Limits {
        deadline: Some(Instant::now() - Duration::from_secs(1)),
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(bad.as_bytes(), &limits),
        Err(VerifyError::DeadlineExceeded)
    );
}

#[test]
fn a_deadline_stops_one_long_multiplication() {
    // The product holds 1200^2 = 1.44 million terms, inside every cap. The
    // deadline passes while the multiplication runs, and the poll inside it
    // returns before the identity does. The identity is false, so a
    // verifier that polled only between entries would answer
    // `OriginIdentity` after finishing the whole product.
    let wide = dense_poly(1200);
    let bytes = format!(
        concat!(
            r#"{{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
            r#""input":[{wide}],"basis":[[[1,[1,0]]]],"origin":[[{wide}]],"#,
            r#""membership":[[[]]],"spairs":[]}}"#
        ),
        wide = wide
    );
    let limits = Limits {
        deadline: Some(Instant::now() + Duration::from_millis(5)),
        ..Limits::default()
    };
    let start = Instant::now();
    assert_eq!(
        verify_with_limits(bytes.as_bytes(), &limits),
        Err(VerifyError::DeadlineExceeded)
    );
    assert!(
        start.elapsed() < Duration::from_secs(5),
        "{:?}",
        start.elapsed()
    );
}

#[test]
fn rejects_an_entry_count_above_the_cap() {
    let limits = Limits {
        max_entries: 3,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(BASE.as_bytes(), &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::Entries,
            limit: 3
        })
    );
}

#[test]
fn stops_an_oversized_entry_count_before_it_allocates() {
    // An entry with no polynomial passes every polynomial cap. The entry
    // cap stops it.
    let entries = 200_000;
    let mut bytes = String::with_capacity(entries * 10);
    bytes.push_str(concat!(
        r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
        r#""input":[],"basis":[],"origin":[],"membership":[],"spairs":["#
    ));
    for index in 0..entries {
        if index > 0 {
            bytes.push(',');
        }
        bytes.push_str("[0,1,[]]");
    }
    bytes.push_str("]}");

    let limits = Limits {
        max_entries: 8,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(bytes.as_bytes(), &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::Entries,
            limit: 8
        })
    );
}

#[test]
fn survives_every_byte_mutation_without_panicking() {
    // Deterministic mutations of accepted certificates: replace one byte,
    // drop one byte, repeat one byte, at every position. The verifier must
    // return a value for each one.
    const REPLACEMENTS: [u8; 12] = [
        0x00, 0x01, 0x20, b'"', b',', b'0', b'9', b'[', b']', b'{', 0x7f, 0xff,
    ];
    let mut mutations = 0;
    for certificate in [SINGLE, UNIT, PAIR, ZERO, NO_VARS, BASE] {
        let original = certificate.as_bytes();
        for position in 0..original.len() {
            for replacement in REPLACEMENTS {
                let mut bytes = original.to_vec();
                bytes[position] = replacement;
                let _ = verify(&bytes);
            }
            let mut dropped = original.to_vec();
            dropped.remove(position);
            let _ = verify(&dropped);
            let mut repeated = original.to_vec();
            repeated.insert(position, original[position]);
            let _ = verify(&repeated);
            mutations += REPLACEMENTS.len() + 2;
        }
    }
    // The corpus is fixed, so the count is stable; the floor guards against
    // a future edit that silently shrinks it.
    assert!(mutations > 10_000, "the corpus ran {mutations} mutations");
}

#[test]
fn survives_every_prefix_without_panicking() {
    for certificate in [SINGLE, UNIT, PAIR, ZERO, NO_VARS, BASE] {
        let original = certificate.as_bytes();
        for end in 0..original.len() {
            let error = verify(&original[..end]).expect_err("a prefix is not a certificate");
            assert!(!error.is_exhaustion());
        }
    }
}
