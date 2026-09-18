//! Computation records preserve claims without granting verification.

use std::time::Duration;

use sylvester::{
    Budget, ClaimedProvenance, ComputeOptions, EnvelopeError, PolynomialRing, RationalOptions,
    ResultEnvelope,
};

#[test]
fn duplicate_fields_are_rejected_in_records_and_domains() {
    let original = r#"{"schema":"sylv-result-v1","order":"grevlex-v1","domain":{"kind":"prime_field","modulus":7},"variables":["x"],"input":["x"],"basis":["x"],"claimed_provenance":"unverified","certificate":null}"#;
    for field in [r#""basis":["1"]"#, r#""certificate":null"#] {
        let duplicate = format!("{},{}{}", &original[..original.len() - 1], field, '}');
        assert!(matches!(
            ResultEnvelope::from_json(duplicate.as_bytes(), Budget::new()),
            Err(EnvelopeError::Format(_))
        ));
    }
    let duplicate = original.replace(r#""modulus":7"#, r#""modulus":7,"modulus":11"#);
    assert!(matches!(
        ResultEnvelope::from_json(duplicate.as_bytes(), Budget::new()),
        Err(EnvelopeError::Format(_))
    ));
}

#[test]
fn a_prime_record_preserves_input_basis_and_certificate_bytes() {
    let ring = PolynomialRing::prime_field(7, ["x", "y"]).unwrap();
    let input = ring
        .ideal([ring.parse_polynomial("x^2 - y").unwrap()])
        .unwrap();
    let certified = input
        .groebner_basis_certified(ComputeOptions::new())
        .unwrap();
    let record = ResultEnvelope::from_certified(&input, &certified, Budget::new()).unwrap();
    let bytes = record.to_json(Budget::new()).unwrap();
    let loaded = ResultEnvelope::from_json(&bytes, Budget::new()).unwrap();
    assert_eq!(loaded, record);
    assert_eq!(loaded.variables(), ["x", "y"]);
    assert_eq!(loaded.input(), ["x^2 + 6*y"]);
    assert_eq!(loaded.claimed_provenance(), ClaimedProvenance::Certified);
    assert_eq!(loaded.certificate(), Some(certified.certificate()));
    let verified = loaded
        .verify_prime(&sylvester::verify::Limits::default(), Budget::new())
        .unwrap();
    assert_eq!(verified.modulus(), 7);
}

#[test]
fn a_valid_certificate_cannot_authenticate_changed_record_data() {
    let ring = PolynomialRing::prime_field(7, ["x"]).unwrap();
    let input = ring.ideal([ring.generator(0).unwrap()]).unwrap();
    let certified = input
        .groebner_basis_certified(ComputeOptions::new())
        .unwrap();
    let record = ResultEnvelope::from_prime(
        &input,
        certified.basis(),
        Some(certified.certificate()),
        Budget::new(),
    )
    .unwrap();
    let original = record.to_json(Budget::new()).unwrap();
    for field in ["input", "basis"] {
        let mut changed: serde_json::Value = serde_json::from_slice(&original).unwrap();
        changed[field] = serde_json::json!(["1"]);
        let loaded =
            ResultEnvelope::from_json(&serde_json::to_vec(&changed).unwrap(), Budget::new())
                .unwrap();
        assert_eq!(
            loaded.verify_prime(&sylvester::verify::Limits::default(), Budget::new()),
            Err(EnvelopeError::CertificateMismatch),
        );
    }
}

#[test]
fn envelope_verification_keeps_the_callers_caps() {
    let ring = PolynomialRing::prime_field(7, ["x"]).unwrap();
    let input = ring.ideal([ring.generator(0).unwrap()]).unwrap();
    let certified = input
        .groebner_basis_certified(ComputeOptions::new())
        .unwrap();
    let record = ResultEnvelope::from_prime(
        &input,
        certified.basis(),
        Some(certified.certificate()),
        Budget::new(),
    )
    .unwrap();
    let caps = sylvester::verify::Limits {
        max_bytes: 0,
        ..Default::default()
    };
    assert!(matches!(record.verify_prime(&caps, Budget::new()),
        Err(EnvelopeError::Verification(error)) if error.is_exhaustion()));
}

#[test]
fn changing_a_claim_does_not_run_a_verifier() {
    let bytes = br#"{
        "schema":"sylv-result-v1","order":"grevlex-v1",
        "domain":{"kind":"prime_field","modulus":7},
        "variables":["x"],"input":["x"],"basis":["1"],
        "claimed_provenance":"certified","certificate":[1,2,3]
    }"#;
    let loaded = ResultEnvelope::from_json(bytes, Budget::new()).unwrap();
    assert_eq!(loaded.claimed_provenance(), ClaimedProvenance::Certified);
    assert!(sylvester::verify::verify(loaded.certificate().unwrap()).is_err());
    assert_eq!(loaded.basis(), ["1"]);
}

#[test]
fn attaching_bytes_does_not_make_a_certified_claim() {
    let ring = PolynomialRing::prime_field(7, ["x"]).unwrap();
    let input = ring.ideal([ring.one()]).unwrap();
    let basis = input.groebner_basis(ComputeOptions::new()).unwrap();
    let record =
        ResultEnvelope::from_prime(&input, &basis, Some(&[1, 2, 3]), Budget::new()).unwrap();
    assert_eq!(record.claimed_provenance(), ClaimedProvenance::Unverified);
}

#[test]
fn rational_records_keep_the_weaker_claim() {
    let ring = PolynomialRing::rationals(["x"]).unwrap();
    let input = ring
        .ideal([ring.parse_polynomial("x^2 - 1/2").unwrap()])
        .unwrap();
    let basis = input.groebner_basis(RationalOptions::new()).unwrap();
    let record = ResultEnvelope::from_rational(&input, &basis, Budget::new()).unwrap();
    let loaded =
        ResultEnvelope::from_json(&record.to_json(Budget::new()).unwrap(), Budget::new()).unwrap();
    assert_eq!(loaded.claimed_provenance(), ClaimedProvenance::Unchanged);
    assert_eq!(loaded.input(), ["x^2 - 1/2"]);
    assert!(loaded.certificate().is_none());
}

#[test]
fn a_rational_equality_record_keeps_the_checked_statement() {
    let ring = PolynomialRing::rationals(["x"]).unwrap();
    let input = ring
        .ideal([ring.parse_polynomial("x^2 - 1").unwrap()])
        .unwrap();
    let basis = input.groebner_basis(RationalOptions::new()).unwrap();
    let checked = input.check_basis_equality(&basis, Budget::new()).unwrap();
    let record = ResultEnvelope::from_checked_rational(&checked, Budget::new()).unwrap();
    let bytes = record.to_json(Budget::new()).unwrap();
    let loaded = ResultEnvelope::from_json(&bytes, Budget::new()).unwrap();
    assert_eq!(loaded.input(), ["x^2 - 1"]);
    assert_eq!(loaded.basis(), ["x^2 - 1"]);
    assert_eq!(loaded.claimed_provenance(), ClaimedProvenance::EqualsInput);
    assert!(loaded.certificate().is_none());
    assert!(
        loaded
            .verify_prime(&sylvester::verify::Limits::default(), Budget::new())
            .is_err()
    );
}

#[test]
fn malformed_or_wrong_contract_records_are_rejected() {
    for bytes in [b"{}".as_slice(), br#"{"schema":"another"}"#, b"[]"] {
        assert!(matches!(
            ResultEnvelope::from_json(bytes, Budget::new()),
            Err(EnvelopeError::Format(_))
        ));
    }
}

#[test]
fn provenance_shape_matches_the_domain_and_attachment() {
    let mut value = serde_json::json!({
        "schema": "sylv-result-v1", "order": "grevlex-v1",
        "domain": {"kind": "prime_field", "modulus": 7},
        "variables": ["x"], "input": ["x"], "basis": ["x"],
        "claimed_provenance": "certified"
    });
    for claim in ["certified", "unchanged", "contains_input", "equals_input"] {
        value["claimed_provenance"] = claim.into();
        assert!(matches!(
            ResultEnvelope::from_json(&serde_json::to_vec(&value).unwrap(), Budget::new()),
            Err(EnvelopeError::Format(_))
        ));
    }
    value["domain"] = serde_json::json!({"kind": "rationals"});
    value["claimed_provenance"] = "certified".into();
    value["certificate"] = serde_json::json!([1, 2, 3]);
    assert!(matches!(
        ResultEnvelope::from_json(&serde_json::to_vec(&value).unwrap(), Budget::new()),
        Err(EnvelopeError::Format(_))
    ));
}

#[test]
fn decoding_observes_limits_before_allocating() {
    assert_eq!(
        ResultEnvelope::from_json(b"{}", Budget::new().memory_limit(0)),
        Err(EnvelopeError::MemoryLimitExceeded)
    );
    assert_eq!(
        ResultEnvelope::from_json(b"{}", Budget::new().timeout(Duration::ZERO)),
        Err(EnvelopeError::Timeout)
    );
}

#[test]
fn malformed_large_json_reports_timeout_during_decoding() {
    let mut bytes = br#"{"schema":"sylv-result-v1","order":"grevlex-v1","domain":{"kind":"prime_field","modulus":7},"variables":["x"],"input":["x"#.to_vec();
    bytes.extend(std::iter::repeat_n(b' ', 8 << 20));
    assert_eq!(
        ResultEnvelope::from_json(&bytes, Budget::new().timeout(Duration::from_millis(1)),),
        Err(EnvelopeError::Timeout)
    );
}

#[test]
fn records_cannot_mix_input_and_basis_rings() {
    let first = PolynomialRing::prime_field(7, ["x"]).unwrap();
    let second = PolynomialRing::prime_field(5, ["x"]).unwrap();
    let input = first.ideal([first.one()]).unwrap();
    let basis = second
        .ideal([second.one()])
        .unwrap()
        .groebner_basis(ComputeOptions::new())
        .unwrap();
    assert_eq!(
        ResultEnvelope::from_prime(&input, &basis, None, Budget::new()),
        Err(EnvelopeError::RingMismatch)
    );
}
