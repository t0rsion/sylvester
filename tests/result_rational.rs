//! Exact equality checks for saved rational computation records.

use std::time::Duration;

use sylvester::{
    Budget, CancellationToken, ClaimedProvenance, EnvelopeError, EqualityCheckError, ResultEnvelope,
};

fn record(input: &[&str], basis: &[&str], claim: ClaimedProvenance) -> ResultEnvelope {
    let value = serde_json::json!({
        "schema": "sylv-result-v1",
        "order": "grevlex-v1",
        "domain": {"kind": "rationals"},
        "variables": ["x", "y"],
        "input": input,
        "basis": basis,
        "claimed_provenance": claim,
    });
    ResultEnvelope::from_json(
        &serde_json::to_vec(&value).expect("the record serializes"),
        Budget::new(),
    )
    .expect("the record header is valid")
}

#[test]
fn a_saved_rational_record_returns_exact_origins() {
    let record = record(&["x + y", "y"], &["x", "y"], ClaimedProvenance::Unchanged);
    let checked = record
        .check_rational(Budget::new())
        .expect("the saved ideals are equal");
    assert_eq!(checked.input().generators().len(), 2);
    assert_eq!(checked.basis().len(), 2);
    assert_eq!(checked.origins()[0][0].to_string(), "1");
    assert_eq!(checked.origins()[0][1].to_string(), "-1");
}

#[test]
fn an_equals_claim_does_not_accept_the_unit_bogus_basis() {
    let record = record(&["x"], &["1"], ClaimedProvenance::EqualsInput);
    assert_eq!(
        record.check_rational(Budget::new()),
        Err(EnvelopeError::EqualityCheck(
            EqualityCheckError::ForwardMembership { basis: 0 },
        ))
    );
}

#[test]
fn exact_rational_check_rejects_a_prime_record() {
    let value = serde_json::json!({
        "schema": "sylv-result-v1",
        "order": "grevlex-v1",
        "domain": {"kind": "prime_field", "modulus": 7},
        "variables": ["x"],
        "input": ["x"],
        "basis": ["x"],
        "claimed_provenance": "unverified",
    });
    let record = ResultEnvelope::from_json(
        &serde_json::to_vec(&value).expect("the record serializes"),
        Budget::new(),
    )
    .expect("the prime record header is valid");
    assert!(matches!(
        record.check_rational(Budget::new()),
        Err(EnvelopeError::Format(_))
    ));
}

#[test]
fn a_rational_certificate_is_a_format_error() {
    let value = serde_json::json!({
        "schema": "sylv-result-v1",
        "order": "grevlex-v1",
        "domain": {"kind": "rationals"},
        "variables": ["x"],
        "input": ["x"],
        "basis": ["x"],
        "claimed_provenance": "unverified",
        "certificate": [],
    });
    assert!(matches!(
        ResultEnvelope::from_json(
            &serde_json::to_vec(&value).expect("the record serializes"),
            Budget::new(),
        ),
        Err(EnvelopeError::Format(_))
    ));
}

#[test]
fn parsing_and_limits_remain_typed() {
    let malformed = record(&["x@"], &["x"], ClaimedProvenance::EqualsInput);
    assert!(matches!(
        malformed.check_rational(Budget::new()),
        Err(EnvelopeError::Expression(_))
    ));

    let valid = record(&["x"], &["x"], ClaimedProvenance::EqualsInput);
    assert_eq!(
        valid.check_rational(Budget::new().timeout(Duration::ZERO)),
        Err(EnvelopeError::Timeout)
    );
    assert_eq!(
        valid.check_rational(Budget::new().memory_limit(1)),
        Err(EnvelopeError::MemoryLimitExceeded)
    );
}

#[test]
fn a_pre_cancelled_budget_stops_before_record_parsing() {
    let record = record(&["x"], &["x"], ClaimedProvenance::EqualsInput);
    let token = CancellationToken::new();
    token.cancel();
    assert_eq!(
        record.check_rational(Budget::new().cancellation(token)),
        Err(EnvelopeError::Timeout)
    );
}

#[test]
fn a_large_valid_variable_name_is_charged_before_parsing() {
    let variable = format!("x{}", "a".repeat(4095));
    let value = serde_json::json!({
        "schema": "sylv-result-v1",
        "order": "grevlex-v1",
        "domain": {"kind": "rationals"},
        "variables": [&variable],
        "input": [],
        "basis": [],
        "claimed_provenance": "unverified",
    });
    let record = ResultEnvelope::from_json(
        &serde_json::to_vec(&value).expect("the record serializes"),
        Budget::new(),
    )
    .expect("the variable name is valid");
    assert_eq!(
        record.check_rational(Budget::new().memory_limit(6_000)),
        Err(EnvelopeError::MemoryLimitExceeded)
    );
}
