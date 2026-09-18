//! Exact rational ideal equality checks.

use std::time::{Duration, Instant};

use sylvester::{
    Budget, EqualityCheckError, GroebnerBasis, Polynomial, PolynomialRing, RationalOptions,
    Rationals,
};

fn check_texts(
    variables: &[&str],
    generators: &[&str],
    basis: &[&str],
) -> Result<sylvester::RationalEqualityCheck, EqualityCheckError> {
    let ring = PolynomialRing::rationals(variables).expect("the variable names are valid");
    let generators: Vec<_> = generators
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the input parses"))
        .collect();
    let basis: Vec<_> = basis
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the candidate parses"))
        .collect();
    let ideal = ring
        .ideal(generators)
        .expect("the generators share the ring");
    let basis = GroebnerBasis::<Rationals>::from_polynomials(&ring, basis, Budget::new())
        .expect("the candidate is a reduced Gröbner basis");
    ideal.check_basis_equality(&basis, Budget::new())
}

#[test]
fn equal_ideals_return_exact_origins() {
    let result = check_texts(&["x", "y"], &["x + y", "y"], &["x", "y"])
        .expect("the candidate has the same ideal");
    assert_eq!(result.origins().len(), 2);
    assert_eq!(result.origins()[0].len(), 2);
    assert_eq!(result.origins()[0][0].to_string(), "1");
    assert_eq!(result.origins()[0][1].to_string(), "-1");
    assert_eq!(result.input().generators().len(), 2);
    assert_eq!(result.basis().len(), 2);
    assert!(result.metrics().final_working_bytes > 0);
}

#[test]
fn a_proper_ideal_rejects_the_unit_candidate() {
    let error = check_texts(&["x"], &["x"], &["1"]).expect_err("unit is outside (x)");
    assert_eq!(error, EqualityCheckError::ForwardMembership { basis: 0 });
}

#[test]
fn distinct_principal_ideals_are_rejected() {
    let error = check_texts(&["x"], &["x"], &["x + 1"]).expect_err("the ideals differ");
    assert_eq!(
        error,
        EqualityCheckError::ReverseMembership { generator: 0 }
    );
}

#[test]
fn zero_ideal_has_no_forward_identities() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are valid");
    let ideal = ring
        .ideal([ring.zero(), ring.zero()])
        .expect("the zero generators share the ring");
    let basis = GroebnerBasis::<Rationals>::from_polynomials(&ring, Vec::new(), Budget::new())
        .expect("the empty list is the zero basis");
    let result = ideal
        .check_basis_equality(&basis, Budget::new())
        .expect("both ideals are zero");
    assert!(result.origins().is_empty());
}

#[test]
fn unit_ideal_returns_a_unit_origin() {
    let ring = PolynomialRing::rationals(["x"]).expect("the name is valid");
    let ideal = ring
        .ideal([ring.one()])
        .expect("the unit generator shares the ring");
    let basis =
        GroebnerBasis::<Rationals>::from_polynomials(&ring, vec![ring.one()], Budget::new())
            .expect("one is a reduced basis");
    let result = ideal
        .check_basis_equality(&basis, Budget::new())
        .expect("the unit ideals are equal");
    assert_eq!(result.origins().len(), 1);
    assert_eq!(result.origins()[0][0], ring.one());
}

#[test]
fn zero_variable_ring_supports_zero_and_unit_ideals() {
    let ring =
        PolynomialRing::rationals(Vec::<&str>::new()).expect("an empty variable list is valid");
    let zero_ideal = ring
        .ideal([ring.zero()])
        .expect("the generator shares the ring");
    let zero_basis = GroebnerBasis::<Rationals>::from_polynomials(&ring, Vec::new(), Budget::new())
        .expect("the empty list is the zero basis");
    zero_ideal
        .check_basis_equality(&zero_basis, Budget::new())
        .expect("the zero ideals are equal");

    let unit_ideal = ring
        .ideal([ring.one()])
        .expect("the generator shares the ring");
    let unit_basis =
        GroebnerBasis::<Rationals>::from_polynomials(&ring, vec![ring.one()], Budget::new())
            .expect("one is a reduced basis");
    unit_ideal
        .check_basis_equality(&unit_basis, Budget::new())
        .expect("the unit ideals are equal");
}

#[test]
fn timeout_and_memory_limits_are_typed() {
    let ring = PolynomialRing::rationals(["x"]).expect("the name is valid");
    let generator = ring.parse_polynomial("x").expect("the input parses");
    let ideal = ring
        .ideal([generator.clone()])
        .expect("the generator shares the ring");
    let basis = GroebnerBasis::<Rationals>::from_polynomials(&ring, vec![generator], Budget::new())
        .expect("the candidate is a reduced basis");
    assert_eq!(
        ideal.check_basis_equality(&basis, Budget::new().timeout(Duration::ZERO)),
        Err(EqualityCheckError::Timeout)
    );
    assert_eq!(
        ideal.check_basis_equality(&basis, Budget::new().memory_limit(1)),
        Err(EqualityCheckError::MemoryLimitExceeded)
    );
}

#[test]
fn the_rational_engine_basis_passes_the_exact_check() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are valid");
    let ideal = ring
        .ideal([
            ring.parse_polynomial("x^2 - y").expect("the input parses"),
            ring.parse_polynomial("x*y - 1").expect("the input parses"),
        ])
        .expect("the generators share the ring");
    let basis = ideal
        .groebner_basis(RationalOptions::new())
        .expect("the system fits");
    let result = ideal
        .check_basis_equality(&basis, Budget::new())
        .expect("the engine basis has equal ideals");
    assert_eq!(result.origins().len(), basis.len());
}

#[test]
#[ignore = "tracked Buchberger feasibility run"]
fn benchmark_named_rational_systems() {
    let selected = std::env::var("RATIONAL_EQUALITY_ONLY").ok();
    for name in ["cyclic-5", "katsura-5", "noon-4", "eco-8"] {
        if selected.as_deref().is_some_and(|selected| selected != name) {
            continue;
        }
        let (ring, generators) = benchmark_input(name);
        let ideal = ring
            .ideal(generators)
            .expect("the benchmark generators share the ring");
        let start = Instant::now();
        let candidate = ideal
            .groebner_basis(RationalOptions::new())
            .expect("the rational candidate fits");
        let engine = start.elapsed();
        let start = Instant::now();
        let result = ideal.check_basis_equality(
            &candidate,
            Budget::new()
                .timeout(Duration::from_secs(300))
                .memory_limit(8 << 30),
        );
        let elapsed = start.elapsed();
        match result {
            Ok(result) => eprintln!(
                "{name}: engine={engine:?}, check={elapsed:?}, raw_basis={}, raw_terms={}, final_working_bytes={}",
                result.metrics().raw_basis_elements,
                result.metrics().raw_basis_terms,
                result.metrics().final_working_bytes,
            ),
            Err(error) => panic!(
                "{name}: engine={engine:?}, check={elapsed:?}, equality check failed: {error}"
            ),
        }
    }
}

fn benchmark_input(name: &str) -> (PolynomialRing<Rationals>, Vec<Polynomial<Rationals>>) {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("benchmarks/gb-comparison/inputs")
        .join(format!("{name}-q.sylq"));
    let text = std::fs::read_to_string(path).expect("the benchmark input is readable");
    let mut lines = text.lines();
    let nvars: usize = lines
        .next()
        .expect("the benchmark names its variables")
        .parse()
        .expect("the variable count is an integer");
    let variables: Vec<String> = (1..=nvars).map(|index| format!("x{index}")).collect();
    let ring = PolynomialRing::rationals(&variables).expect("the generated names are valid");
    let generators = lines
        .map(|line| {
            let terms = line.split(';').map(|term| {
                let mut fields = term.split(',');
                let coefficient: i64 = fields
                    .next()
                    .expect("a term has a coefficient")
                    .parse()
                    .expect("the coefficient is an integer");
                let exponents: Vec<u16> = fields
                    .map(|field| field.parse().expect("the exponent is an integer"))
                    .collect();
                (coefficient, exponents)
            });
            ring.polynomial(terms)
                .expect("the benchmark term widths match")
        })
        .collect();
    (ring, generators)
}
