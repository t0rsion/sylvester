//! The multimodular rational engine.
//!
//! The reference bases under `tests/fixtures/rational` come from Singular
//! 4.4.1, one file per system, generated with
//!
//! ```text
//! option(redSB);
//! ring r = 0, (x1,...,xn), dp;
//! short = 0;
//! ideal I = <generators>;
//! ideal G = simplify(std(I), 1);
//! int k;
//! for (k = size(G); k >= 1; k--) { print(G[k]); }
//! ```
//!
//! `simplify(G, 1)` makes each element monic, and the loop prints the
//! largest leading monomial first, which is the order this crate returns.
//! The generators are the ones `benchmarks/gb-comparison/gen.py` defines
//! for the same families. The comparison is over the full basis,
//! coefficients included.

use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::time::Duration;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use proptest::prelude::*;
use sylvester::{
    ComputeError, ComputeOptions, Established, GroebnerBasis, PolynomialRing, RationalOptions,
    RationalStop, Rationals,
};

/// The largest prime of the sequence, which every run takes first.
const FIRST_PRIME: i64 = (1 << 31) - 1;

/// One reference system: the variables, the generators, and the basis.
struct Reference {
    variables: Vec<String>,
    generators: Vec<String>,
    basis: Vec<String>,
}

/// Read one fixture file.
///
/// The first lines start with `#` and record the Singular invocation. Then
/// come the variables on one line, the generators, an empty line, and the
/// basis.
fn reference(name: &str) -> Reference {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/rational")
        .join(format!("{name}.txt"));
    let text = std::fs::read_to_string(&path).expect("the fixture file is readable");
    let mut lines = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::trim);
    let variables = lines
        .next()
        .expect("the fixture names its variables")
        .split(',')
        .map(str::to_string)
        .collect();
    let generators: Vec<String> = lines
        .by_ref()
        .take_while(|line| !line.is_empty())
        .map(str::to_string)
        .collect();
    let basis: Vec<String> = lines
        .filter(|line| !line.is_empty())
        .map(str::to_string)
        .collect();
    Reference {
        variables,
        generators,
        basis,
    }
}

fn ring_of(reference: &Reference) -> PolynomialRing<Rationals> {
    PolynomialRing::rationals(&reference.variables).expect("the names are variable names")
}

fn basis_of(
    reference: &Reference,
    options: RationalOptions,
) -> Result<GroebnerBasis<Rationals>, ComputeError> {
    let ring = ring_of(reference);
    let generators = reference
        .generators
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the syntax holds"));
    ring.ideal(generators)
        .expect("the generators share the ring")
        .groebner_basis(options)
}

fn printed(basis: &GroebnerBasis<Rationals>) -> Vec<String> {
    basis.iter().map(ToString::to_string).collect()
}

/// The reference basis, read back through the ring so the comparison is
/// over values and not over the text Singular happened to write.
fn expected(reference: &Reference) -> Vec<String> {
    let ring = ring_of(reference);
    reference
        .basis
        .iter()
        .map(|text| {
            ring.parse_polynomial(text)
                .expect("the reference syntax holds")
                .to_string()
        })
        .collect()
}

#[test]
fn the_reference_families_lift_to_the_singular_basis() {
    for name in [
        "cyclic-3",
        "cyclic-4",
        "cyclic-5",
        "katsura-3",
        "katsura-4",
        "katsura-5",
        "noon-3",
    ] {
        let reference = reference(name);
        let basis = basis_of(&reference, RationalOptions::new()).expect("the system fits");
        assert_eq!(printed(&basis), expected(&reference), "{name}");
        let lift = basis.lift().expect("the driver produced the basis");
        assert_eq!(lift.established, Established::Unchanged);
        assert_eq!(
            lift.primes_folded + lift.primes_discarded,
            lift.primes_consumed,
            "{name}"
        );
        assert!(lift.confirming_primes >= 2, "{name}");
        assert!(lift.modulus_bits > 0, "{name}");
    }
}

#[test]
fn a_hand_system_lifts_to_its_reduced_basis() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
    let ideal = ring
        .ideal([
            ring.parse_polynomial("1/2*x^2 - y")
                .expect("the syntax holds"),
            ring.parse_polynomial("1/3*y^2 - x")
                .expect("the syntax holds"),
        ])
        .expect("the generators share the ring");
    let basis = ideal
        .groebner_basis(RationalOptions::new())
        .expect("the system fits");
    // Singular 4.4.1 on this ideal, largest leading monomial first:
    // x^2-2*y, y^2-3*x.
    assert_eq!(printed(&basis), ["x^2 - 2*y", "y^2 - 3*x"]);
}

#[test]
fn the_zero_ideal_has_the_empty_basis() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
    let ideal = ring
        .ideal([ring.zero(), ring.zero()])
        .expect("the generators share the ring");
    let basis = ideal
        .groebner_basis(RationalOptions::new())
        .expect("the empty ideal needs no prime");
    assert!(basis.is_empty());
    let lift = basis.lift().expect("the driver produced the basis");
    assert_eq!(lift.primes_consumed, 0);
    assert_eq!(lift.modulus_bits, 0);

    // No prime ran, so the report says no run started, whatever the
    // caller asked for.
    let (_, report) = ideal
        .groebner_basis_with_report(
            RationalOptions::new().compute(ComputeOptions::new().threads(8)),
        )
        .expect("the empty ideal needs no prime");
    assert_eq!(report.modular_concurrency, Some(0));
}

#[test]
fn the_unit_ideal_lifts_to_one() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
    let ideal = ring
        .ideal([
            ring.parse_polynomial("x").expect("the syntax holds"),
            ring.parse_polynomial("x + 3").expect("the syntax holds"),
        ])
        .expect("the generators share the ring");
    let basis = ideal
        .groebner_basis(RationalOptions::new())
        .expect("the system fits");
    assert_eq!(printed(&basis), ["1"]);
}

#[test]
fn an_unlucky_first_prime_loses_the_vote_and_the_lift_still_holds() {
    // Modulo the first prime of the sequence the two generators are the
    // same polynomial, so that run leads on x alone. Every later prime
    // sees the ideal (x, y). The vote moves to the prevalent category and
    // the accumulator is built again from the runs that hold it.
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
    let ideal = ring
        .ideal([
            ring.parse_polynomial("x + y").expect("the syntax holds"),
            ring.polynomial([(1, [1, 0]), (FIRST_PRIME + 1, [0, 1])])
                .expect("the exponent vectors match the ring"),
        ])
        .expect("the generators share the ring");
    let basis = ideal
        .groebner_basis(RationalOptions::new())
        .expect("the system fits");
    assert_eq!(printed(&basis), ["x", "y"]);
    let lift = basis.lift().expect("the driver produced the basis");
    assert_eq!(lift.primes_discarded, 1);
    assert_eq!(lift.primes_folded, lift.primes_consumed - 1);
    assert_eq!(lift.primes_skipped, 0);
}

#[test]
fn a_prime_that_divides_a_leading_coefficient_is_skipped_before_it_runs() {
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
    let ideal = ring
        .ideal([
            ring.polynomial([(FIRST_PRIME, [2, 0]), (1, [0, 1])])
                .expect("the exponent vectors match the ring"),
            ring.parse_polynomial("y^2 - 1").expect("the syntax holds"),
        ])
        .expect("the generators share the ring");
    let basis = ideal
        .groebner_basis(RationalOptions::new())
        .expect("the system fits");
    let lift = basis.lift().expect("the driver produced the basis");
    assert_eq!(lift.primes_skipped, 1);
    assert_eq!(lift.primes_discarded, 0);
    // The leading coefficient divides out: the basis is monic.
    assert_eq!(printed(&basis)[0], format!("x^2 + 1/{FIRST_PRIME}*y"));
}

#[test]
fn the_basis_and_the_counters_are_the_same_at_every_thread_count() {
    let reference = reference("katsura-4");
    let mut seen: Option<(Vec<String>, sylvester::ModularLift)> = None;
    for threads in [1usize, 2, 8, 16] {
        let options = RationalOptions::new().compute(ComputeOptions::new().threads(threads));
        let basis = basis_of(&reference, options).expect("the system fits");
        let lift = *basis.lift().expect("the driver produced the basis");
        match &seen {
            None => seen = Some((printed(&basis), lift)),
            Some((held_basis, held_lift)) => {
                assert_eq!(&printed(&basis), held_basis, "{threads} threads");
                assert_eq!(&lift, held_lift, "{threads} threads");
            }
        }
    }
}

/// A thread count no machine can run is clamped, and the report says
/// what the driver used.
#[test]
fn a_thread_count_past_the_cap_is_clamped() {
    let reference = reference("cyclic-3");
    let ring = ring_of(&reference);
    let generators = reference
        .generators
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the syntax holds"));
    let ideal = ring
        .ideal(generators)
        .expect("the generators share the ring");
    let options = RationalOptions::new().compute(ComputeOptions::new().threads(usize::MAX));

    let (basis, report) = ideal
        .groebner_basis_with_report(options)
        .expect("the run holds no memory limit");
    assert_eq!(printed(&basis), expected(&reference));
    assert_eq!(report.modular_concurrency, Some(256));
}

#[test]
fn two_runs_of_one_system_return_one_value() {
    let reference = reference("noon-3");
    let first = basis_of(&reference, RationalOptions::new()).expect("the system fits");
    let second = basis_of(&reference, RationalOptions::new()).expect("the system fits");
    assert_eq!(printed(&first), printed(&second));
    assert_eq!(first.lift(), second.lift());
}

#[test]
fn the_contains_input_rule_establishes_what_it_names() {
    let extra = NonZeroUsize::new(2).expect("2 is not zero");
    for name in ["cyclic-3", "katsura-3", "noon-3"] {
        let reference = reference(name);
        let options = RationalOptions::new().stop(RationalStop::ContainsInput { extra });
        let basis = basis_of(&reference, options).expect("the system fits");
        assert_eq!(printed(&basis), expected(&reference), "{name}");
        let lift = basis.lift().expect("the driver produced the basis");
        assert_eq!(lift.established, Established::ContainsInput, "{name}");

        // T1 is what the run checked: every generator belongs to the
        // ideal the basis generates.
        let ring = ring_of(&reference);
        for text in &reference.generators {
            let generator = ring.parse_polynomial(text).expect("the syntax holds");
            assert!(
                basis
                    .contains(&generator, sylvester::Budget::new())
                    .expect("the division fits"),
                "{name}: {text}"
            );
        }
    }
}

#[test]
fn more_confirming_primes_cost_more_primes_and_return_one_basis() {
    let reference = reference("cyclic-3");
    let two = basis_of(&reference, RationalOptions::new()).expect("the system fits");
    let extra = NonZeroUsize::new(5).expect("5 is not zero");
    let five = basis_of(
        &reference,
        RationalOptions::new().stop(RationalStop::Unchanged { extra }),
    )
    .expect("the system fits");
    assert_eq!(printed(&two), printed(&five));
    let two = two.lift().expect("the driver produced the basis");
    let five = five.lift().expect("the driver produced the basis");
    assert_eq!(two.confirming_primes, 2);
    assert_eq!(five.confirming_primes, 5);
    assert!(five.primes_consumed > two.primes_consumed);
}

#[test]
fn a_passed_deadline_is_a_timeout() {
    let reference = reference("katsura-4");
    let options =
        RationalOptions::new().compute(ComputeOptions::new().timeout(Duration::from_nanos(1)));
    assert_eq!(basis_of(&reference, options), Err(ComputeError::Timeout));
}

#[test]
fn an_exhausted_memory_limit_is_reported_as_itself() {
    let reference = reference("katsura-4");
    let options = RationalOptions::new().compute(ComputeOptions::new().memory_limit(64));
    assert_eq!(
        basis_of(&reference, options),
        Err(ComputeError::MemoryLimitExceeded)
    );
}

/// The smallest memory limit one prime run of `reference` fits in.
///
/// The search is over the prime-field path at the first prime of the
/// sequence, which is the prime the rational run takes first. The engine
/// counts tracked bytes, so the value is a function of the code and the
/// system and not of the machine.
fn one_run_needs(reference: &Reference) -> usize {
    let ring = PolynomialRing::prime_field(FIRST_PRIME as u64, &reference.variables)
        .expect("the first prime of the sequence is a ring modulus");
    let ideal = ring
        .ideal(
            reference
                .generators
                .iter()
                .map(|text| ring.parse_polynomial(text).expect("the syntax holds")),
        )
        .expect("the generators share the ring");
    let fits = |bytes: usize| {
        ideal
            .groebner_basis(ComputeOptions::new().memory_limit(bytes).threads(1))
            .is_ok()
    };
    let mut low = 0usize;
    let mut high = 1usize << 26;
    assert!(fits(high), "the search starts above what one run needs");
    while low + 1 < high {
        let middle = (low + high) / 2;
        if fits(middle) {
            high = middle;
        } else {
            low = middle;
        }
    }
    high
}

/// The ledger charges what the driver holds, not the modular bases alone.
///
/// It counts the rational generators, the category of every retained run,
/// the category table, the primes the sequence found, and the second
/// candidate a lift builds next to the held one. Six times what one prime
/// run needs covers the modular bases and the accumulator, and the rest
/// does not fit under it.
#[test]
fn the_ledger_charges_more_than_the_modular_bases() {
    let reference = reference("katsura-5");
    let needs = one_run_needs(&reference);
    let options =
        RationalOptions::new().compute(ComputeOptions::new().memory_limit(needs * 6).threads(1));
    assert_eq!(
        basis_of(&reference, options),
        Err(ComputeError::MemoryLimitExceeded)
    );
}

#[test]
fn a_run_that_fills_its_memory_share_succeeds_on_the_retry() {
    let reference = reference("katsura-5");
    let needs = one_run_needs(&reference);
    let threads = 8;
    let limit = needs * 7;
    // The first wave splits the limit into eight shares, and one share
    // cannot hold one run. The driver cancels the wave and runs the prime
    // alone with the whole residual.
    assert!(limit / threads < needs, "the share holds a whole run");
    let options =
        RationalOptions::new().compute(ComputeOptions::new().memory_limit(limit).threads(threads));
    let basis = basis_of(&reference, options).expect("the retry holds the whole residual");
    assert_eq!(printed(&basis), expected(&reference));
}

#[test]
fn the_lift_agrees_with_a_prime_the_run_did_not_use() {
    // The prime is below the sequence, which starts at 2^31 - 1 and
    // descends, so no run of the lift saw it. This repeats the production
    // heuristic, so it is a regression test and not an oracle.
    const FRESH: u64 = 32003;
    for name in ["cyclic-5", "katsura-5", "noon-3"] {
        let reference = reference(name);
        let lifted = basis_of(&reference, RationalOptions::new()).expect("the system fits");
        let ring =
            PolynomialRing::prime_field(FRESH, &reference.variables).expect("the modulus is prime");
        let image = |f: &sylvester::Polynomial<Rationals>| {
            ring.polynomial(
                f.terms()
                    .map(|(coeff, exps)| (coeff.clone(), exps.to_vec())),
            )
            .expect("the exponent vectors match the ring")
        };
        let basis = ring
            .ideal(
                reference
                    .generators
                    .iter()
                    .map(|text| ring.parse_polynomial(text).expect("the syntax holds")),
            )
            .expect("the generators share the ring")
            .groebner_basis(ComputeOptions::new())
            .expect("the system fits");
        let reduced: Vec<String> = lifted.iter().map(|f| image(f).to_string()).collect();
        let fresh: Vec<String> = basis.iter().map(ToString::to_string).collect();
        assert_eq!(reduced, fresh, "{name}");
    }
}

type Exps = Vec<u16>;

fn grevlex(a: &[u16], b: &[u16]) -> std::cmp::Ordering {
    let da: u64 = a.iter().map(|&e| u64::from(e)).sum();
    let db: u64 = b.iter().map(|&e| u64::from(e)).sum();
    match da.cmp(&db) {
        std::cmp::Ordering::Equal => {}
        unequal => return unequal,
    }
    for (x, y) in a.iter().zip(b).rev() {
        match x.cmp(y) {
            std::cmp::Ordering::Equal => continue,
            order => return order.reverse(),
        }
    }
    std::cmp::Ordering::Equal
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Oracle {
    terms: BTreeMap<Exps, BigRational>,
}

impl Oracle {
    fn zero() -> Self {
        Oracle {
            terms: BTreeMap::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn leading(&self) -> Option<(Exps, BigRational)> {
        self.terms
            .iter()
            .max_by(|(a, _), (b, _)| grevlex(a, b))
            .map(|(mono, coeff)| (mono.clone(), coeff.clone()))
    }

    /// `self -= factor * x^shift * g`.
    fn sub_scaled(&mut self, factor: &BigRational, shift: &[u16], g: &Oracle) {
        for (mono, coeff) in &g.terms {
            let target: Exps = mono.iter().zip(shift).map(|(a, b)| a + b).collect();
            let entry = self.terms.entry(target).or_insert_with(BigRational::zero);
            *entry -= factor * coeff;
        }
        self.terms.retain(|_, coeff| !coeff.is_zero());
    }

    fn monic(&self) -> Oracle {
        let Some((_, lead)) = self.leading() else {
            return Oracle::zero();
        };
        Oracle {
            terms: self
                .terms
                .iter()
                .map(|(mono, coeff)| (mono.clone(), coeff / &lead))
                .collect(),
        }
    }
}

fn divides(divisor: &[u16], multiple: &[u16]) -> bool {
    divisor.iter().zip(multiple).all(|(a, b)| a <= b)
}

fn oracle_normal_form(f: &Oracle, basis: &[Oracle]) -> Oracle {
    let mut work = f.clone();
    let mut remainder = Oracle::zero();
    while let Some((lead_mono, lead_coeff)) = work.leading() {
        let reducer = basis.iter().find_map(|g| {
            let (mono, coeff) = g.leading()?;
            divides(&mono, &lead_mono).then_some((g, mono, coeff))
        });
        match reducer {
            Some((g, mono, coeff)) => {
                let shift: Exps = lead_mono.iter().zip(&mono).map(|(a, b)| a - b).collect();
                work.sub_scaled(&(&lead_coeff / &coeff), &shift, g);
            }
            None => {
                remainder.terms.insert(lead_mono.clone(), lead_coeff);
                work.terms.remove(&lead_mono);
            }
        }
    }
    remainder
}

fn oracle_s_polynomial(f: &Oracle, g: &Oracle) -> Oracle {
    let (f_mono, f_coeff) = f.leading().expect("the operand is not zero");
    let (g_mono, g_coeff) = g.leading().expect("the operand is not zero");
    let lcm: Exps = f_mono.iter().zip(&g_mono).map(|(a, b)| *a.max(b)).collect();
    let f_shift: Exps = lcm.iter().zip(&f_mono).map(|(a, b)| a - b).collect();
    let g_shift: Exps = lcm.iter().zip(&g_mono).map(|(a, b)| a - b).collect();
    let mut s = Oracle::zero();
    s.sub_scaled(&-f_coeff.recip(), &f_shift, f);
    s.sub_scaled(&g_coeff.recip(), &g_shift, g);
    s
}

/// The reduced Gröbner basis of `generators` over `Q`, largest leading
/// monomial first.
///
/// Textbook Buchberger with no criterion: every pair is processed once.
/// Pairs are taken in increasing degree of their least common multiple,
/// which is the normal selection strategy. The order changes the size of
/// the intermediate basis and not the result.
fn oracle_basis(generators: &[Oracle]) -> Vec<Oracle> {
    let mut basis: Vec<Oracle> = generators
        .iter()
        .filter(|f| !f.is_zero())
        .cloned()
        .collect();
    buchberger(&mut basis);
    let mut current = interreduce(minimal_basis(&basis));
    current.sort_by(|a, b| {
        grevlex(
            &b.leading().expect("the element is not zero").0,
            &a.leading().expect("the element is not zero").0,
        )
    });
    current
}

fn lcm_degree(basis: &[Oracle], left: usize, right: usize) -> u64 {
    let (a, _) = basis[left].leading().expect("the element is not zero");
    let (b, _) = basis[right].leading().expect("the element is not zero");
    a.iter().zip(&b).map(|(x, y)| u64::from(*x.max(y))).sum()
}

fn buchberger(basis: &mut Vec<Oracle>) {
    let mut pairs: BTreeMap<u64, Vec<(usize, usize)>> = BTreeMap::new();
    for left in 0..basis.len() {
        for right in (left + 1)..basis.len() {
            pairs
                .entry(lcm_degree(basis, left, right))
                .or_default()
                .push((left, right));
        }
    }
    while let Some((&degree, _)) = pairs.iter().next() {
        let bucket = pairs.get_mut(&degree).expect("the key was just read");
        let (left, right) = bucket.pop().expect("an empty bucket is removed at once");
        if bucket.is_empty() {
            pairs.remove(&degree);
        }
        let s = oracle_s_polynomial(&basis[left], &basis[right]);
        let remainder = oracle_normal_form(&s, basis);
        if !remainder.is_zero() {
            let index = basis.len();
            basis.push(remainder.monic());
            for other in 0..index {
                pairs
                    .entry(lcm_degree(basis, other, index))
                    .or_default()
                    .push((other, index));
            }
        }
    }
}

fn minimal_basis(basis: &[Oracle]) -> Vec<Oracle> {
    let monic: Vec<Oracle> = basis.iter().map(Oracle::monic).collect();
    let mut minimal: Vec<Oracle> = Vec::new();
    for (index, f) in monic.iter().enumerate() {
        let (lead, _) = f.leading().expect("the element is not zero");
        let redundant = monic.iter().enumerate().any(|(other, g)| {
            if other == index {
                return false;
            }
            let (mono, _) = g.leading().expect("the element is not zero");
            divides(&mono, &lead) && (mono != lead || other < index)
        });
        if !redundant {
            minimal.push(f.clone());
        }
    }
    minimal
}

fn interreduce(mut current: Vec<Oracle>) -> Vec<Oracle> {
    loop {
        let mut next: Vec<Oracle> = Vec::new();
        for index in 0..current.len() {
            let others: Vec<Oracle> = current
                .iter()
                .enumerate()
                .filter(|&(other, _)| other != index)
                .map(|(_, g)| g.clone())
                .collect();
            let remainder = oracle_normal_form(&current[index], &others);
            if !remainder.is_zero() {
                next.push(remainder.monic());
            }
        }
        if next == current {
            break;
        }
        current = next;
    }
    current
}

/// The text of one oracle polynomial, in the syntax `Display` writes.
fn oracle_text(f: &Oracle, variables: &[&str]) -> String {
    let ring = PolynomialRing::rationals(variables).expect("the names are variable names");
    let terms = f
        .terms
        .iter()
        .map(|(mono, coeff)| (coeff.clone(), mono.clone()));
    ring.polynomial(terms)
        .expect("the exponent vectors match the ring")
        .to_string()
}

const VARIABLES: [&str; 3] = ["x", "y", "z"];

/// A random system: 2 or 3 generators over 3 variables, each of at most 3
/// terms of degree at most 2, with coefficients in `[-4, 4]`.
fn systems() -> impl Strategy<Value = Vec<Vec<(i64, [u16; 3])>>> {
    let monomial = prop::sample::select(vec![
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [2, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [0, 2, 0],
        [0, 1, 1],
        [0, 0, 2],
    ]);
    let term = (-4i64..=4, monomial).prop_map(|(coeff, exps)| (coeff, exps));
    let generator = prop::collection::vec(term, 1..=3);
    prop::collection::vec(generator, 2..=3)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// The lifted basis is the one the exact oracle computes.
    #[test]
    fn the_lift_agrees_with_an_exact_buchberger_over_q(system in systems()) {
        let ring = PolynomialRing::rationals(VARIABLES).expect("the names are variable names");
        let generators: Vec<_> = system
            .iter()
            .map(|terms| {
                ring.polynomial(terms.iter().map(|(coeff, exps)| (*coeff, *exps)))
                    .expect("the exponent vectors match the ring")
            })
            .collect();
        let ideal = ring.ideal(generators).expect("the generators share the ring");
        let basis = ideal
            .groebner_basis(RationalOptions::new())
            .expect("a small system fits every budget");

        let oracle: Vec<Oracle> = system
            .iter()
            .map(|terms| {
                let mut poly = Oracle::zero();
                for (coeff, exps) in terms {
                    let entry = poly
                        .terms
                        .entry(exps.to_vec())
                        .or_insert_with(BigRational::zero);
                    *entry += BigRational::from(BigInt::from(*coeff));
                }
                poly.terms.retain(|_, coeff| !coeff.is_zero());
                poly
            })
            .collect();
        let expected: Vec<String> = oracle_basis(&oracle)
            .iter()
            .map(|f| oracle_text(f, &VARIABLES))
            .collect();
        prop_assert_eq!(printed(&basis), expected);
    }
}

#[test]
fn a_basis_with_large_coefficients_lifts_over_many_primes() {
    // The reduced basis of this ideal holds a coefficient no single prime
    // of the sequence can carry, so the lift needs several of them.
    let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
    let big: BigInt = BigInt::from(FIRST_PRIME).pow(3u32) + 1;
    let ideal = ring
        .ideal([
            ring.polynomial([(BigInt::one(), [1, 0]), (-big.clone(), [0, 0])])
                .expect("the exponent vectors match the ring"),
            ring.parse_polynomial("y - 1").expect("the syntax holds"),
        ])
        .expect("the generators share the ring");
    let basis = ideal
        .groebner_basis(RationalOptions::new())
        .expect("the system fits");
    assert_eq!(printed(&basis), [format!("x - {big}"), "y - 1".to_string()]);
    let lift = basis.lift().expect("the driver produced the basis");
    assert!(lift.modulus_bits > big.bits(), "{lift:?}");
}
