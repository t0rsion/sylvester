//! What the memory limit and the deadline stop.

use sylvester::{Backend, ComputeError, ComputeOptions, Ideal, PolynomialRing};

/// The same three-variable system, padded into a ring of `nvars`
/// variables.
///
/// The basis and every term count are the same for every width, so a
/// memory count that differs between two widths differs only in what one
/// exponent vector costs.
fn padded_ideal(nvars: usize) -> Ideal {
    let names: Vec<String> = (0..nvars).map(|index| format!("x{index}")).collect();
    let ring = PolynomialRing::prime_field(32003, names).expect("32003 is prime");
    let generators: Vec<_> = ["x0^2*x1 - 3*x2 + 1", "x0*x2 - x1", "x1^3 + x2^2 - x0"]
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the text parses"))
        .collect();
    ring.ideal(generators).expect("one ring")
}

/// The smallest memory limit, in bytes, at which the run finishes.
///
/// A run that finishes at one limit finishes at every larger one, because
/// every check compares a byte count against the limit, so the search is a
/// bisection.
fn smallest_limit_that_finishes(ideal: &Ideal, backend: Backend) -> usize {
    let finishes = |limit: usize| {
        ideal
            .groebner_basis(ComputeOptions::new().backend(backend).memory_limit(limit))
            .is_ok()
    };
    let mut high = 1usize;
    while !finishes(high) {
        high *= 2;
        assert!(high < (1 << 30), "the run needs more than a gigabyte");
    }
    let mut low = high / 2;
    while low + 1 < high {
        let mid = low + (high - low) / 2;
        if finishes(mid) {
            high = mid;
        } else {
            low = mid;
        }
    }
    high
}

/// A ring of more than eleven variables spills every exponent vector to
/// the heap, and the memory count must charge those bytes. Eleven
/// variables still fit inline, so the same system costs strictly less
/// there.
#[test]
fn a_wide_ring_charges_the_heap_exponents() {
    for backend in [Backend::Classic, Backend::F4] {
        let inline = smallest_limit_that_finishes(&padded_ideal(11), backend);
        let spilled = smallest_limit_that_finishes(&padded_ideal(12), backend);
        assert!(
            spilled > inline,
            "backend {backend:?}: 12 variables cost {spilled} bytes, 11 cost {inline}"
        );
    }
}

/// A limit of zero bytes stops every run over a generator, on either
/// backend. An empty generator list holds nothing, so it finishes under
/// any limit.
#[test]
fn a_zero_memory_limit_stops_every_run() {
    let ring = PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime");
    let inputs: Vec<Vec<_>> = vec![
        vec![ring.zero()],
        vec![ring.parse_polynomial("x^2 - y").expect("parses")],
        vec![
            ring.parse_polynomial("x^2 - y").expect("parses"),
            ring.parse_polynomial("x*y - 1").expect("parses"),
        ],
    ];
    for generators in inputs {
        let ideal = ring.ideal(generators.clone()).expect("one ring");
        for backend in [Backend::Classic, Backend::F4] {
            assert_eq!(
                ideal.groebner_basis(ComputeOptions::new().backend(backend).memory_limit(0)),
                Err(ComputeError::MemoryLimitExceeded),
                "backend {backend:?}, {} generators",
                generators.len()
            );
        }
    }

    let empty = ring.ideal(Vec::new()).expect("one ring");
    for backend in [Backend::Classic, Backend::F4] {
        assert!(
            empty
                .groebner_basis(ComputeOptions::new().backend(backend).memory_limit(0))
                .is_ok(),
            "backend {backend:?}"
        );
    }
}
