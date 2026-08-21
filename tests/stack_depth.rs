//! The classic backend runs eco-9 in 1 MiB of stack.
//!
//! eco-9 does not finish inside the time limit. The test therefore accepts
//! `Timeout` as well as a basis: it checks stack use, not completion. The
//! 1 MiB worker stack is one eighth of the 8 MiB glibc default (`ulimit -s`
//! on the recorded machine), so a run that grows its stack with the size of
//! the input fails here first.

use std::time::Duration;
use sylvester::{Backend, ComputeError, ComputeOptions, Ideal, PolynomialRing};

const ECO_9: &str = include_str!("inputs/eco-9.syl");

const MODULUS: u64 = 1_073_741_827;

const WORKER_STACK_BYTES: usize = 1 << 20;

const TIMEOUT: Duration = Duration::from_secs(3);

/// Parse the runner input format: line 1 is the variable count, each
/// further line is one polynomial as `;`-separated `coeff,e1,...,en` terms.
fn parse_syl(text: &str, modulus: u64) -> Ideal {
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let nvars: usize = lines
        .next()
        .expect("input must start with a variable count")
        .trim()
        .parse()
        .expect("variable count must be an integer");

    let names: Vec<String> = (0..nvars).map(|index| format!("x{index}")).collect();
    let ring = PolynomialRing::prime_field(modulus, names).expect("the modulus is prime");

    let generators: Vec<_> = lines
        .map(|line| {
            let terms: Vec<(i64, Vec<u16>)> = line
                .trim()
                .split(';')
                .map(|term| {
                    let mut fields = term.split(',');
                    let coeff: i64 = fields
                        .next()
                        .expect("term must start with a coefficient")
                        .parse()
                        .expect("coefficient must be an integer");
                    let exps: Vec<u16> = fields
                        .map(|exp| exp.parse().expect("exponent must be an integer"))
                        .collect();
                    (coeff, exps)
                })
                .collect();
            ring.polynomial(terms)
                .expect("every exponent vector must have nvars entries")
        })
        .collect();
    ring.ideal(generators)
        .expect("the polynomials share a ring")
}

#[test]
fn classic_backend_runs_eco_9_in_one_mib_of_stack() {
    let ideal = parse_syl(ECO_9, MODULUS);
    assert_eq!(ideal.generators().len(), 9, "eco-9 has nine generators");

    let worker = std::thread::Builder::new()
        .stack_size(WORKER_STACK_BYTES)
        .spawn(move || {
            ideal.groebner_basis(
                ComputeOptions::new()
                    .backend(Backend::Classic)
                    .timeout(TIMEOUT),
            )
        })
        .expect("worker thread must start");

    match worker.join().expect("worker thread must not panic") {
        Ok(_) | Err(ComputeError::Timeout) => {}
        Err(other) => panic!("classic backend on eco-9 reported {other}"),
    }
}
