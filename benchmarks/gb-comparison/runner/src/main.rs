//! One-instance-per-invocation runner for the sylvester Gröbner engine.
//!
//! Usage: sylv-runner <input.syl> <f4|classic|certified|f4-certified>
//! <threads> <modulus>
//!
//! `threads` goes to `ComputeOptions::threads`. A count of 1 keeps the
//! whole computation on the calling thread. Above 1, F4 reduces the rows
//! of a large batch in a pool of that size.
//!
//! `certified` runs the classic backend through
//! [`sylvester::Ideal::groebner_basis_certified`] and writes a
//! `sylv-gb-cert-v1` certificate. `f4-certified` runs the F4 backend
//! through the same call and writes a `sylv-gb-cert-v2` certificate.
//! Both report the certificate bytes and a standalone re-verification
//! time. `TIME_S` covers the engine, the certificate emitter, and the
//! bundled verifier. It is not comparable to the raw `TIME_S` of the
//! same backend. `driver.py` computes `split` from the two cells.
//!
//! Input format: line 1 is nvars. Each further line is one polynomial as
//! terms separated by `;`, each term `coeff,e1,e2,...,en`.
//!
//! Output: `TIME_S` (wall seconds of the compute call), `SIZE`, one
//! `LM e1 e2 ... en` line per basis element (grevlex leading monomial,
//! x1 > x2 > ... > xn, independent of engine term order), then one
//! `POLY <term count>` block per element with terms sorted grevlex
//! descending. Coefficients are the crate's reduced values, not monic.
//! `driver.py` applies the same monic normalization to every engine.
//! Raw modes add one `COUNTERS <key>=<value> ...` line from
//! `Ideal::groebner_basis_with_report`. Certified modes add `CERT_BYTES`
//! and `VERIFY_S`. Every mode prints `PEAK_RSS_KB` last, from
//! `/proc/self/status`.

use std::cmp::Ordering;
use std::time::{Duration, Instant};

use sylvester::{Backend, ComputeError, ComputeOptions, ComputeReport, Polynomial, PolynomialRing};

fn grevlex_cmp(a: &[u16], b: &[u16]) -> Ordering {
    let da: u32 = a.iter().map(|&e| e as u32).sum();
    let db: u32 = b.iter().map(|&e| e as u32).sum();
    match da.cmp(&db) {
        Ordering::Equal => {
            for (x, y) in a.iter().zip(b.iter()).rev() {
                match x.cmp(y) {
                    Ordering::Equal => continue,
                    ord => return ord.reverse(),
                }
            }
            Ordering::Equal
        }
        ord => ord,
    }
}

fn print_basis(polys: &[Polynomial]) {
    for poly in polys {
        let lm = poly
            .terms()
            .map(|(_, exps)| exps)
            .max_by(|a, b| grevlex_cmp(a, b))
            .expect("nonzero poly");
        let s: Vec<String> = lm.iter().map(|e| e.to_string()).collect();
        println!("LM {}", s.join(" "));
    }
    for poly in polys {
        let mut terms: Vec<(u64, &[u16])> = poly.terms().collect();
        terms.sort_by(|(_, a), (_, b)| grevlex_cmp(b, a));
        println!("POLY {}", terms.len());
        for (coeff, exps) in terms {
            let s: Vec<String> = exps.iter().map(|e| e.to_string()).collect();
            println!("{coeff} {}", s.join(" "));
        }
    }
}

/// Peak resident set size in kibibytes, from the `VmHWM` line of
/// `/proc/self/status`. `None` off Linux or if the line is missing.
fn peak_rss_kb() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            let digits: String = rest.chars().filter(|c| c.is_ascii_digit()).collect();
            return digits.parse().ok();
        }
    }
    None
}

/// One `COUNTERS <key>=<value> ...` line from the run's report.
///
/// A backend that counts nothing omits those keys. The keys are the field
/// names of `sylvester::F4Counters`. `threads_used` is a field of
/// `sylvester::ComputeReport`.
fn print_counters(report: &ComputeReport) {
    let backend = match report.backend {
        Backend::F4 => "f4",
        Backend::Classic => "classic",
    };
    let mut line = format!("COUNTERS backend={backend}");
    line.push_str(&format!(" threads_used={}", report.threads_used));
    if let Some(c) = report.counters {
        for (key, value) in [
            ("batches", c.batches),
            ("matrix_rows", c.matrix_rows),
            ("matrix_columns", c.matrix_columns),
            ("matrix_nonzeros", c.matrix_nonzeros),
            ("zero_rows", c.zero_rows),
            ("new_pivots", c.new_pivots),
            ("lane_restarts", u64::from(c.lane_restarts)),
            ("batch_retries", c.batch_retries),
            ("basis_monomials", c.basis_monomials as u64),
            ("pairs_generated", c.pairs_generated),
            ("pairs_discarded_product", c.pairs_discarded_product),
            ("pairs_discarded_b", c.pairs_discarded_b),
            ("pairs_discarded_m", c.pairs_discarded_m),
            ("pairs_discarded_f", c.pairs_discarded_f),
        ] {
            line.push_str(&format!(" {key}={value}"));
        }
    }
    println!("{line}");
}

fn print_peak_rss() {
    match peak_rss_kb() {
        Some(kb) => println!("PEAK_RSS_KB {kb}"),
        None => println!("PEAK_RSS_KB null"),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!(
            "usage: sylv-runner <input.syl> <f4|classic|certified|f4-certified> <threads> <modulus>"
        );
        std::process::exit(2);
    }
    let text = std::fs::read_to_string(&args[1]).expect("read input");
    let mode = args[2].as_str();
    let threads: usize = args[3].parse().expect("threads");
    let modulus: u64 = args[4].parse().expect("modulus");

    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let nvars: usize = lines
        .next()
        .expect("nvars line")
        .trim()
        .parse()
        .expect("nvars");
    let names: Vec<String> = (0..nvars).map(|index| format!("x{index}")).collect();
    let ring = PolynomialRing::prime_field(modulus, names).expect("prime modulus");

    let mut polys = Vec::new();
    for line in lines {
        let mut terms = Vec::new();
        for term in line.trim().split(';') {
            let mut it = term.split(',');
            let coeff: i64 = it.next().expect("coeff").parse().expect("coeff");
            let exps: Vec<u16> = it.map(|e| e.parse().expect("exp")).collect();
            assert_eq!(exps.len(), nvars, "bad exponent vector length");
            terms.push((coeff, exps));
        }
        polys.push(ring.polynomial(terms).expect("exponent vector width"));
    }
    let ideal = ring.ideal(polys).expect("one ring");

    let (backend, certified) = match mode {
        "f4" => (Backend::F4, false),
        "classic" => (Backend::Classic, false),
        "certified" => (Backend::Classic, true),
        "f4-certified" => (Backend::F4, true),
        other => {
            eprintln!("unknown mode {other}");
            std::process::exit(2);
        }
    };

    if certified {
        let options = ComputeOptions::new()
            .backend(backend)
            .threads(threads)
            .timeout(Duration::from_secs(120));
        let start = Instant::now();
        let result = ideal.groebner_basis_certified(options);
        let elapsed = start.elapsed().as_secs_f64();
        match result {
            Ok(run) => {
                let (basis, bytes) = run.into_parts();
                println!("TIME_S {elapsed}");
                println!("SIZE {}", basis.len());
                print_basis(&basis);
                println!("CERT_BYTES {}", bytes.len());
                // The default caps guard against hostile bytes. This
                // certificate came from the call above, so only the work
                // cap is opened; every structural cap stays.
                let limits = sylvester::verify::Limits {
                    max_work_units: u64::MAX,
                    ..sylvester::verify::Limits::default()
                };
                let vstart = Instant::now();
                sylvester::verify::verify_with_limits(&bytes, &limits)
                    .expect("emitted certificate must verify");
                let verify_s = vstart.elapsed().as_secs_f64();
                println!("VERIFY_S {verify_s}");
                print_peak_rss();
            }
            Err(sylvester::CertifyError::Engine(ComputeError::Timeout)) => {
                println!("STATUS TIMEOUT");
                println!("TIME_S {elapsed}");
                print_peak_rss();
            }
            Err(sylvester::CertifyError::VerifierExhausted(e)) if e.is_exhaustion() => {
                println!("STATUS TIMEOUT");
                println!("TIME_S {elapsed}");
                print_peak_rss();
            }
            Err(e) => {
                println!("STATUS ERROR {e}");
                std::process::exit(1);
            }
        }
        return;
    }

    let options = ComputeOptions::new()
        .backend(backend)
        .threads(threads)
        .timeout(Duration::from_secs(120));
    let start = Instant::now();
    let result = ideal.groebner_basis_with_report(options);
    let elapsed = start.elapsed().as_secs_f64();

    match result {
        Ok((gb, report)) => {
            println!("TIME_S {elapsed}");
            println!("SIZE {}", gb.len());
            print_basis(&gb);
            print_counters(&report);
            print_peak_rss();
        }
        Err(ComputeError::Timeout) => {
            println!("STATUS TIMEOUT");
            println!("TIME_S {elapsed}");
            print_peak_rss();
        }
        Err(e) => {
            println!("STATUS ERROR {e}");
            std::process::exit(1);
        }
    }
}
