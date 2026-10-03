use std::fs;
use std::time::Instant;

use sylvester::verify::{Limits, verify_with_limits};

const CERT_DIR: &str = "baseline";
const RUNS: usize = 21;
const CELLS: [&str; 5] = ["cyclic-5", "cyclic-6", "katsura-5", "eco-8", "noon-5"];

#[test]
#[ignore]
fn record_v2_samples() {
    let limits = Limits {
        max_work_units: u64::MAX,
        ..Limits::default()
    };
    for cell in CELLS {
        let bytes = fs::read(format!("{CERT_DIR}/{cell}.cert")).expect("the certificate exists");
        for repeat in 0..RUNS {
            let start = Instant::now();
            verify_with_limits(&bytes, &limits).expect("the certificate holds");
            println!("sample,{cell},{repeat},{}", start.elapsed().as_secs_f64());
        }
    }
}
