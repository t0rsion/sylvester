//! Cancellation exhausts verification without asserting invalidity.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use sylvester::verify::{Limits, VerifyError, verify_with_limits};

#[test]
fn either_certificate_format_observes_a_shared_stop_flag() {
    let v1 = br#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":0,"input":[],"basis":[],"origin":[],"membership":[],"spairs":[]}"#;
    let v2 = include_bytes!("fixtures/cert-v2/tiny.cert");
    for bytes in [v1.as_slice(), v2.as_slice()] {
        let flag = Arc::new(AtomicBool::new(false));
        let limits = Limits {
            cancellation: Some(Arc::clone(&flag)),
            ..Limits::default()
        };
        assert!(verify_with_limits(bytes, &limits).is_ok());
        flag.store(true, Ordering::Relaxed);
        assert_eq!(
            verify_with_limits(bytes, &limits),
            Err(VerifyError::DeadlineExceeded)
        );
    }
}

#[test]
fn cancellation_outranks_an_unknown_format() {
    let limits = Limits {
        cancellation: Some(Arc::new(AtomicBool::new(true))),
        ..Limits::default()
    };
    for bytes in [b"".as_slice(), b"invalid certificate"] {
        assert_eq!(
            verify_with_limits(bytes, &limits),
            Err(VerifyError::DeadlineExceeded)
        );
    }
}
