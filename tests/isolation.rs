//! The module import isolation scan.
//!
//! Both verifier trees are code-isolated from the engines. The scan reads
//! every file under `src/verify` and fails on an import the isolation
//! rule forbids. `AGENTS.md` states the rule; this test is what holds it.
//!
//! The scan holds two more rules, both from `docs/rational-design.md`. The
//! verifier gains no dependency (section 10): `crate::ring::rational` and
//! `crate::compute::modular` read rational coefficients through
//! `num-bigint` and `num-rational`, and no file under `src/verify` names
//! either crate. The verifiers keep their own division (section 5): none
//! of them calls `crate::normal_form`. Nothing in them reads a leading
//! monomial ideal through `crate::hilbert` either.

#[test]
fn the_verifier_shares_no_code_with_the_engines() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/verify");
    let forbidden = [
        "crate::ring",
        "crate::poly",
        "crate::compute",
        "crate::cert",
        "crate::ideal",
        "crate::normal_form",
        "crate::hilbert",
        "crate::compute::modular",
        "use super::super",
        "serde",
        "num_bigint",
        "num_rational",
        "num_integer",
        "num_traits",
    ];
    let mut checked = 0;
    let mut v2_files = 0;
    let mut stack = vec![root.clone()];
    while let Some(path) = stack.pop() {
        for entry in std::fs::read_dir(&path).expect("the verifier directory is readable") {
            let entry = entry.expect("the directory entry is readable").path();
            if entry.is_dir() {
                stack.push(entry);
                continue;
            }
            let text = std::fs::read_to_string(&entry).expect("the file is readable");
            for needle in forbidden {
                assert!(!text.contains(needle), "{} holds {needle}", entry.display());
            }
            if entry.parent() == Some(root.join("v2").as_path()) {
                v2_files += 1;
            }
            checked += 1;
        }
    }
    assert!(
        checked >= 5,
        "expected the verifier module tree, found {checked} files"
    );
    assert!(
        v2_files >= 5,
        "expected the v2 verifier files in the scan, found {v2_files}"
    );
}
