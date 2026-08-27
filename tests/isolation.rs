//! Verifier and engine import isolation.

fn source_files(path: &std::path::Path) -> Vec<(String, String)> {
    let mut files = Vec::new();
    let mut stack = vec![path.to_path_buf()];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir).expect("the directory is readable") {
            let entry = entry.expect("the directory entry is readable").path();
            if entry.is_dir() {
                stack.push(entry);
            } else if entry.extension().is_some_and(|ext| ext == "rs") {
                let text = std::fs::read_to_string(&entry).expect("the file is readable");
                files.push((entry.display().to_string(), text));
            }
        }
    }
    files
}

#[test]
fn the_verifier_shares_no_code_with_the_engines() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/verify");
    let forbidden = [
        "crate::ring",
        "crate::poly",
        "crate::compute",
        "crate::cert",
        "crate::ideal",
        "use super::super",
        "serde",
        "num_bigint",
        "num-bigint",
        "num_integer",
        "num-integer",
        "num_rational",
        "num-rational",
        "num_traits",
        "num-traits",
    ];
    let files = source_files(&root);
    for (path, source) in &files {
        for needle in forbidden {
            assert!(!source.contains(needle), "{path} holds {needle}");
        }
    }
    assert!(
        files.len() >= 5,
        "expected the verifier module tree, found {} files",
        files.len()
    );
}

#[test]
fn the_engines_do_not_reach_into_the_verifier() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let verifier = root.join("verify");
    let allowed = ["compute/mod.rs", "certificate.rs"];
    let mut checked = 0;
    for (path, source) in source_files(&root) {
        if std::path::Path::new(&path).starts_with(&verifier)
            || allowed
                .iter()
                .any(|tail| std::path::Path::new(&path).ends_with(tail))
        {
            continue;
        }
        for (number, line) in source.lines().enumerate() {
            let code = line.trim_start();
            if !code.starts_with("//") {
                assert!(
                    !code.contains("crate::verify"),
                    "{path}:{} reaches into the verifier",
                    number + 1
                );
            }
        }
        checked += 1;
    }
    assert!(
        checked >= 8,
        "expected the engine tree, found {checked} files"
    );
}
