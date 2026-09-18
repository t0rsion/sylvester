#[path = "../src/io.rs"]
mod cli_io;

use std::io::{self, Cursor};
use std::path::Path;
use std::time::{Duration, Instant};

use cli_io::{ReadError, ReadLimits, TargetConflict};

#[test]
fn a_read_stops_at_the_source_limit() {
    let error = cli_io::read_bytes(Cursor::new(b"abcdef"), ReadLimits::new(None, None, Some(5)))
        .expect_err("the sixth byte exceeds the limit");
    assert!(matches!(error, ReadError::TooLong { limit: 5 }));
}

#[test]
fn a_read_stops_at_the_memory_limit() {
    let error = cli_io::read_bytes(Cursor::new(b"abcdef"), ReadLimits::new(None, Some(5), None))
        .expect_err("the sixth byte exceeds the limit");
    assert!(matches!(error, ReadError::MemoryLimitExceeded));
}

#[test]
fn a_read_checks_an_expired_deadline() {
    let error = cli_io::read_bytes(
        Cursor::new(b"x"),
        ReadLimits::new(Some(Instant::now() - Duration::from_secs(1)), None, None),
    )
    .expect_err("the deadline has passed");
    assert!(matches!(error, ReadError::Timeout));
}

#[test]
fn a_text_read_reports_invalid_utf8() {
    let error = cli_io::read_text(Cursor::new([0xff]), ReadLimits::default())
        .expect_err("invalid UTF-8 is a source error");
    assert!(matches!(error, ReadError::InvalidUtf8));
}

#[test]
fn a_path_read_checks_metadata_and_decodes_text() {
    let directory = temporary_directory("path-read");
    let source = directory.join("input.txt");
    std::fs::write(&source, b"hello").expect("the source is writable");
    let error = cli_io::read_path_bytes(&source, ReadLimits::new(None, None, Some(4)))
        .expect_err("metadata exceeds the limit");
    assert!(matches!(error, ReadError::TooLong { limit: 4 }));
    assert_eq!(
        cli_io::read_path_text(&source, ReadLimits::default()).expect("the source is UTF-8"),
        "hello"
    );
    remove_directory(&directory);
}

#[test]
fn an_atomic_write_leaves_the_old_file_on_render_error() {
    let directory = temporary_directory("write-error");
    let destination = directory.join("result.txt");
    std::fs::write(&destination, b"old").expect("the destination is writable");
    let error =
        cli_io::write_atomic_with(&destination, |_file| Err(io::Error::other("render failed")))
            .expect_err("the renderer failed");
    assert_eq!(error.kind(), io::ErrorKind::Other);
    assert_eq!(
        std::fs::read(&destination).expect("the destination remains"),
        b"old"
    );
    assert_no_temporary_files(&directory);
    remove_directory(&directory);
}

#[test]
fn an_atomic_write_replaces_the_destination() {
    let directory = temporary_directory("write-success");
    let destination = directory.join("result.txt");
    std::fs::write(&destination, b"old").expect("the destination is writable");
    cli_io::write_atomic(&destination, b"new").expect("the replacement succeeds");
    assert_eq!(
        std::fs::read(&destination).expect("the destination exists"),
        b"new"
    );
    assert_no_temporary_files(&directory);
    remove_directory(&directory);
}

#[test]
fn lexical_aliases_are_rejected_as_one_target() {
    let directory = temporary_directory("lexical");
    let destination = directory.join("result.txt");
    let alias = directory.join("nested").join("..").join("result.txt");
    std::fs::create_dir(directory.join("nested")).expect("the nested directory is writable");
    std::fs::write(&destination, b"old").expect("the destination exists");
    assert!(matches!(
        cli_io::reject_same_targets(Some(&destination), Some(&alias)),
        Err(TargetConflict { .. })
    ));
    remove_directory(&directory);
}

#[cfg(unix)]
#[test]
fn existing_symlink_aliases_are_rejected_as_one_target() {
    use std::os::unix::fs::symlink;

    let directory = temporary_directory("symlink");
    let destination = directory.join("result.txt");
    let alias = directory.join("alias.txt");
    std::fs::write(&destination, b"old").expect("the destination exists");
    symlink(&destination, &alias).expect("the alias is writable");
    assert!(matches!(
        cli_io::reject_same_targets(Some(&destination), Some(&alias)),
        Err(TargetConflict { .. })
    ));
    remove_directory(&directory);
}

#[cfg(unix)]
#[test]
fn missing_files_under_one_existing_symlink_parent_are_rejected() {
    use std::os::unix::fs::symlink;

    let directory = temporary_directory("missing-symlink-parent");
    let real = directory.join("real");
    let alias = directory.join("alias");
    std::fs::create_dir(&real).expect("the real directory is writable");
    symlink(&real, &alias).expect("the alias is writable");
    let first = real.join("result.txt");
    let second = alias.join("result.txt");
    assert!(matches!(
        cli_io::reject_same_targets(Some(&first), Some(&second)),
        Err(TargetConflict { .. })
    ));
    remove_directory(&directory);
}

#[cfg(unix)]
#[test]
fn existing_hardlink_aliases_are_rejected_as_one_target() {
    let directory = temporary_directory("hardlink");
    let destination = directory.join("result.txt");
    let alias = directory.join("alias.txt");
    std::fs::write(&destination, b"old").expect("the destination exists");
    std::fs::hard_link(&destination, &alias).expect("the alias is writable");
    assert!(matches!(
        cli_io::reject_same_targets(Some(&destination), Some(&alias)),
        Err(TargetConflict { .. })
    ));
    remove_directory(&directory);
}

#[cfg(windows)]
#[test]
fn windows_case_aliases_are_rejected_as_one_target() {
    let directory = temporary_directory("case");
    let first = directory.join("Foo.txt");
    let second = directory.join("foo.txt");
    std::fs::write(&first, b"old").expect("the destination exists");
    assert!(matches!(
        cli_io::reject_same_targets(Some(&first), Some(&second)),
        Err(TargetConflict { .. })
    ));
    remove_directory(&directory);
}

fn temporary_directory(label: &str) -> std::path::PathBuf {
    let directory =
        std::env::temp_dir().join(format!("sylvester-cli-io-{label}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&directory);
    std::fs::create_dir(&directory).expect("the temporary directory is writable");
    directory
}

fn assert_no_temporary_files(directory: &Path) {
    let entries = std::fs::read_dir(directory)
        .expect("the temporary directory is readable")
        .collect::<Result<Vec<_>, _>>()
        .expect("the directory entries are readable");
    assert!(
        entries
            .iter()
            .all(|entry| { !entry.file_name().to_string_lossy().contains(".sylv-") })
    );
}

fn remove_directory(directory: &Path) {
    std::fs::remove_dir_all(directory).expect("the temporary directory is removable");
}
