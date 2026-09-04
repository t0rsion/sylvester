//! Golden tests for the `sylv` binary.
//!
//! Each test runs the built binary through `CARGO_BIN_EXE_sylv` on a
//! committed input under `tests/inputs`, and asserts the standard output
//! and the exit code. The codes are the ones `docs/rational-design.md` section
//! 8.4 lists.

use std::env;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// What one run of `sylv` printed.
struct Run {
    code: i32,
    out: String,
    err: String,
}

/// Run `sylv` with no standard input.
fn sylv(args: &[&str]) -> Run {
    run(args, None)
}

/// Run `sylv` on text it reads from standard input.
fn sylv_stdin(args: &[&str], input: &str) -> Run {
    run(args, Some(input))
}

fn run(args: &[&str], input: Option<&str>) -> Run {
    let mut child = Command::new(env!("CARGO_BIN_EXE_sylv"))
        .args(args)
        .stdin(if input.is_some() {
            Stdio::piped()
        } else {
            Stdio::null()
        })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("the binary starts");
    if let Some(text) = input {
        child
            .stdin
            .take()
            .expect("a pipe on standard input")
            .write_all(text.as_bytes())
            .expect("the binary reads standard input");
    }
    let output = child.wait_with_output().expect("the binary exits");
    Run {
        code: output.status.code().expect("the binary exits with a code"),
        out: String::from_utf8(output.stdout).expect("UTF-8 on standard output"),
        err: String::from_utf8(output.stderr).expect("UTF-8 on standard error"),
    }
}

/// A committed input.
fn input(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/inputs")
        .join(name)
        .to_str()
        .expect("the path is UTF-8")
        .to_string()
}

/// A file this test writes. Every test names its own.
fn scratch(name: &str) -> PathBuf {
    let directory = env::temp_dir().join("sylv-cli-tests");
    fs::create_dir_all(&directory).expect("the scratch directory");
    let path = directory.join(name);
    let _ = fs::remove_file(&path);
    path
}

/// The path as the command line takes it.
fn text_of(path: &Path) -> String {
    path.to_str().expect("the path is UTF-8").to_string()
}

/// The reduced basis of cyclic-3 over `F_32003`, in the ring of
/// `cyclic3.text`.
const CYCLIC3: &str = "\
# vars: x, y, z
# modulus: 32003
z^3 + 32002
y^2 + y*z + z^2
x + y + z
";

/// The reduced basis of cyclic-3 over `Q`.
const CYCLIC3_Q: &str = "\
# vars: x, y, z
# coefficients: rationals
z^3 - 1
y^2 + y*z + z^2
x + y + z
";

/// The reduced basis of cyclic-3 over `F_32003`, in the ring the `ms` and
/// `syl` inputs name.
const CYCLIC3_SYNTHETIC: &str = "\
# vars: x1, x2, x3
# modulus: 32003
x3^3 + 32002
x2^2 + x2*x3 + x3^2
x1 + x2 + x3
";

#[test]
fn gb_reads_a_text_system() {
    let run = sylv(&["gb", &input("cyclic3.text")]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3);
}

#[test]
fn gb_reads_a_text_system_over_the_rationals() {
    let run = sylv(&["gb", &input("cyclic3-q.text")]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3_Q);
}

#[test]
fn gb_reads_an_ms_system() {
    let run = sylv(&["gb", &input("cyclic3.ms")]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3_SYNTHETIC);
}

#[test]
fn gb_reads_an_ms_system_over_the_rationals() {
    let run = sylv(&["gb", &input("cyclic3-q.ms")]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(
        run.out,
        "\
# vars: x1, x2, x3
# coefficients: rationals
x3^3 - 1
x2^2 + x2*x3 + x3^2
x1 + x2 + x3
"
    );
}

#[test]
fn gb_reads_a_syl_system_under_a_named_modulus() {
    let run = sylv(&["gb", &input("cyclic3.syl"), "--modulus", "32003"]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3_SYNTHETIC);
}

#[test]
fn gb_reads_a_benchmark_syl_input() {
    let run = sylv(&["gb", &input("cyclic4.syl"), "--modulus", "32003"]);
    assert_eq!(run.code, 0, "{}", run.err);
    let lines: Vec<&str> = run.out.lines().collect();
    assert_eq!(lines[0], "# vars: x1, x2, x3, x4");
    assert_eq!(lines[1], "# modulus: 32003");
    assert_eq!(lines.len(), 9, "cyclic-4 has 7 basis elements");
}

#[test]
fn gb_reads_standard_input_as_text() {
    let run = sylv_stdin(&["gb"], CYCLIC3);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3, "the basis of a basis is the basis");
}

#[test]
fn gb_writes_the_ms_format() {
    let run = sylv(&["gb", &input("cyclic3.text"), "--out-format", "ms"]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(
        run.out,
        "\
x,y,z
32003
z^3+32002,
y^2+y*z+z^2,
x+y+z
"
    );
}

#[test]
fn gb_writes_the_syl_format() {
    let run = sylv(&["gb", &input("cyclic3.text"), "--out-format", "syl"]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(
        run.out,
        "\
3
1,0,0,3;32002,0,0,0
1,0,2,0;1,0,1,1;1,0,0,2
1,1,0,0;1,0,1,0;1,0,0,1
"
    );
}

#[test]
fn gb_writes_rational_coefficients_in_the_syl_format() {
    let run = sylv(&["gb", &input("halves.text"), "--out-format", "syl"]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(
        run.out,
        "\
2
1,2,0;-2,0,1
1,1,1;-1/3,0,0
1,0,2;-1/6,1,0
"
    );
}

#[test]
fn ms_output_reads_back_as_the_same_basis() {
    let path = scratch("round-trip.ms");
    let run = sylv(&[
        "gb",
        &input("cyclic3.text"),
        "--out-format",
        "ms",
        "-o",
        &text_of(&path),
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    let back = sylv(&["gb", &text_of(&path)]);
    assert_eq!(back.code, 0, "{}", back.err);
    assert_eq!(back.out, CYCLIC3);
}

#[test]
fn gb_reports_the_counters_on_standard_error() {
    let run = sylv(&["gb", &input("cyclic3.text"), "--report"]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3);
    assert!(run.err.contains("backend: f4"), "{}", run.err);
    assert!(run.err.contains("pairs_generated: "), "{}", run.err);
}

#[test]
fn gb_reports_the_lift_over_the_rationals() {
    let run = sylv(&["gb", &input("cyclic3-q.text"), "--report"]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3_Q);
    assert!(
        run.err.contains("established: contains-input"),
        "{}",
        run.err
    );
}

#[test]
fn gb_runs_the_classic_backend() {
    let run = sylv(&["gb", &input("cyclic3.text"), "--backend", "classic"]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3);
}

#[test]
fn certify_writes_a_certificate_the_verifier_accepts() {
    let path = scratch("f4.cert");
    let run = sylv(&[
        "certify",
        &input("cyclic3.text"),
        "--certificate",
        &text_of(&path),
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3);
    let bytes = fs::read(&path).expect("the certificate file");
    assert_eq!(bytes[0], b'S', "the F4 backend writes sylv-gb-cert-v2");

    let verified = sylv(&["verify", &text_of(&path)]);
    assert_eq!(verified.code, 0, "{}", verified.err);
    assert_eq!(verified.out, CYCLIC3_SYNTHETIC);
    assert!(
        verified.err.contains("synthetic names x1 to x3"),
        "{}",
        verified.err
    );
}

#[test]
fn certify_writes_a_v1_certificate_under_the_classic_backend() {
    let path = scratch("classic.cert");
    let run = sylv(&[
        "certify",
        &input("cyclic3.text"),
        "--backend",
        "classic",
        "--certificate",
        &text_of(&path),
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3);
    let bytes = fs::read(&path).expect("the certificate file");
    assert_eq!(bytes[0], b'{', "the classic backend writes sylv-gb-cert-v1");

    let verified = sylv(&["verify", &text_of(&path)]);
    assert_eq!(verified.code, 0, "{}", verified.err);
    assert_eq!(verified.out, CYCLIC3_SYNTHETIC);
}

#[test]
fn gb_certified_writes_no_file() {
    let run = sylv(&["gb", &input("cyclic3.text"), "--certified"]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3);
}

#[test]
fn verify_rejects_a_flipped_byte() {
    let path = scratch("flipped.cert");
    let run = sylv(&[
        "certify",
        &input("cyclic3.text"),
        "--certificate",
        &text_of(&path),
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    let mut bytes = fs::read(&path).expect("the certificate file");
    let last = bytes.len() - 5;
    bytes[last] ^= 1;
    fs::write(&path, &bytes).expect("the flipped certificate");

    let verified = sylv(&["verify", &text_of(&path)]);
    assert_eq!(verified.code, 1, "{}", verified.err);
    assert!(verified.out.is_empty());
    assert!(verified.err.contains("rejected"), "{}", verified.err);
}

#[test]
fn verify_prints_nothing_when_quiet() {
    let path = scratch("quiet.cert");
    let run = sylv(&[
        "certify",
        &input("cyclic3.text"),
        "--certificate",
        &text_of(&path),
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    let verified = sylv(&["verify", &text_of(&path), "--quiet"]);
    assert_eq!(verified.code, 0, "{}", verified.err);
    assert!(verified.out.is_empty());
    assert!(verified.err.is_empty());
}

#[test]
fn verify_rejects_a_file_past_the_byte_cap() {
    let path = scratch("capped.cert");
    let run = sylv(&[
        "certify",
        &input("cyclic3.text"),
        "--certificate",
        &text_of(&path),
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    let verified = sylv(&["verify", &text_of(&path), "--max-bytes", "10"]);
    assert_eq!(verified.code, 3, "{}", verified.err);
}

#[test]
fn certify_over_the_rationals_is_a_usage_error() {
    let run = sylv(&["certify", &input("cyclic3-q.text")]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("no certified path"), "{}", run.err);
}

#[test]
fn normal_form_prints_the_remainder() {
    let run = sylv(&[
        "normal-form",
        "--basis",
        &input("cyclic3-basis.text"),
        "--poly",
        "x^2",
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, "# vars: x, y, z\n# modulus: 32003\ny*z\n");
}

#[test]
fn normal_form_prints_the_remainder_over_the_rationals() {
    let run = sylv(&[
        "nf",
        "--basis",
        &input("cyclic3-basis-q.text"),
        "--poly",
        "x^2",
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, "# vars: x, y, z\n# coefficients: rationals\ny*z\n");
}

#[test]
fn normal_form_takes_the_ring_from_the_source_that_names_it() {
    let run = sylv(&[
        "normal-form",
        "--basis",
        &input("cyclic3-basis.syl"),
        "--poly-file",
        &input("poly.text"),
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(
        run.out, "# vars: x, y, z\n# modulus: 32003\ny*z\n",
        "the syl basis takes the names and the modulus of the text polynomial"
    );
}

#[test]
fn two_sources_must_name_one_ring() {
    let run = sylv(&[
        "normal-form",
        "--basis",
        &input("cyclic3-basis.text"),
        "--poly-file",
        &input("poly-q.text"),
    ]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("names"), "{}", run.err);
}

#[test]
fn a_named_ring_rules_out_the_ring_options() {
    let run = sylv(&["gb", &input("cyclic3.text"), "--modulus", "7"]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("do not apply"), "{}", run.err);
}

#[test]
fn a_syl_source_needs_a_named_ring() {
    let run = sylv(&["gb", &input("cyclic3.syl")]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("--modulus"), "{}", run.err);
}

#[test]
fn only_one_source_reads_standard_input() {
    let run = sylv_stdin(
        &["normal-form", "--basis", "-", "--poly-file", "-"],
        CYCLIC3,
    );
    assert_eq!(run.code, 2);
    assert!(run.err.contains("standard input"), "{}", run.err);
}

#[test]
fn a_polynomial_file_holds_one_polynomial() {
    let run = sylv(&[
        "normal-form",
        "--basis",
        &input("cyclic3-basis.text"),
        "--poly-file",
        &input("cyclic3.text"),
    ]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("exactly one"), "{}", run.err);
}

#[test]
fn member_reports_a_member() {
    let run = sylv(&[
        "member",
        "--basis",
        &input("cyclic3-basis.text"),
        "--poly",
        "x*y*z - 1",
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, "true\n");
}

#[test]
fn member_reports_a_polynomial_outside_the_ideal() {
    let run = sylv(&[
        "member",
        "--basis",
        &input("cyclic3-basis.text"),
        "--poly",
        "x",
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, "false\n");
}

#[test]
fn member_reports_a_member_over_the_rationals() {
    let run = sylv(&[
        "member",
        "--basis",
        &input("cyclic3-basis-q.text"),
        "--poly",
        "x*y*z - 1",
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, "true\n");
}

#[test]
fn member_rejects_a_list_that_is_not_a_basis() {
    let run = sylv(&["member", "--basis", &input("cyclic3.text"), "--poly", "x"]);
    assert_eq!(run.code, 2, "a supplied basis that fails a check is input");
}

#[test]
fn hilbert_prints_the_series() {
    let run = sylv(&["hilbert", &input("cyclic3.text")]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(
        run.out,
        "\
series: 1, -1, -1, 0, 1, 1, -1
denominator_power: 3
dimension: 0
multiplicity: 6
"
    );
}

#[test]
fn hilbert_reads_the_same_series_over_the_rationals() {
    let run = sylv(&["hilbert", &input("cyclic3-q.text")]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert!(
        run.out.contains("series: 1, -1, -1, 0, 1, 1, -1\n"),
        "{}",
        run.out
    );
}

/// A value the rational driver produced carries what the run
/// established. A basis the command reads is the caller's own, so it
/// carries no such line.
#[test]
fn a_rational_hilbert_value_states_what_the_run_established() {
    let default = sylv(&["hilbert", &input("cyclic3-q.text")]);
    assert_eq!(default.code, 0, "{}", default.err);
    assert!(
        default.out.starts_with("# established: contains-input\n"),
        "{}",
        default.out
    );

    let unchanged = sylv(&["dim", &input("cyclic3-q.text"), "--stop", "unchanged"]);
    assert_eq!(unchanged.code, 0, "{}", unchanged.err);
    assert!(
        unchanged.out.starts_with("# established: unchanged\n"),
        "{}",
        unchanged.out
    );

    let from_basis = sylv(&["dim", &input("cyclic3-basis-q.text"), "--from-basis"]);
    assert_eq!(from_basis.code, 0, "{}", from_basis.err);
    assert_eq!(from_basis.out, "dimension: 0\n");
}

#[test]
fn hilbert_reads_a_basis() {
    let run = sylv(&["hilbert", &input("cyclic3-basis.text"), "--from-basis"]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert!(run.out.contains("dimension: 0"), "{}", run.out);
}

#[test]
fn hilbert_from_a_basis_runs_no_computation() {
    let run = sylv(&[
        "hilbert",
        &input("cyclic3-basis.text"),
        "--from-basis",
        "--backend",
        "classic",
    ]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("--backend"), "{}", run.err);
}

/// A basis whose colon chain is 65,535 nodes deep reports its memory
/// limit. It does not abort: the worklist of `src/hilbert.rs` runs on the
/// heap.
#[test]
fn hilbert_of_a_deep_colon_chain_reports_its_limit() {
    let run = sylv(&[
        "hilbert",
        &input("deep-chain.text"),
        "--from-basis",
        "--memory",
        "67108864",
    ]);
    assert_eq!(run.code, 3, "{}", run.err);
    assert!(run.err.contains("memory limit"), "{}", run.err);
}

#[test]
fn dim_prints_the_dimension() {
    let run = sylv(&["dim", &input("cyclic3.text")]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, "dimension: 0\n");
}

#[test]
fn dim_of_the_unit_ideal_is_none() {
    let run = sylv(&["dim", &input("unit.text")]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, "dimension: none\n");
}

#[test]
fn hilbert_of_the_unit_ideal_states_the_zero_numerator() {
    let run = sylv(&["hilbert", &input("unit.text")]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(
        run.out,
        "series: 0\ndenominator_power: 2\ndimension: none\nmultiplicity: none\n"
    );
}

#[test]
fn a_polynomial_that_does_not_parse_is_input() {
    let run = sylv(&["gb", &input("broken.text")]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("unexpected character"), "{}", run.err);
}

#[test]
fn an_unknown_extension_is_a_usage_error() {
    let path = scratch("system.unknown");
    fs::write(&path, CYCLIC3).expect("the input file");
    let run = sylv(&["gb", &text_of(&path)]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("format"), "{}", run.err);
}

#[test]
fn a_deadline_that_is_out_stops_the_command() {
    let run = sylv(&[
        "gb",
        &input("cyclic4.syl"),
        "--modulus",
        "32003",
        "--timeout",
        "0",
    ]);
    assert_eq!(run.code, 3);
    assert!(run.err.contains("deadline"), "{}", run.err);
}

#[test]
fn a_memory_limit_that_is_out_stops_the_command() {
    let run = sylv(&[
        "gb",
        &input("cyclic4.syl"),
        "--modulus",
        "32003",
        "--memory",
        "1",
    ]);
    assert_eq!(run.code, 3);
}

/// A timeout no duration holds is a usage error, not a panic inside the
/// conversion.
#[test]
fn a_timeout_past_what_a_duration_holds_is_a_usage_error() {
    let run = sylv(&["gb", &input("cyclic3.text"), "--timeout", "1e300"]);
    assert_eq!(run.code, 2, "{}", run.err);
    assert!(run.err.contains("--timeout"), "{}", run.err);
}

/// A `syl` header naming more variables than a ring holds is a usage
/// error. The count is checked before the reader reserves for it.
#[test]
fn a_syl_variable_count_past_the_ring_is_a_usage_error() {
    let path = scratch("wide.syl");
    fs::write(
        &path,
        "18446744073709551615
1,0
",
    )
    .expect("the input file");
    let run = sylv(&["gb", &text_of(&path), "--modulus", "32003"]);
    assert_eq!(run.code, 2, "{}", run.err);
    assert!(run.err.contains("largest supported count"), "{}", run.err);
}

#[test]
fn a_rational_option_on_a_prime_field_is_a_usage_error() {
    let run = sylv(&["gb", &input("cyclic3.text"), "--stop", "contains-input"]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("--stop"), "{}", run.err);
}

#[test]
fn the_rational_stopping_rule_runs() {
    let run = sylv(&[
        "gb",
        &input("cyclic3-q.text"),
        "--stop",
        "contains-input",
        "--extra-primes",
        "1",
        "--report",
    ]);
    assert_eq!(run.code, 0, "{}", run.err);
    assert_eq!(run.out, CYCLIC3_Q);
    assert!(
        run.err.contains("established: contains-input"),
        "{}",
        run.err
    );
}

#[test]
fn a_file_that_cannot_be_read_is_an_input_failure() {
    let run = sylv(&["gb", &input("absent.text")]);
    assert_eq!(run.code, 5);
    assert!(run.err.contains("cannot read"), "{}", run.err);
}

#[test]
fn a_file_that_cannot_be_written_is_an_output_failure() {
    let run = sylv(&[
        "gb",
        &input("cyclic3.text"),
        "-o",
        "/nonexistent-directory/basis.text",
    ]);
    assert_eq!(run.code, 5);
    assert!(run.err.contains("cannot write"), "{}", run.err);
}

#[test]
fn help_names_every_subcommand() {
    let run = sylv(&["--help"]);
    assert_eq!(run.code, 0, "{}", run.err);
    for name in [
        "gb",
        "certify",
        "normal-form",
        "member",
        "hilbert",
        "dim",
        "verify",
    ] {
        assert!(run.out.contains(name), "--help names {name}: {}", run.out);
    }
}

#[test]
fn a_division_command_names_its_polynomial() {
    let run = sylv(&["member", "--basis", &input("cyclic3-basis.text")]);
    assert_eq!(run.code, 2);
    assert!(run.err.contains("--poly"), "{}", run.err);
}
