# Protocol

The measured clean worktree is commit `396aa0d6cc581078ce9889c78317264501cc8eb9`.
The public v0.3 baseline is commit `5a32b213cb83e509887f38a7b389fe1bc86d9edf`.
The patched worktree has the v2 verifier changes and the measurement test.
The engine and verifier sources at both baseline commits are identical;
`baseline-source-equivalence.txt` records the Git tree hashes and comparison.
The patched source files and the `tests/verifier_raw.rs` raw harness logic
are archived in `patched-source.tar`. The portable harness is also archived
as `verifier_raw.rs`. Its only post-run normalization replaces the local
certificate directory prefix with the relative `baseline` directory. The
pre-normalization source hash is in `measured-harness-sha256.txt`; both timed
worktrees used the same source. `source-current-comparison.csv` compares each
patched source hash with the current worktree. The baseline equivalence
evidence is in `baseline-source-equivalence.txt`.

The compiler build command ran in each worktree:

```text
(cd <clean-worktree> && taskset -c 16-31 cargo +1.92 build --release -p sylvester-cli --bin sylv)
(cd <clean-worktree> && taskset -c 16-31 cargo +1.92 test --release --test verifier_raw --no-run)
(cd <patched-worktree> && taskset -c 16-31 cargo +1.92 build --release -p sylvester-cli --bin sylv)
(cd <patched-worktree> && taskset -c 16-31 cargo +1.92 test --release --test verifier_raw --no-run)
```

Each timing command ran while holding the shared benchmark lock:

```text
(cd <record-scratch> && taskset -c 0-3,12-15 <clean-worktree>/target/release/deps/verifier_raw-* --ignored --nocapture)
(cd <record-scratch> && taskset -c 0-3,12-15 <patched-worktree>/target/release/deps/verifier_raw-* --ignored --nocapture)
```

The raw harness in `verifier_raw.rs` reads each certificate from the
`baseline` directory below its working directory. The recorded runs used
`<record-scratch>/baseline` as that directory. The path normalization does not
change the timed verification code. The harness uses
`Limits { max_work_units: u64::MAX, ..Limits::default() }`, and measures only
the `verify_with_limits` call with `std::time::Instant`. It performs 21
repeats for each cell. The five cells use the same baseline certificate bytes
for both verifier builds. Their generated patched copies were compared byte
for byte before timing.

Certificate generation used the clean worktree CLI and these arguments for
each `<repo>/benchmarks/gb-comparison/inputs/<cell>.syl`. The output
directories are `<record-scratch>/baseline` and `<record-scratch>/candidate`:

```text
taskset -c 0-3,12-15 <clean-worktree>/target/release/sylv certify <repo>/benchmarks/gb-comparison/inputs/<cell>.syl --modulus 1073741827 --threads 1 --certificate <record-scratch>/baseline/<cell>.cert
taskset -c 0-3,12-15 <patched-worktree>/target/release/sylv certify <repo>/benchmarks/gb-comparison/inputs/<cell>.syl --modulus 1073741827 --threads 1 --certificate <record-scratch>/candidate/<cell>.cert
```

The source, compiler, CPU, commands, certificate hashes, raw values, and
computed medians are recorded beside this file. `SHA256SUMS` covers every
record file except itself.
