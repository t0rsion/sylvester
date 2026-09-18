# Rational equality exploratory probe, 2026-09-18

`results.json.zst`, `REPORT.md`, and `results.csv` are generated from the
captured output of four exact rational equality checks. The JSON includes the
measured source commit, source and input hashes, portable commands, limits,
toolchain, CPU affinity, and sanitized command output. No source archive is
included because the recorded commit and source hashes identify the measured
code in the release development history.

Each input ran once under the typed 300 second and 8 GiB checker limits. The
probe is exploratory and has no performance gate. The P-core benchmark lock
serialized the four commands with the release benchmark. The lock does not
provide exclusive host CPU access, and an unrelated CPU-bound process with
affinity spanning the sampled CPUs was observed during the run.

Inspect the JSON with `unzstd -c results.json.zst | python -m json.tool` and
verify the tracked files with `sha256sum -c SHA256SUMS`.

Future runs must set `SYLVESTER_PCORE_BENCH_LOCK` to the shared lock file and
hold it while running the probe:

```text
flock -n "$SYLVESTER_PCORE_BENCH_LOCK" -c 'RATIONAL_EQUALITY_ONLY=<name> CARGO_TARGET_DIR=<isolated-target> RUST_TEST_THREADS=1 taskset -c 0-3,12-15 cargo +1.92 test --release --test rational_equality -- --ignored --nocapture'
```
