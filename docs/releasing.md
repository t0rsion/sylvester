# Release checks

Run every check from the workspace root. Keep release work on the local
development branch until all checks pass.

## Structure

Measure cyclomatic complexity with Lizard 1.24.0:

```sh
uvx --from lizard==1.24.0 lizard -l rust -C 10 -w \
  src tests benches examples benchmarks/gb-comparison \
  crates/sylvester-cli/src crates/sylvester-cli/tests \
  crates/sylvester-py/src crates/sylvester-py/tests
```

The command must print no warning. Apply these limits to every function:

- CCN 1 through 5: leave it.
- CCN 6 through 10: refactor it when the same code changes.
- CCN 11 through 15: split it before release.
- CCN above 15: split it before any other release work.

Read the Lizard NLOC total before and after the release. Remove duplicate
paths, unused interfaces, aliases, dead code, and narrative comments. A
feature can add code. Its shared mechanisms must reduce the code needed per
interface.

Inspect each public boundary. A boundary has one input type, one error
contract, and one term for each concept. Before API stability, remove a weak boundary
instead of keeping a compatibility alias.

## Rust

```sh
taskset -c 16-31 cargo +1.92 fmt --all --check
taskset -c 16-31 cargo +1.92 clippy --workspace --all-targets --all-features -- -D warnings
taskset -c 16-31 cargo +1.92 test --workspace
taskset -c 16-31 env RUSTDOCFLAGS="-D warnings" cargo +1.92 doc --workspace --no-deps
taskset -c 16-31 cargo +1.88 test --workspace
taskset -c 16-31 cargo deny check
```

Run the slow oracle sweeps:

```sh
taskset -c 16-31 cargo +1.92 test --release --test differential -- --ignored
taskset -c 16-31 cargo +1.92 test --release --test f4_engine -- --ignored
```

## Python and benchmark tools

From an active virtual environment:

```sh
cd crates/sylvester-py
maturin develop --release
python -m pytest tests
```

Then run:

```sh
ruff check benchmarks/gb-comparison
ruff format --check benchmarks/gb-comparison
python -m unittest discover -s benchmarks/gb-comparison -p 'test_*.py'
```

## Packages and privacy

Build the library package and test its extracted source. Cargo cannot package
the CLI until the matching library version exists on crates.io. Test the CLI
in the workspace before publication, then package it after the library reaches
the index. Build and install the Python wheel and source distribution in fresh
environments.

Store benchmark commands through `portable_command`. Track the compressed
JSON record, its digest, and the CSV report. Do not track the raw JSON.
Scan tracked text and decompressed records for home directories, secrets,
and private contact data. Names and addresses in license or citation files
must be intentional public attribution.

Do not publish from a development branch. Cut the public release branch so
removed private data is absent from its history. Tag the tested commit.
