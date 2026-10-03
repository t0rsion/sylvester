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
taskset -c 16-31 cargo deny --workspace check
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
environments. The release workflow repeats the wheel smoke test on every
supported operating system.

Pull requests and manual runs execute the release builds and checks without
publishing. A pushed version tag runs those checks, then publishes the crates
and creates the GitHub release. The tagged commit must be an ancestor of `main`.

The release workflow emits these GitHub release assets. `<tag>` includes the
leading `v`, and `<version>` does not:

- `sylv-<tag>-x86_64-unknown-linux-gnu.tar.gz`
- `sylv-<tag>-aarch64-unknown-linux-gnu.tar.gz`
- `sylv-<tag>-x86_64-apple-darwin.tar.gz`
- `sylv-<tag>-aarch64-apple-darwin.tar.gz`
- `sylv-<tag>-x86_64-pc-windows-msvc.zip`
- `sylv-<tag>-aarch64-pc-windows-msvc.zip`
- `sylvester-<version>-*.whl`
- `sylvester-<version>.tar.gz`
- `sylvester-<version>.crate` and `sylvester-cli-<version>.crate`
- `SHA256SUMS`

Each standalone archive contains `sylv`, the CLI README, and both license
files. Extract the archive for the host target, then run `sylv --help`.
Verify downloaded assets with `sha256sum -c --ignore-missing SHA256SUMS` on
Linux or `shasum -a 256 -c SHA256SUMS` on macOS. Use `Get-FileHash` with the
`SHA256` algorithm on Windows, for example,
`Get-FileHash path/to/file -Algorithm SHA256`.

CLI and wheel release builds run `scripts/release-paths.py`. It carries
existing Rust flags through `CARGO_ENCODED_RUSTFLAGS`, then remaps the
workspace, Cargo home, and build home prefixes. Each build scans its binary
for Unix and Windows home paths, including UTF-16 paths.

The wheel smoke test creates a new virtual environment and installs the wheel
with no index and no dependencies. It then compares `sylvester.__version__`
with the installed distribution metadata and computes a one variable basis.
The native Windows ARM runner uses Python 3.11 because the runner provides no
Python 3.10 ARM64 build. The extension uses the CPython `abi3` for Python 3.10
and later, so the wheel remains compatible with its declared minimum.
Repeat that check locally with:

```sh
wheel_dir=$(mktemp -d)
python -m venv "$wheel_dir"
"$wheel_dir/bin/python" -m pip install --no-index --no-deps \
  path/to/sylvester-<version>-<wheel-tags>.whl
cd "$wheel_dir"
"$wheel_dir/bin/python" -c 'import sylvester; print(sylvester.__version__)'
```

The source distribution needs a Rust toolchain during installation because its
build backend compiles the extension. Install it in a separate environment
before testing the source archive. The `sympy` optional extra declares the
adapter dependency for source installs. From the source tree, install it with
`python -m pip install ".[sympy]"`. A wheel downloaded from a GitHub release
can install `sympy` separately when the adapter is used.

Store benchmark commands through `portable_command`. Track the compressed
JSON record, its digest, and the CSV report. Do not track the raw JSON.
Scan tracked text and decompressed records for Unix and Windows home paths,
secrets, and private contact data. Check the explicit private document names
`AGENTS.md`, `CLAUDE.md`, and `PLAN.md`. Names and addresses in license or
citation files must be intentional public attribution.

Run `cargo deny --workspace check` for the library and both interfaces.
The pinned Python dependency has two advisory exceptions in `deny.toml`.
The binding does not call `PyString::from_object` or
`PyCFunction::new_closure`, the affected APIs. Reassess these exceptions
when changing string conversion, callable registration, or the dependency pin.

Do not publish from a development branch. Cut the public release branch so
removed private data is absent from its history. Tag the tested commit.
