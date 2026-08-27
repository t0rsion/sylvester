# sylvester

Gröbner bases over finite prime fields in Rust.

Version 0.2.0 contains:

- a default F4 engine with sparse matrix reduction;
- a repaired classic F5 engine used as an independent oracle;
- independent verifiers for classic and F4 certificates;
- typed time, memory, exponent, and table limits;
- deterministic differential tests on complete bases.

This is an alpha. The API may change before 1.0.

## Status

Both engines are empirically checked, not proven. The default suite compares
complete reduced bases with an independent Buchberger oracle. Fixed
counterexamples block failures that output-only checks miss.

An accepted certificate establishes the result. The verifier is trusted,
unproven code. No proof in this repository is machine-checked.

See [KNOWN_ISSUES.md][known-issues] for the counterexamples, repairs,
termination arguments, and evidence limits.

## Domain

`PolynomialRing` fixes:

- a prime \(p\), where \(2\le p\le2^{31}-1\);
- at most 256 variable names;
- the variable order;
- grevlex, named `grevlex-v1` in certificates.

There is no runtime monomial order. One exponent is at most 65,535.

## Example

```rust
use sylvester::{ComputeOptions, PolynomialRing};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"])?;
    let ideal = ring.ideal([
        ring.parse_polynomial("x + y + z")?,
        ring.parse_polynomial("x*y + y*z + z*x")?,
        ring.parse_polynomial("x*y*z - 1")?,
    ])?;

    let basis = ideal.groebner_basis(ComputeOptions::new())?;
    for polynomial in &basis {
        println!("{polynomial}");
    }

    let certified = ideal.groebner_basis_certified(ComputeOptions::new())?;
    println!("{} certificate bytes", certified.certificate().len());
    Ok(())
}
```

F4 is the default. Select classic explicitly when needed:

```rust
use sylvester::{Backend, ComputeOptions};

let options = ComputeOptions::new()
    .backend(Backend::Classic)
    .memory_limit(512 << 20)
    .threads(1);
# let _ = options;
```

`ComputeOptions` sets the backend, timeout, memory limit, and F4 thread count.
Zero threads use Rayon's global pool. One thread uses the calling thread.
Classic and certified calls remain serial.

`Ideal::groebner_basis_with_report` also returns elapsed time and F4 counters.

## Partiality

The computation returns typed errors. It never truncates silently.

| Error | Meaning |
| --- | --- |
| `Timeout` | The deadline passed. |
| `MemoryLimitExceeded` | The tracked live data passed the byte budget. |
| `DegreeLimit` | Classic needed a monomial past the supported total degree. |
| `ExponentLimit` | F4 needed one exponent above 65,535. |
| `TableFull` | An F4 monomial table exhausted its identifier space. |

The memory budget charges owned engine and certificate data. It does not cap
process RSS.

## Certificates

`Ideal::groebner_basis_certified` follows the selected backend:

| Backend | Contract | Encoding |
| --- | --- | --- |
| Classic | `sylv-gb-cert-v1` | Canonical JSON with explicit cofactors |
| F4 | `sylv-gb-cert-v2` | Canonical binary operation trace |

The method verifies its bytes before returning. The returned basis is decoded
from those accepted bytes.

Acceptance establishes:

\[
\langle F\rangle=\langle G\rangle,
\qquad
G=\operatorname{rGB}_{\mathrm{grevlex}}(\langle F\rangle).
\]

`verify::verify` checks stored v1 or v2 bytes. `verify::verify_with_limits`
adds caller-supplied caps and a deadline. The verifier code is isolated from
both engines.

The frozen contracts are [certificate-v1.md][certificate-v1] and
[certificate-v2.md][certificate-v2].

## Evidence

The default suite checks both engines against a criterion-free Buchberger
oracle over five prime fields. The release suite adds 166,375 exhaustive
systems over \(\mathbb F_2\), 1,000 degree-reversal systems over
\(\mathbb F_3\), and an `eco-8` F4-to-classic comparison.

The external harness compares complete bases with Singular, msolve,
Macaulay2, and Groebner.jl. Timings apply only to the recorded machine,
versions, inputs, and limits. The v0.2 one-thread F4 gate passes at a
1.13x geometric mean and a 1.48x worst cell against msolve. See
[the comparison record][comparison-record].

## Build

The minimum supported Rust version is 1.88. Development pins Rust 1.92.

```text
cargo +1.92 test --workspace
cargo +1.92 test --release --test differential -- --ignored
cargo +1.92 test --release --test f4_engine -- --ignored
cargo +1.92 clippy --workspace --all-targets --all-features -- -D warnings
cargo +1.92 fmt --all --check
RUSTDOCFLAGS="-D warnings" cargo +1.92 doc --workspace --no-deps
```

The crate has no feature flags.

## License

MIT OR Apache-2.0, at your option. See [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).

[certificate-v1]: https://github.com/t0rsion/sylvester/blob/v0.2.0/docs/certificate-v1.md
[certificate-v2]: https://github.com/t0rsion/sylvester/blob/v0.2.0/docs/certificate-v2.md
[comparison-record]: https://github.com/t0rsion/sylvester/blob/v0.2.0/benchmarks/gb-comparison/REPORT.md
[known-issues]: https://github.com/t0rsion/sylvester/blob/v0.2.0/KNOWN_ISSUES.md
