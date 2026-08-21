# sylvester

A signature-based Gröbner basis library over finite prime fields. Every
value comes from a `PolynomialRing`, which fixes the prime (`2 <= p <= 2^31
- 1`), the variable names (at most 256), and the variable order. The
monomial order is grevlex, and there is no other. Two backends compute the
reduced basis: classic F5, one critical pair at a time in signature order,
and a Macaulay-matrix backend that batches pairs by degree. The certified
path returns a basis only after an independent verifier accepts a
certificate for it, and that verifier shares no code with the engines.

This is a first alpha. The API will change.

## Status

Repaired, checked against oracles, not proven. Two counterexamples proved
both backends wrong as first extracted; this tree carries the repair,
checked by the independent checker in `tests/known_defects.rs` and a
randomized differential suite against a self-contained Buchberger oracle.
Termination has a pen-and-paper proof by Dickson's lemma. No proof here is
machine-checked, and the verifier itself is a trusted, unproven base.

Certification is classic-only: the matrix backend emits no certificate, so
a basis from it is unchecked. This is a `0.1` alpha; treat it as such.

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
    println!("verifier accepted {} bytes", certified.certificate().len());
    Ok(())
}
```

`ComputeOptions` picks the `backend`, `timeout`, and `memory_limit`. An
exhausted budget is a typed `ComputeError::Timeout` or
`MemoryLimitExceeded`; a monomial past total degree 65,535 is
`ComputeError::DegreeLimit`. Nothing truncates silently.

## Certificates

`groebner_basis_certified` returns the basis together with the bytes the
verifier accepted, decoded from those bytes rather than taken from the
engine. Acceptance establishes that the input and the returned basis
generate the same ideal and that the basis is the reduced Gröbner basis of
that ideal, over F_p under grevlex. It does not make an engine correct: an
engine can still hang, exhaust its budget, or write a certificate the
verifier rejects.

The bytes are canonical JSON under the `sylv-gb-cert-v1` contract, so
`verify::verify` can check a stored certificate later, applying every
resource cap before it allocates. A certified run costs several times an
uncertified one.

## Performance

sylvester is slower than msolve and Groebner.jl and does not finish every
instance they do. Measured comparisons come with v0.2.

## Roadmap

- v0.2: speed. An F4 engine, with certificates built from an execution trace.
- v0.3: Python bindings, rational coefficients, and more ideal operations.

## Building and testing

Rust 1.85 or later; `rust-toolchain.toml` pins 1.92 for development.

```
cargo +1.92 test
cargo +1.92 test --features parallel
cargo +1.92 test --release --test differential -- --ignored    # slow sweeps
```

`parallel` is the only feature: it adds rayon row elimination to the matrix
backend.

## License

MIT OR Apache-2.0, at your option. See [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).
