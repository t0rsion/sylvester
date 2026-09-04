# sylvester

Gröbner bases over prime fields and the rationals: an F4 engine, a
multimodular rational driver, a classic F5 oracle, and independent
certificate verifiers, with a command line interface and Python bindings.
Every value comes from a `PolynomialRing<D>`, which fixes the coefficient
domain, the variable names, and the variable order. `D` is `PrimeField`
(the default) or `Rationals`. The monomial order is grevlex over the
variable order, sealed and named `grevlex-v1` in every certificate. There
is no user-supplied comparator.

Over `F_p`, two backends compute the reduced basis. F4 is the default: it
batches critical pairs by degree, builds one sparse matrix per batch over
interned monomials, and reduces it, with Gebauer-Moller pair management.
Classic F5 takes critical pairs one at a time in signature order. It is
much slower than F4, and it stays in the crate as the certified
`sylv-gb-cert-v1` oracle and as the reference the differential test suite
checks F4 against. A certified path writes a certificate for the result
and returns a basis only after an independent verifier accepts the
certificate bytes. That verifier shares no code with the engines.

Over `Q`, one multimodular driver clears denominators, runs F4 or classic
F5 modulo many 31-bit primes, combines the runs by the Chinese remainder
theorem, and lifts the result to rational numbers. This path is a
heuristic: no isolated verifier checks it, and there is no certified path
over `Q`.

The name is for James Joseph Sylvester, who coined the term "syzygy".

## Status

F4 over `F_p`: checked against the classic backend, against a
self-contained Buchberger oracle, and against msolve on the comparison
record under `benchmarks/gb-comparison`. Not proven.

Classic over `F_p`: repaired, empirically certified, not proven. Two
counterexamples proved both original backends wrong as extracted, and
this tree carries the repair. `KNOWN_ISSUES.md` records the systems, the
wrong outputs, the witnesses, and the root causes. The counterexample
tests in `tests/known_defects.rs` check the correct bases through an
independent checker, and a randomized differential suite compares every
backend against a self-contained Buchberger oracle. Termination of the
classic backend has a pen-and-paper proof by Dickson's lemma. No proof
here is machine-checked.

Over `Q`: heuristic, not certified, and not proven. `Ideal::<Rationals>::groebner_basis`
returns a basis after a stopping rule observes that the lift stopped
changing over further primes; `GroebnerBasis::lift` names what the rule
established through `Established`. `Established::Unchanged`, the Rust
default outcome, establishes nothing about the ideal: it says the reconstructed
basis agreed with itself over the confirming primes. `Established::ContainsInput`
adds two exact tests over `Q` and establishes that the returned basis is
the reduced Gröbner basis of an ideal that contains the input ideal, not
that it equals the input ideal, and the basis `{1}` passes both tests. There
is no path to a stronger claim over `Q` in this release; section 4 of the
[rational design](docs/rational-design.md) states why a per-prime certificate
would not close the gap.

The `sylvester` library and `sylvester-cli` binary are on crates.io. The
Python distribution remains a release artifact and is not uploaded to PyPI.

## Example

```rust
use sylvester::{Budget, ComputeOptions, PolynomialRing, RationalOptions};

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
    assert_eq!(certified.basis().len(), basis.len());

    let f = ring.parse_polynomial("x^2*y + x*y*z")?;
    let remainder = basis.normal_form(&f, Budget::new())?;
    println!("{f} reduces to {remainder}");

    let dimension = basis.krull_dimension(Budget::new())?;
    println!("krull dimension {dimension:?}");

    // Over Q the basis is a heuristic: the multimodular driver computed it,
    // and `lift()` carries what its stopping rule observed.
    let qring = PolynomialRing::rationals(["x", "y", "z"])?;
    let qideal = qring.ideal([
        qring.parse_polynomial("x + y + z")?,
        qring.parse_polynomial("x*y + y*z + z*x")?,
        qring.parse_polynomial("x*y*z - 1")?,
    ])?;
    let qbasis = qideal.groebner_basis(RationalOptions::new())?;
    println!("{}", qbasis[0]);
    println!("{:?}", qbasis.lift().map(|lift| lift.established));

    Ok(())
}
```

`groebner_basis` returns the engine's unverified claim.
`groebner_basis_with_report` returns the same basis next to a
`ComputeReport`. `groebner_basis_certified`, defined on `Ideal<PrimeField>`
alone, returns a basis only after an independent verifier accepts the
certificate bytes the run wrote. `GroebnerBasis::normal_form`,
`GroebnerBasis::contains`, `GroebnerBasis::hilbert_series`, and
`GroebnerBasis::krull_dimension` work in both domains, under a `Budget`.
Each of them is about the ideal the basis generates. Over `Q` what
relates that ideal to the input is what `Established` says: under
`Unchanged` nothing does, and under `ContainsInput` a `contains` answer
of `false` carries over to the input ideal while `true` does not.
`GroebnerBasis::from_polynomials` builds a basis from a list the caller
already holds, checked rather than computed. Runnable versions of the
prime-field parts of this code are in `examples/`.

## Command line

`sylv`, in the `sylvester-cli` crate, computes, certifies, verifies, and
reduces from the shell:

```
sylv gb cyclic-7.ms
sylv gb system.syl --modulus 1073741827 --backend classic --threads 8
sylv gb system.text --rationals --stop contains-input --extra-primes 3
sylv certify system.ms --certificate system.cert
sylv verify system.cert
sylv normal-form --basis basis.text --poly "x^2*y - 1"
sylv member --basis basis.text --poly-file f.text
sylv hilbert system.ms
sylv dim system.ms --from-basis
```

`hilbert` and `dim` over `Q` write an `# established:` comment above the
value, because the value is of the ideal the lifted basis generates.

It reads `ms` (msolve's format), `syl` (the benchmark format, `.sylq` for
a rational instance), and `text` (the crate's own syntax), resolving one
ring per command over every file it reads. `crates/sylvester-cli/README.md`
documents every subcommand, the ring resolution rules, and the six exit
codes.

## Python

`sylvester`, in the `sylvester-py` crate, binds the same library through
pyo3:

```python
import sylvester

ring = sylvester.PolynomialRing.prime_field(32003, ["x", "y", "z"])
ideal = ring.ideal([
    ring.parse("x^2*y - 3*z + 1"),
    ring.parse("x*z - y"),
])

basis = ideal.groebner_basis()
print([str(f) for f in basis])
print(basis.krull_dimension())

certified = ideal.groebner_basis_certified()
verified = sylvester.verify(certified.certificate)
assert verified.modulus == 32003
```

Over the rational numbers the coefficients cross as `fractions.Fraction`,
and `basis.lift()` carries the counters of the multimodular run:

```python
ring = sylvester.PolynomialRing.rationals(["x", "y", "z"])
ideal = ring.ideal([
    ring.parse("x + y + z"),
    ring.parse("x*y + y*z + z*x"),
    ring.parse("x*y*z - 1"),
])
basis = ideal.groebner_basis()
assert str(basis[0]) == "z^3 - 1"
assert basis.lift()["established"] == "contains_input"
```

`crates/sylvester-py/README.md` documents the full surface and the error
classes.

## Options and errors

`ComputeOptions` builds the backend, the budget, and the thread count of
one call over `F_p`: `backend` (`Backend::F4`, the default, or
`Backend::Classic`), a `Budget` (`timeout` and `memory_limit`, taken
together or set through the `ComputeOptions::timeout` and
`ComputeOptions::memory_limit` shortcuts), and `threads`. `threads(0)` is
the rayon global pool; `threads(1)` keeps the computation on the calling
thread; a count above 1 runs it in a rayon pool of that size. The classic
backend, and every run that writes a certificate, computes on one thread
whatever the count says.

`RationalOptions` wraps a `ComputeOptions`, whose `threads` sets the
number of prime runs the multimodular driver starts at once (each prime
run itself uses one thread), and adds `stop`: `RationalStop::Unchanged {
extra }` (the Rust default, `extra` = 2) or `RationalStop::ContainsInput {
extra }`, the count of further confirming primes each rule needs. The CLI
and Python binding default to `ContainsInput`.

`Budget` alone, with no backend or thread count, drives `normal_form`,
`contains`, `hilbert_series`, `krull_dimension`, and `from_polynomials`,
in both domains.

An exhausted budget is a typed error: `Timeout` or
`MemoryLimitExceeded`, on every operation that takes one. Nothing
truncates silently. `ComputeError::ExponentLimit` (F4: one exponent past
65,535) and `ComputeError::DegreeLimit` (classic: a critical pair's
leading-monomial least common multiple past 65,535 in total degree) report
an input past what this release supports. `ComputeError::TableFull`
reports an F4 monomial table that interned more entries than one table
holds. `ComputeError::PrimesExhausted` reports a rational run that walked
the whole descending prime sequence without stopping. `NormalFormError`
and `BasisError` add `ExponentLimit` for the same reason division can
overflow an exponent that a critical pair never would. `Felt`, the
coefficient a prime-field term carries, is public; `Felt::value` reads the
residue back as a `u64`.

## Certificates

The certificate format follows the backend, over `F_p` alone. The classic
backend writes `sylv-gb-cert-v1`, canonical JSON built from the origin
cofactors it tracks through reduction, specified in
`docs/certificate-v1.md`. The F4 backend writes `sylv-gb-cert-v2`, a
canonical binary trace of the operations the run performed, specified in
`docs/certificate-v2.md`. There is no format option; the backend decides,
and there is no certificate over `Q` in this release.
`groebner_basis_certified` is a method of `Ideal<PrimeField>` alone, so
calling it over `Q` is a compile error in Rust, a usage error (exit 2) in
`sylv`, and a `ValueError` in Python.

Both formats are checked by isolated verifiers that share no code with
the engines. `verify::verify` (`sylv verify`, `sylvester.verify` in
Python) dispatches on the certificate's first byte. Acceptance proves two
facts over `F_p` under grevlex: the input and the returned basis generate
the same ideal, and the basis is the reduced Gröbner basis of that ideal.
Acceptance excludes a wrong answer. It does not make an engine correct:
an engine can still hang, exhaust its budget, or write a certificate the
verifier rejects.

`benchmarks/gb-comparison/REPORT.md` holds the current record: the
protocol, the timings against Singular, msolve, Macaulay2, and
Groebner.jl over `F_p`, against Singular, msolve, and Groebner.jl over
`Q`, the correctness cross-check in both domains, and the gate result
against msolve measured in-process. The gate is over `F_p` alone; the
rational cells are measured and recorded, not gated, because the
per-prime cost this release runs is not the pipeline the design intends.
The record is regenerated at the release commit; read it rather than any
number in this file.

`verify::verify` checks stored v1 or v2 bytes. `verify::verify_with_limits`
adds caller-supplied caps and a deadline.

Candidates for a later release, not commitments: general monomial orders,
elimination orders, FGLM, and the radical.

## Evidence

The default suite checks both engines against a criterion-free Buchberger
oracle over five prime fields. The release suite adds 166,375 exhaustive
systems over \(\mathbb F_2\), 1,000 degree-reversal systems over
\(\mathbb F_3\), and an `eco-8` F4-to-classic comparison.

The external harness compares complete bases with Singular, msolve,
Macaulay2, and Groebner.jl. Timings apply only to the recorded machine,
versions, inputs, and limits. The current one-thread F4 gate passes at a
1.00x geometric mean and a 1.09x worst cell against msolve. See
[the comparison record][comparison-record].

## Build

The manifest states the minimum Rust version. `rust-toolchain.toml` pins the
development toolchain.

```text
cargo test --workspace
cargo test --release --test differential -- --ignored
cargo test --release --test f4_engine -- --ignored
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo fmt --all --check
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
```

The workspace run covers the library, `sylvester-cli`'s golden tests, and
every doctest.
The full release checklist is in [docs/releasing.md](docs/releasing.md).

Python, from a virtual environment:

```
cd crates/sylvester-py
maturin develop --release
python -m pytest tests
```

## License

MIT OR Apache-2.0, at your option. See [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).

[certificate-v1]: docs/certificate-v1.md
[certificate-v2]: docs/certificate-v2.md
[comparison-record]: benchmarks/gb-comparison/REPORT.md
[known-issues]: KNOWN_ISSUES.md
