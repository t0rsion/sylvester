# sylvester-py

Python bindings for the `sylvester` crate: reduced Gröbner bases under
grevlex, over a prime field or over the rational numbers, with independent
certificate verifiers. The Python package is named `sylvester`.

Release tags build wheels and a source distribution. They are GitHub
release artifacts, not PyPI releases. `publish = false` guards the Rust
crate, and the private classifier blocks a PyPI upload.

## Build and test

Build into the active virtualenv, then run the tests:

```sh
cd crates/sylvester-py
maturin develop --release
python -m pytest tests
```

`maturin develop` needs [maturin](https://www.maturin.rs/) 1.11 or later,
Python 3.10 or later, and the pinned Rust toolchain of
`rust-toolchain.toml`. The wheel is abi3, so one build runs on every
CPython from 3.10 up. To build a wheel instead of installing in place, run
`maturin build --release`.

## The surface

The coefficient domain is a type parameter in Rust and a runtime choice
here, so one Python class covers both domains and the ring constructor
fixes which one.

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

Over the rational numbers the coefficients cross as `fractions.Fraction`:

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

`PolynomialRing` builds every value: `prime_field(p, variables)`,
`rationals(variables)`, `parse(text)`, `polynomial(terms)`, and
`ideal(generators)`. `Ideal` computes: `groebner_basis(**options)`,
`groebner_basis_with_report(**options)`, and
`groebner_basis_certified(**options)`. `GroebnerBasis` is a sequence and
answers questions about the ideal: `normal_form(f)`, `contains(f)`,
`hilbert_series()`, `krull_dimension()`, `is_homogeneous()`, `lift()`, and
the static `from_polynomials(ring, polynomials)`. `sylvester.verify(data)`
checks certificate bytes on their own.

The compute methods take keywords, each defaulting to `None`: `backend`
(`"f4"` or `"classic"`), `timeout` in seconds, `memory_limit` in bytes,
`threads`, and, over the rational numbers alone, `stop` (`"unchanged"` or
`"contains_input"`) and `extra_primes`. Omitted values use the Rust defaults,
except that `stop=None` selects `"contains_input"`. Passing a rational
keyword to a prime-field ideal raises `ValueError`.

Under the sequence protocol `f in basis` reads as "is `f` one of these
polynomials", which is not ideal membership. Membership is
`basis.contains(f)`. `__contains__` is not implemented, because a wrong
answer there would be silent.

## What the results claim

`groebner_basis` returns the engine's unverified result.
`groebner_basis_certified` returns a value only after an independent
verifier accepts the certificate bytes. The basis is decoded from those
bytes.

There is no certified path over the rational numbers: the vote, the Chinese
remainder step, and the rational reconstruction are outside every
certificate. `groebner_basis_certified` on a rational ideal raises
`ValueError`. `GroebnerBasis.lift()` carries the counters of the run and
`established`, which names what the stopping rule observed. `"unchanged"`
means the lift did not change over the confirming primes and establishes
nothing about the ideal. `"contains_input"` adds two exact tests over `Q`,
after which the basis is the reduced Gröbner basis of an ideal that
contains the input ideal. That is not equality: the basis `{1}` passes both
tests.

## Errors

Every exception subclasses `sylvester.SylvesterError`. Bad input also
subclasses `ValueError`, and an exhausted budget, a structural limit, or a
defect also subclasses `RuntimeError`.

| class | raised for |
|---|---|
| `RingError` | a ring or a polynomial the domain rejects |
| `ParseError` | text the ring cannot read |
| `BasisError` | a supplied list that is not a reduced Gröbner basis |
| `CertificateInvalid` | certificate bytes the verifier rejected |
| `Timeout`, `MemoryLimitExceeded` | an exhausted budget, both `BudgetExhausted` |
| `LimitExceeded` | a degree, an exponent, a table, a prime sequence, or a certificate cap |
| `InternalDefect` | the crate contradicted itself |

## Threads

Every computation releases the GIL around the Rust work, so other Python
threads keep running while it runs. A Ctrl-C during a call takes effect
when the call returns: the engines have no cancellation point a signal
reaches. Use the `timeout` keyword to bound a call instead.
