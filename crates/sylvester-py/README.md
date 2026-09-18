# sylvester-py

Python bindings for the `sylvester` crate. The package computes reduced
Gröbner bases under grevlex over a prime field or the rational numbers. It
also exposes finite quotient algebras, independent certificate verifiers,
and an optional structured SymPy adapter.

The monomial order is grevlex over the variable order passed to
the ring constructor. There is no Python option for another order.

## Install from a GitHub release

Download the wheel for the Python version and platform from the project's
[GitHub Releases](https://github.com/t0rsion/sylvester/releases). Replace
`WHEEL_FILE` with its filename:

```sh
python -m pip install WHEEL_FILE
```

The wheel requires CPython 3.10 or later and uses the CPython `abi3`.

To build from a checkout, create a virtual environment and run:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install "maturin>=1.12,<2" pytest
cd crates/sylvester-py
maturin develop --release
python -m pytest tests
```

In PowerShell, replace the activation command with `.venv\Scripts\Activate.ps1`.

`maturin build --release` writes a wheel instead of installing the extension
in the active environment.

## Quick start

`PolynomialRing.gens` returns the ring generators in variable order. Python
operators build new polynomials with the ring's coefficient arithmetic.

```python
import sylvester

ring = sylvester.PolynomialRing.prime_field(32003, ["x", "y", "z"])
x, y, z = ring.gens
f = x**2 * y - 3 * z + 1
ideal = ring.ideal([f, x * z - y])

basis = ideal.groebner_basis()
assert basis.contains(f)
print([str(polynomial) for polynomial in basis])

certified = ideal.groebner_basis_certified()
verified = sylvester.verify(certified.certificate)
assert verified.modulus == 32003
```

The same constructors work over the rational numbers:

```python
from fractions import Fraction

ring = sylvester.PolynomialRing.rationals(["x", "y"])
x, y = ring.gens
f = Fraction(1, 2) * x + y
assert f == ring.parse("1/2*x + y")
```

## Rings and polynomials

`PolynomialRing.prime_field(p, variables)` constructs a ring over `F_p`.
`p` must be prime. `PolynomialRing.rationals(variables)` constructs a ring
over `Q`. Each ring owns its variables, grevlex order, and coefficient domain.

The ring methods are:

- `variables`, `modulus`, and `gens` expose ring metadata. `modulus` is
  `None` over `Q`.
- `zero` and `one` return the constant polynomials of the ring.
- `parse(text)` reads the ring's polynomial syntax.
- `polynomial(terms)` accepts `(coefficient, exponents)` pairs. A coefficient
  is an `int`, a `fractions.Fraction`, or a `(numerator, denominator)` pair.
- `ideal(generators)` stores a list of polynomials from this ring.

`Polynomial.terms()` returns `(coefficient, exponents)` pairs with the largest
monomial first. `degree()` returns the total degree, or `None` for zero.
`is_zero()` tests the zero polynomial. Polynomials compare and hash by value.

`+`, `-`, `*`, unary `-`, and `**` work between compatible polynomials and
with scalar coefficients. These operators use the default unbounded budget.
Use `add`, `sub`, `mul`, or `pow` when the operation needs an explicit budget.

## Compute with explicit budgets

Compute methods accept keyword budgets. `timeout` is in seconds and
`memory_limit` is in bytes. Omitted values use the library defaults.

```python
basis = ideal.groebner_basis(
    backend="f4",
    timeout=2.0,
    memory_limit=256 << 20,
    threads=4,
)

bounded_sum = x.add(y, timeout=1.0, memory_limit=1 << 20)
```

`backend` is `"f4"` or `"classic"`. `threads` controls the prime-field run.
For a rational ideal, it controls the number of prime runs started at once.
Rational calls also accept `stop="unchanged"` or
`stop="contains_input"`, and `extra_primes` sets the number of confirming
primes. Passing rational-only keywords to a prime-field ideal raises
`ValueError`.

`Ideal.groebner_basis_with_report()` returns `(basis, report)`. The report
contains the backend, elapsed time, thread count, and backend counters.
`GroebnerBasis.from_polynomials(ring, polynomials)` checks a supplied reduced
Gröbner basis under the same `timeout` and `memory_limit` keywords.

## Bases, normal forms, and division

`GroebnerBasis` implements the sequence protocol. Indexing, slicing, length,
and iteration use the stored basis order. `basis.contains(f)` tests ideal
membership by reducing `f`. The expression `f in basis` scans the listed
polynomials and does not test ideal membership.

`basis.normal_form(f)` returns the remainder of division by the basis. It
keeps the remainder-only path and does not allocate quotient polynomials.
`basis.is_zero_dimensional()`, `basis.hilbert_series()`,
`basis.krull_dimension()`, and `basis.is_homogeneous()` inspect the basis.

`basis.divide(f)` returns `(quotients, remainder)`. If
`G = list(basis) = [g_0, ..., g_{k-1}]`, the result has one quotient for each
element of `G` and satisfies

```text
f = quotients[0] * g_0 + ... + quotients[k-1] * g_{k-1} + remainder
```

The quotient list follows the order of `G`, including when some quotients are
zero. The original generator list `F` used to compute that basis does not
occur in the identity.
`basis.divmod(f)` is an alias for `basis.divide(f)`. Both methods accept
`timeout` and `memory_limit`.

```python
quotients, remainder = basis.divide(f, timeout=2.0, memory_limit=256 << 20)
assert len(quotients) == len(basis)
rebuilt = remainder
for quotient, generator in zip(quotients, basis):
    rebuilt = rebuilt + quotient * generator
assert rebuilt == f
```

## Finite quotient algebras

`basis.finite_quotient()` builds the quotient algebra when the leading
monomial ideal contains a pure power of every variable. It raises
`ValueError` for a non-finite quotient. The returned
`FiniteQuotient.standard_monomials` lists exponent vectors in grevlex order,
and `dimension` gives their count.

```python
ring = sylvester.PolynomialRing.prime_field(101, ["x", "y"])
x, y = ring.gens
basis = ring.ideal([x**2, y**2]).groebner_basis()
quotient = basis.finite_quotient(timeout=2.0, memory_limit=16 << 20)

assert quotient.dimension == 4
residue = quotient.reduce(x**3 + x * y + 1)
assert residue == x * y + 1
assert len(quotient.coordinates(residue)) == quotient.dimension

matrix = quotient.multiplication_matrix(x)
assert matrix.dimension == quotient.dimension
assert str(quotient.characteristic_polynomial(x)) == "t^4"
assert str(quotient.minimal_polynomial(x)) == "t^2"
```

`FiniteQuotient.add`, `multiply`, and `pow` reduce quotient residues. Their
matrix and polynomial results are `MultiplicationMatrix` and
`UnivariatePolynomial` values. Matrix `entries` are nested row lists in
standard-monomial order, with ordinary zero coefficients. `flat_entries` is
row-major with the same zero values. `entry(row, column)` returns `None` only
for an out-of-range index. Univariate coefficients are constant-term first,
with ordinary zero coefficients. Every quotient operation accepts `timeout`
and `memory_limit`.

The quotient is defined by the basis supplied to `finite_quotient`. Over `Q`,
the rational lift determines what that basis establishes about the input
ideal.

## Rational results and provenance

The rational driver is multimodular and heuristic. It has no certified path.
`basis.lift()` returns the counters and the stopping claim when the basis came
from a rational computation. `"unchanged"` records a lift that stayed the
same over its confirming primes and establishes nothing about the input ideal.
`"contains_input"` records exact containment checks. It establishes that the
basis is the reduced Gröbner basis of an ideal containing the input ideal.
It does not establish equality with the input ideal.

`ideal.check_basis_equality(basis)` performs the exact rational check and
returns a `RationalEqualityCheck` with the input, basis, origin cofactors, and
metrics. `ResultEnvelope.from_checked_rational` records its `"equals_input"`
claim for later inspection.

`ResultEnvelope` stores input and basis expressions, ring metadata, an
untrusted provenance claim, and optional prime-field certificate bytes. It
provides `from_prime`, `from_certified`, `from_rational`,
`from_checked_rational`, `from_json`, and `to_json`.

```python
record = sylvester.ResultEnvelope.from_prime(
    ideal,
    basis,
    certified.certificate,
)
saved = record.to_json(timeout=1.0, memory_limit=1 << 20)
loaded = sylvester.ResultEnvelope.from_json(saved)
assert loaded.claimed_provenance == "unverified"
checked = loaded.verify_prime(max_bytes=1 << 20, max_work_units=1 << 20)
```

Loading a record keeps `claimed_provenance` as data. `verify_prime` verifies
the stored certificate under configurable verifier caps and checks that the
accepted input and basis match the record. Rational records have no
certificate contract. The record exposes `domain`, `variables`, `input`,
`basis`, `claimed_provenance`, and `certificate`.

## Certificates and verifier caps

`Ideal.groebner_basis_certified()` is available over a prime field. It returns
`CertifiedGroebnerBasis` only after an independent verifier accepts the
certificate. The value exposes `basis` and `certificate`.

`sylvester.verify(data)` accepts certificate bytes without using an engine.
The optional caps bound the verifier:

```python
verified = sylvester.verify(
    certified.certificate,
    timeout=2.0,
    max_bytes=16 << 20,
    max_work_units=1 << 20,
    max_intermediate_bytes=32 << 20,
    max_live_bytes=32 << 20,
)
```

The byte caps cover certificate input, intermediate data, and live verifier
data. `max_work_units` bounds verifier work. An exhausted cap raises
`LimitExceeded` or `Timeout` according to the cap. Certificate schemas store
variable positions rather than names, so `VerifiedGroebnerBasis` prints its
polynomials with synthetic names `x1`, `x2`, and so on.

## Threads and cancellation

Long-running Rust computations release the Python GIL. Other Python threads
can run while the call is active. Every budgeted method runs its Rust work in
a joined worker and enforces `timeout` and `memory_limit` when provided.

Calls originating on the main Python thread poll pending signals. Ctrl-C
requests cancellation, the binding joins the worker, and then raises
`KeyboardInterrupt`. `check_signals` is a no-op for calls originating on a
background Python thread, so Ctrl-C is not propagated through those calls.
Use `timeout` to bound a background call.

## Optional SymPy adapter

The adapter is separate from `sylvester` and does not import SymPy at package
import time. Install the optional dependency with `python -m pip install sympy`,
or install the `sympy` extra from a source checkout:

```sh
cd crates/sylvester-py
python -m pip install '.[sympy]'
```

`from_sympy` and `to_sympy` exchange structured `sympy.Poly` values:

```python
import sympy as sp
from sylvester.sympy import from_sympy, to_sympy

x_sym, y_sym = sp.symbols("x y")
source = sp.Poly(x_sym**2 + 2 * y_sym, x_sym, y_sym, modulus=101)
ring = sylvester.PolynomialRing.prime_field(101, ["x", "y"])
polynomial = from_sympy(source, ring)
assert to_sympy(polynomial) == source
```

The adapter accepts `QQ`, `ZZ`, and prime `GF(p)` domains. It preserves exact
coefficients. It rejects inexact coefficients, expressions, extensions, and
composite finite fields. SymPy generators must be plain symbols with names
and order matching the Sylvester ring. The adapter does not translate
monomial orders.

## Errors

Named library exceptions subclass `sylvester.SylvesterError`. Input errors
also subclass `ValueError`. Budget exhaustion and structural limits subclass
`RuntimeError`.

| class | raised for |
|---|---|
| `RingError` | a ring or a polynomial the domain rejects |
| `ParseError` | text the ring cannot read |
| `BasisError` | a supplied list that is not a reduced Gröbner basis |
| `CertificateInvalid` | certificate bytes the verifier rejects |
| `Timeout` | a deadline passed before the operation finished |
| `MemoryLimitExceeded` | live data passed `memory_limit` |
| `LimitExceeded` | a degree, exponent, monomial table, prime sequence, or verifier cap reached its bound |
| `InternalDefect` | the crate contradicted itself |

`Timeout` and `MemoryLimitExceeded` subclass `BudgetExhausted`.
`BasisError` exposes `index`, `left`, and `right`; the
attributes identify the first element or critical pair reported by basis
validation. A polynomial from another ring raises `RingError`.
