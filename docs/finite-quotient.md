# Finite quotient algebras

A reduced Gröbner basis defines a finite quotient when its leading monomial
ideal contains a pure power of every variable. `GroebnerBasis.finite_quotient`
checks this condition and stores the standard monomials.

`FiniteQuotient.standard_monomials` lists exponent vectors in ascending
grevlex order. `FiniteQuotient.dimension` counts those vectors.
The quotient also reduces residues, returns their coordinates, and builds
multiplication matrices. Operations that compute or allocate accept
`timeout` and `memory_limit` keywords.

## A verified quotient over a prime field

The certificate path is available over `F_p`. This example uses the
structured SymPy adapter to build the input. The adapter keeps the explicit
generator order and the exact `GF(p)` domain.

```python
import sympy as sp
import sylvester
from sylvester.sympy import from_sympy, to_sympy

p = 101
ring = sylvester.PolynomialRing.prime_field(p, ["x", "y"])
x_sym, y_sym = sp.symbols("x y")


def prime(expr):
    return from_sympy(sp.Poly(expr, x_sym, y_sym, modulus=p), ring)


x = prime(x_sym)
y = prime(y_sym)
ideal = ring.ideal([prime(x_sym**2), prime(y_sym**2)])

certified = ideal.groebner_basis_certified()
basis = certified.basis
verified = sylvester.verify(certified.certificate)
assert verified.modulus == p
assert to_sympy(x) == sp.Poly(x_sym, x_sym, y_sym, modulus=p)

quotient = basis.finite_quotient()
assert quotient.dimension == 4
assert {tuple(exponents) for exponents in quotient.standard_monomials} == {
    (0, 0),
    (1, 0),
    (0, 1),
    (1, 1),
}
```

The four standard monomials are `1`, `y`, `x`, and `x*y`. The quotient
relation `x^2 = y^2 = 0` reduces higher powers. `reduce` returns a polynomial
residue, and `coordinates` follows the standard monomial list.

```python
residue = quotient.reduce(prime(x_sym**3 + x_sym*y_sym + 1))
assert residue == prime(x_sym*y_sym + 1)
coordinates = quotient.coordinates(residue)
assert len(coordinates) == quotient.dimension
```

## Nilpotents separate characteristic and minimal polynomials

The multiplication operator for `x` acts on a four-dimensional space. It is
nilpotent and satisfies `x^2 = 0`, while `x` itself is nonzero. Its minimal
polynomial is therefore `t^2`. The characteristic polynomial records the
whole four-dimensional operator and is `t^4`.

```python
matrix = quotient.multiplication_matrix(x)
assert matrix.dimension == 4

characteristic = quotient.characteristic_polynomial(x)
minimal = quotient.minimal_polynomial(x)
assert str(characteristic) == "t^4"
assert str(minimal) == "t^2"
```

The matrix columns follow `quotient.standard_monomials`. `characteristic_polynomial`
and `minimal_polynomial` take the residue class represented by a polynomial.

## Rational computations carry a caveat

The rational driver is multimodular and heuristic. It has no certificate
path. A `contains_input` lift record establishes that the returned basis is
the reduced Gröbner basis of an ideal containing the input ideal. It does not
establish equality with the input ideal.

The following small system has an exact rational minimal polynomial for
multiplication by `x`.

```python
qring = sylvester.PolynomialRing.rationals(["x", "y"])
qx_sym, qy_sym = sp.symbols("x y")


def rational(expr):
    return from_sympy(sp.Poly(expr, qx_sym, qy_sym, domain=sp.QQ), qring)


qideal = qring.ideal(
    [
        rational(qx_sym**2 + qy_sym**2 - 1),
        rational(4*qx_sym*qy_sym - 1),
    ]
)
qbasis = qideal.groebner_basis(stop="contains_input")
qquotient = qbasis.finite_quotient()
qminimal = qquotient.minimal_polynomial(rational(qx_sym))

assert str(qminimal) == "t^4 - t^2 + 1/16"
assert qbasis.lift()["established"] == "contains_input"

checked = qideal.check_basis_equality(
    qbasis,
    timeout=10.0,
    memory_limit=64 * 1024 * 1024,
)
assert checked.input.generators == qideal.generators
assert [str(value) for value in checked.basis] == [str(value) for value in qbasis]
```

The explicit
`qideal.check_basis_equality(qbasis, timeout=..., memory_limit=...)` call checks
both ideal inclusions over `Q` with exact rational arithmetic. It returns a
`RationalEqualityCheck` with the candidate and its origin cofactors. This is a
library check. Independent certificate verification is available over `F_p`.

## SymPy conversion limits

`from_sympy` accepts a `sympy.Poly` over `QQ`, `ZZ`, or a prime `GF(p)`.
`ZZ` maps to exact rational coefficients. A rational coefficient can map to
an explicit prime ring when its denominator has an inverse modulo that prime.
The adapter rejects inexact, expression, composite finite-field, and
extension domains.

The SymPy generator list gives the variable order. With an explicit Sylvester
ring, names and order must match exactly. Symbols with assumptions are
rejected because Sylvester variables carry names only. The adapter converts
terms and coefficients directly. It does not evaluate expression strings or
translate Gröbner orders.
