"""Structured conversion between Sylvester and SymPy polynomials.

SymPy is an optional dependency. Importing :mod:`sylvester` does not load it.
The adapter accepts ``sympy.Poly`` values over ``QQ``, ``ZZ``, and ``GF(p)``.
It does not translate expression strings, coefficient extensions, or
monomial orders.
"""

from __future__ import annotations

from fractions import Fraction
from importlib import import_module
from typing import Any

from . import Polynomial, PolynomialRing

__all__ = ["from_sympy", "to_sympy"]


def _load_sympy() -> Any:
    """Load SymPy and report the optional dependency requirement."""

    try:
        return import_module("sympy")
    except ImportError as error:
        raise ImportError(
            "SymPy support requires the optional dependency. "
            "Install it with `pip install sympy`."
        ) from error


def _generator_names(sympy: Any, poly: Any) -> list[str]:
    """Read plain, assumption-free SymPy symbols in generator order."""

    names: list[str] = []
    for generator in poly.gens:
        if type(generator) is not sympy.Symbol:
            raise TypeError("SymPy polynomial generators must be plain Symbols")
        if generator.assumptions0 != {"commutative": True}:
            raise ValueError("SymPy polynomial generators cannot carry assumptions")
        names.append(generator.name)

    if len(names) != len(set(names)):
        raise ValueError("SymPy polynomial generators must have distinct names")
    return names


def _source_domain(sympy: Any, poly: Any) -> tuple[str, int | None]:
    """Classify a supported SymPy coefficient domain."""

    domain = poly.domain
    if domain.is_FiniteField:
        modulus = domain.mod
        if not isinstance(modulus, int) or not sympy.isprime(modulus):
            raise ValueError(
                f"unsupported SymPy finite-field modulus {modulus}; use a prime p"
            )
        return "prime", modulus
    if domain.is_RationalField or domain.is_IntegerRing:
        return "rational", None
    raise ValueError(
        f"unsupported SymPy coefficient domain {domain}; use QQ, ZZ, or GF(p)"
    )


def _coefficient(sympy: Any, value: Any) -> int | tuple[int, int]:
    """Read an exact integer or rational coefficient from SymPy."""

    if isinstance(value, sympy.Rational):
        return (int(value.p), int(value.q))
    raise TypeError(f"unsupported SymPy coefficient {value!r}")


def _target_ring(
    source_kind: str,
    source_modulus: int | None,
    names: list[str],
    ring: PolynomialRing | None,
) -> PolynomialRing:
    """Infer a ring or validate the supplied ring against the source."""

    if ring is None:
        if source_kind == "prime":
            return PolynomialRing.prime_field(source_modulus, names)
        return PolynomialRing.rationals(names)
    if not isinstance(ring, PolynomialRing):
        raise TypeError("ring must be a sylvester.PolynomialRing")
    if ring.variables != names:
        raise ValueError(
            "SymPy generator names and order must match the Sylvester ring"
        )

    if source_kind == "prime":
        if ring.modulus is None:
            raise ValueError(
                "a SymPy GF(p) polynomial is incompatible with a rational "
                "Sylvester ring"
            )
        if ring.modulus != source_modulus:
            raise ValueError("SymPy and Sylvester prime moduli do not match")
    return ring


def from_sympy(
    poly: Any,
    ring: PolynomialRing | None = None,
) -> Polynomial:
    """Convert a SymPy ``Poly`` to a Sylvester ``Polynomial``.

    The SymPy generators define the variable names and order when ``ring``
    is omitted. With ``ring``, names and order must match exactly. ``QQ`` and
    ``ZZ`` map to a rational ring. ``GF(p)`` maps to a prime-field ring.
    Rational coefficients reduce modulo a supplied prime when their
    denominators are invertible.
    """

    sympy = _load_sympy()
    if not isinstance(poly, sympy.Poly):
        raise TypeError("from_sympy expects a sympy.Poly value")

    names = _generator_names(sympy, poly)
    source_kind, source_modulus = _source_domain(sympy, poly)
    ring = _target_ring(source_kind, source_modulus, names, ring)

    terms = [
        (_coefficient(sympy, coefficient), list(exponents))
        for exponents, coefficient in poly.terms()
    ]
    return ring.polynomial(terms)


def _output_coefficient(sympy: Any, value: Any) -> Any:
    """Convert a Sylvester coefficient to an exact SymPy value."""

    if isinstance(value, Fraction):
        return sympy.Rational(value.numerator, value.denominator)
    if isinstance(value, int):
        return sympy.Integer(value)
    raise TypeError(f"unsupported Sylvester coefficient {value!r}")


def to_sympy(poly: Polynomial) -> Any:
    """Convert a Sylvester ``Polynomial`` to a SymPy ``Poly``.

    The result uses fresh assumption-free symbols in the Sylvester ring's
    variable order. Rational values use ``QQ`` and prime-field values use
    ``GF(p)``. Zero-variable rings are unsupported because SymPy cannot
    construct a ``Poly`` without generators.
    """

    sympy = _load_sympy()
    if not isinstance(poly, Polynomial):
        raise TypeError("to_sympy expects a sylvester.Polynomial value")

    ring = poly.ring
    if not ring.variables:
        raise ValueError("to_sympy does not support zero-variable rings")
    symbols = tuple(sympy.Symbol(name) for name in ring.variables)
    data = {
        tuple(exponents): _output_coefficient(sympy, coefficient)
        for coefficient, exponents in poly.terms()
    }

    if ring.modulus is None:
        return sympy.Poly.from_dict(data, *symbols, domain=sympy.QQ)
    return sympy.Poly.from_dict(data, *symbols, modulus=ring.modulus)
