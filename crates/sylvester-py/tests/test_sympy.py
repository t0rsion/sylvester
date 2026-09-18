"""The optional structured SymPy adapter."""

from fractions import Fraction

import pytest

sympy = pytest.importorskip("sympy")

import sylvester
import sylvester.sympy as adapter
from sylvester.sympy import from_sympy, to_sympy


def test_rational_poly_round_trips_with_exact_coefficients():
    x, y = sympy.symbols("x y")
    source = sympy.Poly(
        10**80 * x**2 - sympy.Rational(7, 13) * y + sympy.Rational(5, 11),
        x,
        y,
        domain=sympy.QQ,
    )

    converted = from_sympy(source)

    assert converted.ring == sylvester.PolynomialRing.rationals(["x", "y"])
    assert converted.terms() == [
        (Fraction(10**80, 1), [2, 0]),
        (Fraction(-7, 13), [0, 1]),
        (Fraction(5, 11), [0, 0]),
    ]
    assert to_sympy(converted) == source
    assert to_sympy(converted).domain == sympy.QQ


def test_prime_poly_round_trips_with_the_prime_field_domain():
    x, y = sympy.symbols("x y")
    source = sympy.Poly(-(x**2) - 2 * y + 1, x, y, modulus=5)

    converted = from_sympy(source)

    assert converted.ring == sylvester.PolynomialRing.prime_field(5, ["x", "y"])
    assert converted.terms() == [(4, [2, 0]), (3, [0, 1]), (1, [0, 0])]
    assert to_sympy(converted) == source
    assert to_sympy(converted).get_modulus() == 5


def test_rational_coefficients_reduce_into_an_explicit_prime_ring():
    x, y = sympy.symbols("x y")
    source = sympy.Poly(
        sympy.Rational(1, 2) * x + sympy.Rational(3, 5) * y,
        x,
        y,
        domain=sympy.QQ,
    )
    ring = sylvester.PolynomialRing.prime_field(7, ["x", "y"])

    converted = from_sympy(source, ring)

    assert converted.terms() == [(4, [1, 0]), (2, [0, 1])]
    assert to_sympy(converted) == sympy.Poly(-3 * x + 2 * y, x, y, modulus=7)


def test_a_denominator_divisible_by_the_prime_is_rejected():
    x = sympy.Symbol("x")
    source = sympy.Poly(sympy.Rational(1, 5) * x, x, domain=sympy.QQ)
    ring = sylvester.PolynomialRing.prime_field(5, ["x"])

    with pytest.raises(sylvester.RingError, match="denominator"):
        from_sympy(source, ring)


def test_an_explicit_ring_requires_the_same_generator_order():
    x, y = sympy.symbols("x y")
    source = sympy.Poly(x + 2 * y, x, y, domain=sympy.QQ)
    ring = sylvester.PolynomialRing.rationals(["y", "x"])

    with pytest.raises(ValueError, match="names and order"):
        from_sympy(source, ring)


def test_the_poly_generators_define_the_order_when_no_ring_is_given():
    x, y = sympy.symbols("x y")
    source = sympy.Poly(x + 2 * y, y, x, domain=sympy.QQ)

    converted = from_sympy(source)

    assert converted.ring.variables == ["y", "x"]
    assert converted.terms() == [(Fraction(2, 1), [1, 0]), (Fraction(1, 1), [0, 1])]


def test_zero_and_constant_keep_their_generators_and_domain():
    x, y = sympy.symbols("x y")
    for source in (
        sympy.Poly(0, x, y, domain=sympy.QQ),
        sympy.Poly(17, x, y, modulus=11),
    ):
        converted = from_sympy(source)
        result = to_sympy(converted)
        assert result.gens == (x, y)
        assert result.is_zero == source.is_zero
        assert result == source


def test_an_integer_ring_domain_maps_to_exact_rationals():
    x = sympy.Symbol("x")
    source = sympy.Poly(10**80 * x + 1, x, domain=sympy.ZZ)

    converted = from_sympy(source)

    assert converted.ring == sylvester.PolynomialRing.rationals(["x"])
    assert converted.terms() == [(Fraction(10**80, 1), [1]), (Fraction(1, 1), [0])]
    assert to_sympy(converted) == sympy.Poly(10**80 * x + 1, x, domain=sympy.QQ)


def test_a_composite_finite_field_domain_is_rejected():
    x = sympy.Symbol("x")
    source = sympy.Poly(x + 1, x, modulus=4)

    with pytest.raises(ValueError, match="finite-field modulus"):
        from_sympy(source)


def test_a_zero_variable_ring_is_rejected_by_the_sympy_adapter():
    ring = sylvester.PolynomialRing.prime_field(5, [])

    with pytest.raises(ValueError, match="zero-variable"):
        to_sympy(ring.polynomial([]))


@pytest.mark.parametrize(
    "source",
    [
        lambda x: sympy.Poly(sympy.Float("0.5") * x, x),
        lambda x: sympy.Poly(x + sympy.sqrt(2), x, extension=True),
        lambda x: sympy.Poly(x + sympy.sqrt(2), x, domain=sympy.EX),
    ],
)
def test_inexact_and_extension_domains_are_rejected(source):
    x = sympy.Symbol("x")
    with pytest.raises(ValueError, match="unsupported SymPy coefficient domain"):
        from_sympy(source(x))


def test_symbols_with_assumptions_are_rejected():
    x = sympy.Symbol("x", positive=True)
    source = sympy.Poly(x + 1, x, domain=sympy.QQ)

    with pytest.raises(ValueError, match="cannot carry assumptions"):
        from_sympy(source)


def test_missing_sympy_dependency_reports_the_install_name(monkeypatch):
    def missing(_name):
        raise ModuleNotFoundError("No module named 'sympy'", name="sympy")

    monkeypatch.setattr(adapter, "import_module", missing)

    with pytest.raises(ImportError, match="pip install sympy"):
        adapter._load_sympy()


def test_incompatible_prime_and_rational_domains_are_rejected():
    x = sympy.Symbol("x")
    source = sympy.Poly(x + 1, x, modulus=5)
    ring = sylvester.PolynomialRing.rationals(["x"])

    with pytest.raises(ValueError, match="incompatible"):
        from_sympy(source, ring)


def test_incompatible_prime_moduli_are_rejected():
    x = sympy.Symbol("x")
    source = sympy.Poly(x + 1, x, modulus=5)
    ring = sylvester.PolynomialRing.prime_field(7, ["x"])

    with pytest.raises(ValueError, match="moduli"):
        from_sympy(source, ring)


def test_the_adapter_rejects_expressions_and_non_polynomials():
    x = sympy.Symbol("x")
    with pytest.raises(TypeError, match="sympy.Poly"):
        from_sympy(x + 1)
    with pytest.raises(TypeError, match="sylvester.Polynomial"):
        to_sympy(x + 1)
