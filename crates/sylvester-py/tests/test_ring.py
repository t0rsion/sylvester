"""The ring and the polynomial: construction, text, terms, and equality."""

from fractions import Fraction

import pytest

import sylvester


def test_prime_field_ring_reports_its_modulus_and_variables():
    ring = sylvester.PolynomialRing.prime_field(32003, ["x", "y", "z"])
    assert ring.modulus == 32003
    assert ring.variables == ["x", "y", "z"]


def test_rational_ring_has_no_modulus():
    ring = sylvester.PolynomialRing.rationals(["x", "y"])
    assert ring.modulus is None
    assert ring.variables == ["x", "y"]


def test_two_constructions_of_one_ring_are_equal_and_hash_equally():
    first = sylvester.PolynomialRing.prime_field(7, ["x", "y"])
    second = sylvester.PolynomialRing.prime_field(7, ["x", "y"])
    assert first == second
    assert hash(first) == hash(second)
    assert first != sylvester.PolynomialRing.prime_field(11, ["x", "y"])
    assert first != sylvester.PolynomialRing.rationals(["x", "y"])


@pytest.mark.parametrize(
    "modulus, variables",
    [
        (4, ["x"]),
        (1 << 32, ["x"]),
        (7, ["x", "x"]),
        (7, ["2x"]),
    ],
)
def test_a_rejected_ring_raises_ring_error(modulus, variables):
    with pytest.raises(sylvester.RingError):
        sylvester.PolynomialRing.prime_field(modulus, variables)


def test_ring_error_is_a_value_error():
    with pytest.raises(ValueError):
        sylvester.PolynomialRing.rationals(["x", "x"])


def test_parse_and_str_round_trip_over_a_prime_field(prime_ring):
    for text in ["x^2*y + 32000*z + 1", "x*z + 32002*y", "1"]:
        assert str(prime_ring.parse(text)) == text


def test_parse_and_str_round_trip_over_the_rationals(rational_ring):
    for text in ["1/2*x - y", "2/3*x^2 - y*z + 1", "-1/2", "x + y + z"]:
        assert str(rational_ring.parse(text)) == text


def test_a_prime_field_coefficient_reduces_and_a_rational_one_stays_exact():
    ring = sylvester.PolynomialRing.prime_field(7, ["x"])
    assert str(ring.parse("-1")) == "6"
    assert ring.parse("8*x") == ring.parse("x")

    rational = sylvester.PolynomialRing.rationals(["x"])
    assert rational.parse("2/4*x").terms() == [(Fraction(1, 2), [1])]


def test_terms_run_from_the_largest_monomial_down(prime_ring):
    f = prime_ring.parse("x^2*y - 3*z + 1")
    assert f.terms() == [(1, [2, 1, 0]), (32000, [0, 0, 1]), (1, [0, 0, 0])]
    assert f.degree() == 3
    assert not f.is_zero()


def test_the_zero_polynomial_has_no_degree(prime_ring):
    zero = prime_ring.polynomial([])
    assert zero.is_zero()
    assert zero.degree() is None
    assert zero.terms() == []


def test_polynomial_takes_ints_fractions_and_pairs(rational_ring):
    f = rational_ring.polynomial(
        [(Fraction(1, 2), [1, 0, 0]), ((3, 4), [0, 1, 0]), (5, [0, 0, 1])]
    )
    assert str(f) == "1/2*x + 3/4*y + 5*z"


def test_repeated_monomials_add_up_and_a_zero_term_drops_out(prime_ring):
    f = prime_ring.polynomial([(1, [1, 0, 0]), (2, [1, 0, 0]), (0, [0, 1, 0])])
    assert str(f) == "3*x"


def test_a_wrong_exponent_count_is_a_ring_error(prime_ring):
    with pytest.raises(sylvester.RingError):
        prime_ring.polynomial([(1, [1, 0])])


def test_an_exponent_outside_the_width_of_one_exponent_is_rejected(prime_ring):
    with pytest.raises(ValueError):
        prime_ring.polynomial([(1, [65536, 0, 0])])
    with pytest.raises(ValueError):
        prime_ring.polynomial([(1, [-1, 0, 0])])


def test_a_denominator_the_prime_divides_is_a_ring_error():
    ring = sylvester.PolynomialRing.prime_field(3, ["x"])
    with pytest.raises(sylvester.RingError):
        ring.parse("1/3*x")
    assert str(ring.parse("3/3*x")) == "x"


def test_polynomials_compare_and_hash_by_content(prime_ring):
    f = prime_ring.parse("x^2*y - 3*z + 1")
    g = prime_ring.polynomial([(1, [2, 1, 0]), (-3, [0, 0, 1]), (1, [0, 0, 0])])
    assert f == g
    assert hash(f) == hash(g)
    assert f != prime_ring.parse("x")
    assert f != "x^2*y - 3*z + 1"
    assert {f, g} == {f}


def test_a_polynomial_carries_its_ring(prime_ring):
    assert prime_ring.parse("x").ring == prime_ring


def test_an_ideal_takes_polynomials_of_its_own_ring_alone(prime_ring, rational_ring):
    with pytest.raises(sylvester.RingError):
        prime_ring.ideal([rational_ring.parse("x")])
