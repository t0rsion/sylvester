"""The multimodular path over Q: the basis, the lift record, and the stop."""

from fractions import Fraction

import pytest

import sylvester


def cyclic_three(ring):
    return ring.ideal(
        [
            ring.parse("x + y + z"),
            ring.parse("x*y + y*z + z*x"),
            ring.parse("x*y*z - 1"),
        ]
    )


def test_cyclic_three_over_the_rationals(rational_ring):
    basis = cyclic_three(rational_ring).groebner_basis()
    assert [str(f) for f in basis] == [
        "z^3 - 1",
        "y^2 + y*z + z^2",
        "x + y + z",
    ]


def test_a_rational_basis_holds_exact_fractions():
    ring = sylvester.PolynomialRing.rationals(["x", "y"])
    basis = ring.ideal(
        [ring.parse("2*x - 3*y"), ring.parse("y^2 - 5")]
    ).groebner_basis()
    assert [str(f) for f in basis] == ["y^2 - 5", "x - 3/2*y"]
    assert basis[1].terms() == [(Fraction(1, 1), [1, 0]), (Fraction(-3, 2), [0, 1])]


def test_the_lift_counters_add_up(rational_ring):
    basis = cyclic_three(rational_ring).groebner_basis()
    lift = basis.lift()
    assert lift["primes_folded"] + lift["primes_discarded"] == lift["primes_consumed"]
    assert lift["primes_consumed"] >= 3
    assert lift["confirming_primes"] == 2
    assert lift["modulus_bits"] > 0
    assert lift["established"] == "contains_input"


def test_extra_primes_sets_the_number_of_confirmations(rational_ring):
    basis = cyclic_three(rational_ring).groebner_basis(extra_primes=4)
    assert basis.lift()["confirming_primes"] == 4


def test_the_exact_tests_report_what_they_establish(rational_ring):
    basis = cyclic_three(rational_ring).groebner_basis(stop="unchanged", extra_primes=1)
    assert basis.lift()["established"] == "unchanged"


def test_the_report_carries_the_lift_and_no_f4_counters(rational_ring):
    basis, report = cyclic_three(rational_ring).groebner_basis_with_report()
    assert report.counters is None
    assert report.threads_used == 1
    assert report.modular == basis.lift()
    assert report.modular_concurrency >= 1


def test_a_thread_count_changes_no_byte_of_the_basis(rational_ring):
    ideal = cyclic_three(rational_ring)
    one = ideal.groebner_basis(threads=1)
    eight = ideal.groebner_basis(threads=8)
    assert [str(f) for f in one] == [str(f) for f in eight]
    assert one.lift() == eight.lift()


def test_a_checked_rational_basis_has_no_lift(rational_ring):
    computed = cyclic_three(rational_ring).groebner_basis()
    checked = sylvester.GroebnerBasis.from_polynomials(rational_ring, list(computed))
    assert checked.lift() is None


def test_a_prime_field_basis_has_no_lift_at_all(prime_ring):
    basis = prime_ring.ideal([prime_ring.parse("x")]).groebner_basis()
    with pytest.raises(ValueError):
        basis.lift()


def test_the_rationals_have_no_certified_path(rational_ring):
    with pytest.raises(ValueError):
        cyclic_three(rational_ring).groebner_basis_certified()


def test_normal_form_over_the_rationals(rational_ring):
    basis = cyclic_three(rational_ring).groebner_basis()
    assert basis.contains(rational_ring.parse("x*y*z - 1"))
    assert str(basis.normal_form(rational_ring.parse("z^3"))) == "1"
