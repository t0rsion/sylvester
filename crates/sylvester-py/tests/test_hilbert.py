"""The Hilbert series and the dimension read off it."""

import sylvester


def test_the_series_of_a_monomial_ideal(prime_ring):
    ideal = prime_ring.ideal([prime_ring.parse("x^2"), prime_ring.parse("x*y")])
    series = ideal.groebner_basis().hilbert_series()
    assert series.numerator == [1, 0, -2, 1]
    assert series.denominator_power == 3
    assert str(series) == "(1 - 2*t^2 + t^3)/(1 - t)^3"
    assert series.dimension() == 2


def test_the_coefficients_count_the_standard_monomials(prime_ring):
    ideal = prime_ring.ideal([prime_ring.parse("x^2"), prime_ring.parse("y^2")])
    series = ideal.groebner_basis().hilbert_series()
    assert [series.coefficient(degree) for degree in range(4)] == [1, 3, 4, 4]


def test_a_zero_dimensional_quotient_reports_its_multiplicity(prime_ring):
    ideal = prime_ring.ideal(
        [
            prime_ring.parse("x^2 - 1"),
            prime_ring.parse("y^2 - 1"),
            prime_ring.parse("z^2 - 1"),
        ]
    )
    series = ideal.groebner_basis().hilbert_series()
    assert series.dimension() == 0
    assert series.multiplicity() == 8


def test_the_unit_ideal_has_no_dimension(prime_ring):
    basis = prime_ring.ideal([prime_ring.parse("1")]).groebner_basis()
    series = basis.hilbert_series()
    assert series.numerator == []
    assert series.dimension() is None
    assert series.multiplicity() is None
    assert basis.krull_dimension() is None


def test_the_zero_ideal_has_the_dimension_of_the_ring(prime_ring):
    basis = prime_ring.ideal([]).groebner_basis()
    assert basis.krull_dimension() == 3
    assert basis.hilbert_series().numerator == [1]


def test_krull_dimension_agrees_with_the_series(prime_ring, cyclic):
    basis = prime_ring.ideal(cyclic(prime_ring, 3)).groebner_basis()
    assert basis.krull_dimension() == basis.hilbert_series().dimension()
    assert basis.krull_dimension() == 0


def test_the_series_reads_leading_monomials_alone(rational_ring, prime_ring):
    over_q = rational_ring.ideal(
        [rational_ring.parse("x + y + z"), rational_ring.parse("x*y*z - 1")]
    ).groebner_basis()
    over_p = prime_ring.ideal(
        [prime_ring.parse("x + y + z"), prime_ring.parse("x*y*z - 1")]
    ).groebner_basis()
    assert over_q.hilbert_series().numerator == over_p.hilbert_series().numerator
    assert over_q.krull_dimension() == over_p.krull_dimension()
