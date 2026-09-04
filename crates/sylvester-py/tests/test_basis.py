"""The ideal and the basis: computation, the sequence protocol, division."""

import pytest

import sylvester


def test_cyclic_four_over_a_prime_field_has_seven_elements(cyclic):
    ring = sylvester.PolynomialRing.prime_field(32003, ["x", "y", "z", "w"])
    basis = ring.ideal(cyclic(ring, 4)).groebner_basis()
    assert len(basis) == 7
    assert str(basis[-1]) == "x + y + z + w"
    assert str(basis[0]) == "z^2*w^4 + y*z + 32002*y*w + z*w + 32001*w^2"


def test_the_two_backends_compute_one_basis(prime_ring, cyclic):
    ideal = prime_ring.ideal(cyclic(prime_ring, 3))
    f4 = ideal.groebner_basis(backend="f4")
    classic = ideal.groebner_basis(backend="classic")
    assert [str(f) for f in f4] == [str(f) for f in classic]


def test_the_basis_is_a_sequence(prime_ring, cyclic):
    basis = prime_ring.ideal(cyclic(prime_ring, 3)).groebner_basis()
    assert len(basis) == 3
    assert basis[0] == basis[-3]
    assert basis[-1] == list(basis)[2]
    assert basis[0:2] == [basis[0], basis[1]]
    assert basis[::-1] == [basis[2], basis[1], basis[0]]
    assert [str(f) for f in basis] == [str(f) for f in list(basis)]
    with pytest.raises(IndexError):
        basis[3]


def test_the_in_operator_reads_the_list_and_not_the_ideal(prime_ring):
    """`__contains__` is not implemented, so `in` scans the listed
    polynomials. The zero polynomial belongs to every ideal and is in no
    list, which separates the two questions."""
    basis = prime_ring.ideal(
        [prime_ring.parse("x^2 - 1"), prime_ring.parse("x*y - 1")]
    ).groebner_basis()
    zero = prime_ring.polynomial([])
    assert basis[0] in basis
    assert zero not in basis
    assert basis.contains(zero)


def test_normal_form_is_zero_exactly_on_a_member(prime_ring):
    f = prime_ring.parse("x^2 - 1")
    g = prime_ring.parse("x*y - 1")
    basis = prime_ring.ideal([f, g]).groebner_basis()
    assert basis.normal_form(f).is_zero()
    assert basis.contains(f)
    assert basis.contains(g)

    outside = prime_ring.parse("z")
    assert not basis.contains(outside)
    remainder = basis.normal_form(outside)
    assert str(remainder) == "z"
    assert basis.normal_form(remainder) == remainder


def test_every_generator_belongs_to_the_ideal_it_generates(prime_ring, cyclic):
    generators = cyclic(prime_ring, 3)
    basis = prime_ring.ideal(generators).groebner_basis()
    for generator in generators:
        assert basis.contains(generator)
    assert basis.contains(prime_ring.polynomial([]))


def test_normal_form_rejects_a_polynomial_of_another_ring(prime_ring, rational_ring):
    basis = prime_ring.ideal([prime_ring.parse("x")]).groebner_basis()
    with pytest.raises(sylvester.RingError):
        basis.normal_form(rational_ring.parse("x"))
    other = sylvester.PolynomialRing.prime_field(7, ["x", "y", "z"])
    with pytest.raises(sylvester.RingError):
        basis.contains(other.parse("x"))


def test_from_polynomials_accepts_a_computed_basis(prime_ring, cyclic):
    computed = prime_ring.ideal(cyclic(prime_ring, 3)).groebner_basis()
    checked = sylvester.GroebnerBasis.from_polynomials(prime_ring, list(computed))
    assert [str(f) for f in checked] == [str(f) for f in computed]


@pytest.mark.parametrize(
    "texts, index",
    [
        (["2*x^2", "y"], 0),
        (["x", "y^2"], 1),
    ],
)
def test_from_polynomials_names_the_element_that_failed(prime_ring, texts, index):
    polynomials = [prime_ring.parse(text) for text in texts]
    with pytest.raises(sylvester.BasisError) as caught:
        sylvester.GroebnerBasis.from_polynomials(prime_ring, polynomials)
    assert caught.value.index == index


def test_from_polynomials_names_the_pair_that_is_not_closed(prime_ring):
    polynomials = [prime_ring.parse("x^2 - y"), prime_ring.parse("x*y - z")]
    with pytest.raises(sylvester.BasisError) as caught:
        sylvester.GroebnerBasis.from_polynomials(prime_ring, polynomials)
    assert (caught.value.left, caught.value.right) == (0, 1)


def test_the_report_carries_the_backend_and_the_f4_counters(prime_ring, cyclic):
    ideal = prime_ring.ideal(cyclic(prime_ring, 3))
    basis, report = ideal.groebner_basis_with_report()
    assert [str(f) for f in basis] == [str(f) for f in ideal.groebner_basis()]
    assert report.backend == "f4"
    assert report.elapsed > 0.0
    assert report.threads_used >= 1
    assert report.modular is None
    assert report.modular_concurrency is None
    assert report.counters["batches"] > 0
    assert report.counters["matrix_rows"] > 0


def test_the_classic_backend_counts_nothing(prime_ring, cyclic):
    _, report = prime_ring.ideal(cyclic(prime_ring, 3)).groebner_basis_with_report(
        backend="classic"
    )
    assert report.backend == "classic"
    assert report.counters is None
    assert report.threads_used == 1


def test_homogeneity_is_read_off_the_generators_and_decided_on_the_basis(prime_ring):
    inhomogeneous = prime_ring.ideal(
        [prime_ring.parse("x^2 - z^2"), prime_ring.parse("x*y - z^2")]
    )
    assert inhomogeneous.has_homogeneous_generators()
    assert inhomogeneous.groebner_basis().is_homogeneous()

    mixed = prime_ring.ideal([prime_ring.parse("x^2 + y")])
    assert not mixed.has_homogeneous_generators()
    assert not mixed.groebner_basis().is_homogeneous()
