"""Shared input systems for the sylvester test suite."""

import pytest

import sylvester


def cyclic_generators(ring, n):
    """The cyclic-n system in the first n variables of the ring.

    The first n - 1 generators are the cyclic sums of products of d
    consecutive variables, and the last one is the product of all of them
    minus 1.
    """
    generators = []
    for degree in range(1, n):
        terms = []
        for start in range(n):
            exponents = [0] * n
            for step in range(degree):
                exponents[(start + step) % n] += 1
            terms.append((1, exponents))
        generators.append(ring.polynomial(terms))
    generators.append(ring.polynomial([(1, [1] * n), (-1, [0] * n)]))
    return generators


@pytest.fixture
def cyclic():
    """A builder for the cyclic-n system over a given ring."""
    return cyclic_generators


@pytest.fixture
def prime_ring():
    """F_32003[x, y, z], the ring most prime-field cases use."""
    return sylvester.PolynomialRing.prime_field(32003, ["x", "y", "z"])


@pytest.fixture
def rational_ring():
    """Q[x, y, z], the ring most rational cases use."""
    return sylvester.PolynomialRing.rationals(["x", "y", "z"])
