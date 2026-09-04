"""The error mapping of `docs/rational-design.md` section 9.3, row by row."""

import pytest

import sylvester


def hard_ideal(cyclic):
    """Cyclic-7 over F_32003, which no budget below stops in time."""
    ring = sylvester.PolynomialRing.prime_field(
        32003, [f"x{index}" for index in range(1, 8)]
    )
    return ring.ideal(cyclic(ring, 7))


def test_a_parse_failure_is_a_parse_error(prime_ring):
    with pytest.raises(sylvester.ParseError) as caught:
        prime_ring.parse("x^^2")
    assert isinstance(caught.value, ValueError)
    assert isinstance(caught.value, sylvester.SylvesterError)


def test_an_unknown_variable_is_a_parse_error(prime_ring):
    with pytest.raises(sylvester.ParseError):
        prime_ring.parse("q + 1")


def test_an_exhausted_deadline_is_a_timeout(cyclic):
    with pytest.raises(sylvester.Timeout) as caught:
        hard_ideal(cyclic).groebner_basis(timeout=0.001)
    assert isinstance(caught.value, sylvester.BudgetExhausted)
    assert isinstance(caught.value, RuntimeError)
    assert isinstance(caught.value, sylvester.SylvesterError)


def test_an_exhausted_memory_limit_is_its_own_class(cyclic):
    with pytest.raises(sylvester.MemoryLimitExceeded) as caught:
        hard_ideal(cyclic).groebner_basis(memory_limit=1000)
    assert isinstance(caught.value, sylvester.BudgetExhausted)


def test_the_certified_path_reports_an_exhausted_budget_the_same_way(cyclic):
    """The budget covers the engine, the certificate writer, and the
    verifier, and each of the three reports through the same two classes."""
    with pytest.raises(sylvester.Timeout):
        hard_ideal(cyclic).groebner_basis_certified(timeout=0.001)
    with pytest.raises(sylvester.MemoryLimitExceeded):
        hard_ideal(cyclic).groebner_basis_certified(memory_limit=1000)


def test_a_structural_limit_carries_its_bound(prime_ring):
    """Dividing x^65535*y by x^65535 + y^65535 needs the tail multiple
    y^65536, which is past the width of one exponent."""
    ring = sylvester.PolynomialRing.prime_field(32003, ["x", "y"])
    divisor = ring.polynomial([(1, [65535, 0]), (1, [0, 65535])])
    basis = sylvester.GroebnerBasis.from_polynomials(ring, [divisor])
    with pytest.raises(sylvester.LimitExceeded) as caught:
        basis.normal_form(ring.polynomial([(1, [65535, 1])]))
    assert caught.value.limit == 65535
    assert isinstance(caught.value, RuntimeError)


def test_a_supplied_list_that_is_no_basis_is_a_basis_error(prime_ring):
    with pytest.raises(sylvester.BasisError) as caught:
        sylvester.GroebnerBasis.from_polynomials(
            prime_ring, [prime_ring.parse("2*x")]
        )
    assert isinstance(caught.value, ValueError)
    assert caught.value.index == 0


def test_a_basis_error_carries_all_three_attributes(prime_ring):
    """The stub declares `index`, `left`, and `right` on every instance,
    so reading one never raises AttributeError. The failure the error
    names decides which of them is not None."""
    with pytest.raises(sylvester.BasisError) as shape:
        sylvester.GroebnerBasis.from_polynomials(prime_ring, [prime_ring.parse("2*x")])
    assert shape.value.index == 0
    assert shape.value.left is None
    assert shape.value.right is None

    # The two are monic, sorted, and interreduced, and the S-polynomial
    # of the pair leaves the remainder x - y^2.
    with pytest.raises(sylvester.BasisError) as pair:
        sylvester.GroebnerBasis.from_polynomials(
            prime_ring,
            [prime_ring.parse("x^2 - y"), prime_ring.parse("x*y - 1")],
        )
    assert pair.value.index is None
    assert pair.value.left == 0
    assert pair.value.right == 1


def test_a_rational_keyword_on_a_prime_field_ideal_is_rejected(prime_ring, cyclic):
    ideal = prime_ring.ideal(cyclic(prime_ring, 3))
    for keywords in [{"stop": "unchanged"}, {"extra_primes": 3}]:
        with pytest.raises(ValueError):
            ideal.groebner_basis(**keywords)


@pytest.mark.parametrize(
    "keywords",
    [
        {"backend": "buchberger"},
        {"timeout": -1.0},
        {"timeout": float("nan")},
        {"timeout": float("inf")},
        # A finite value past what a Duration holds. The conversion
        # panics on it unless the binding checks first.
        {"timeout": 1e300},
    ],
)
def test_a_wrong_keyword_value_is_rejected(prime_ring, keywords):
    ideal = prime_ring.ideal([prime_ring.parse("x")])
    with pytest.raises(ValueError):
        ideal.groebner_basis(**keywords)


@pytest.mark.parametrize(
    "keywords",
    [{"stop": "certified"}, {"extra_primes": 0}],
)
def test_a_wrong_rational_keyword_value_is_rejected(rational_ring, keywords):
    ideal = rational_ring.ideal([rational_ring.parse("x")])
    with pytest.raises(ValueError):
        ideal.groebner_basis(**keywords)


def test_an_unknown_keyword_is_a_type_error(prime_ring):
    ideal = prime_ring.ideal([prime_ring.parse("x")])
    with pytest.raises(TypeError):
        ideal.groebner_basis(rewrite=True)


def test_every_class_of_the_package_shares_one_base():
    classes = [
        sylvester.RingError,
        sylvester.ParseError,
        sylvester.BasisError,
        sylvester.CertificateInvalid,
        sylvester.BudgetExhausted,
        sylvester.Timeout,
        sylvester.MemoryLimitExceeded,
        sylvester.LimitExceeded,
        sylvester.InternalDefect,
    ]
    for cls in classes:
        assert issubclass(cls, sylvester.SylvesterError)
    for cls in [
        sylvester.RingError,
        sylvester.ParseError,
        sylvester.BasisError,
        sylvester.CertificateInvalid,
    ]:
        assert issubclass(cls, ValueError)
    for cls in [
        sylvester.BudgetExhausted,
        sylvester.LimitExceeded,
        sylvester.InternalDefect,
    ]:
        assert issubclass(cls, RuntimeError)
