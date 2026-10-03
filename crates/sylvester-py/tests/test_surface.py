"""The public operators, quotient algebra, and result records."""

import subprocess
import sys
from fractions import Fraction

import pytest
import sylvester


def test_ring_generators_and_polynomial_operators(prime_ring):
    x, y, z = prime_ring.gens
    assert prime_ring.zero.is_zero()
    assert prime_ring.one == 1 + prime_ring.zero
    assert x + y == y + x
    assert x - y == -(y - x)
    assert 3 * x == x * 3
    assert (x + y) ** 2 == x**2 + 2 * x * y + y**2
    assert (x + y).add(z, timeout=1.0) == x + y + z
    assert (x + y).sub(1) == x + y - 1
    assert (x + y).mul(2) == 2 * x + 2 * y
    assert (x + y).neg() == -(x + y)
    assert (x + y).pow(0) == prime_ring.one


def test_rational_operator_coercion_is_exact():
    ring = sylvester.PolynomialRing.rationals(["x", "y"])
    x, y = ring.gens
    expression = Fraction(1, 2) * x + y
    assert expression == ring.parse("1/2*x + y")
    assert expression.mul(Fraction(2, 3)) == ring.parse("1/3*x + 2/3*y")
    assert expression + (1, 2) == ring.parse("1/2*x + y + 1/2")


@pytest.mark.parametrize("method", ["add", "sub", "mul"])
def test_explicit_arithmetic_rejects_unsupported_operands(prime_ring, method):
    with pytest.raises(TypeError):
        getattr(prime_ring.gens[0], method)(object())


def test_operator_arithmetic_keeps_python_reflected_dispatch(prime_ring):
    class RightOperand:
        def __radd__(self, left):
            return left

    assert prime_ring.gens[0] + RightOperand() == prime_ring.gens[0]


def test_explicit_arithmetic_budget_reports_timeout(prime_ring):
    with pytest.raises(sylvester.Timeout):
        prime_ring.parse("x^2 + y", timeout=0.0)
    with pytest.raises(sylvester.Timeout):
        prime_ring.parse("x^2 + y").mul(prime_ring.parse("x + y"), timeout=0.0)


def test_division_returns_one_quotient_per_basis_element(prime_ring):
    x, y, z = prime_ring.gens
    basis = prime_ring.ideal([x**2, y**2, z**2]).groebner_basis()
    value = x**3 + x * y + z + 1
    quotients, remainder = basis.divide(value)
    assert len(quotients) == len(basis)
    rebuilt = remainder
    for quotient, divisor in zip(quotients, basis):
        rebuilt = rebuilt + quotient * divisor
    assert rebuilt == value
    assert basis.divmod(value) == (quotients, remainder)


def test_nilpotent_finite_quotient_exposes_coordinates_and_spectra():
    ring = sylvester.PolynomialRing.prime_field(5, ["x", "y"])
    x, y = ring.gens
    basis = ring.ideal([x**2, y**2]).groebner_basis()
    quotient = basis.finite_quotient()

    assert basis.is_zero_dimensional()
    assert quotient.dimension == 4
    assert quotient.dimension == len(quotient.standard_monomials)
    assert quotient.basis == quotient.standard_monomials

    residue = quotient.reduce(x**3 + x * y + 1)
    assert residue == x * y + 1
    coordinates = quotient.coordinates(residue)
    assert len(coordinates) == quotient.dimension
    assert all(isinstance(value, int) for value in coordinates)
    assert quotient.multiply(x, y) == x * y
    assert quotient.pow(x, 2).is_zero()
    assert quotient.add(x, y) == x + y

    matrix = quotient.multiplication_matrix(x)
    assert matrix.dimension == 4
    assert len(matrix.entries) == matrix.dimension
    assert all(len(row) == matrix.dimension for row in matrix.entries)
    assert all(value == 0 or isinstance(value, int) for row in matrix.entries for value in row)
    assert matrix.flat_entries == [value for row in matrix.entries for value in row]
    assert matrix.entry(-1, 0) is None
    assert matrix.entry(10**100, 0) is None
    assert matrix.entry(matrix.dimension, 0) is None
    assert str(quotient.characteristic_polynomial(x)) == "t^4"
    minimal = quotient.minimal_polynomial(x)
    assert str(minimal) == "t^2"
    assert minimal.coefficient(-1) is None
    assert minimal.coefficient(10**100) is None


def test_nonfinite_quotient_is_rejected(prime_ring):
    x, y, z = prime_ring.gens
    basis = prime_ring.ideal([x * y, z]).groebner_basis()
    with pytest.raises(ValueError):
        basis.finite_quotient()


def test_large_term_exponent_is_a_value_error(prime_ring):
    with pytest.raises(ValueError):
        prime_ring.polynomial([(1, [2**63, 0, 0])])


def test_rational_equality_check_retains_origins_and_envelope(rational_ring):
    ideal = rational_ring.ideal([rational_ring.parse("x"), rational_ring.parse("y")])
    basis = ideal.groebner_basis()
    checked = ideal.check_basis_equality(basis)

    assert checked.input.generators == ideal.generators
    assert [str(value) for value in checked.basis] == [str(value) for value in basis]
    assert len(checked.origins) == len(basis)
    assert all(len(row) == len(ideal.generators) for row in checked.origins)
    assert checked.metrics["raw_basis_elements"] >= len(basis)

    record = sylvester.ResultEnvelope.from_checked_rational(checked)
    assert record.claimed_provenance == "equals_input"
    loaded = sylvester.ResultEnvelope.from_json(record.to_json())
    assert loaded.claimed_provenance == "equals_input"
    assert loaded.domain == {"kind": "rationals"}
    recovered = loaded.check_rational()
    assert recovered.input.generators == ideal.generators
    assert [str(value) for value in recovered.basis] == [str(value) for value in basis]


def test_prime_envelope_revalidates_its_certificate_and_verifier_caps(prime_ring, cyclic):
    ideal = prime_ring.ideal(cyclic(prime_ring, 3))
    accepted = ideal.groebner_basis_certified()
    record = sylvester.ResultEnvelope.from_prime(
        ideal, accepted.basis, accepted.certificate
    )
    assert record.claimed_provenance == "unverified"

    certified_record = sylvester.ResultEnvelope.from_certified(ideal, accepted)

    verified = certified_record.verify_prime(
        max_bytes=1 << 20, max_work_units=1 << 20
    )
    assert verified.modulus == prime_ring.modulus
    assert certified_record.claimed_provenance == "certified"

    with pytest.raises(sylvester.LimitExceeded) as error:
        sylvester.verify(accepted.certificate, max_bytes=1)
    assert error.value.limit == 1


def test_result_envelope_checks_json_budget_before_decoding(prime_ring):
    with pytest.raises(sylvester.MemoryLimitExceeded):
        sylvester.ResultEnvelope.from_json(b"{}" * 128, memory_limit=128)


def test_quotient_keeps_rational_lift_provenance(rational_ring):
    ideal = rational_ring.ideal(
        [
            rational_ring.parse("x^2 - 1"),
            rational_ring.parse("y^2 - 1"),
            rational_ring.parse("z^2 - 1"),
        ]
    )
    basis = ideal.groebner_basis()
    quotient = basis.finite_quotient()
    assert quotient.source_basis.lift() == basis.lift()


def test_keyboard_interrupt_cancels_and_joins_the_worker():
    script = r'''
import _thread
import threading
import sylvester

ring = sylvester.PolynomialRing.prime_field(32003, [f"x{i}" for i in range(1, 9)])
generators = []
for degree in range(1, 8):
    terms = []
    for start in range(8):
        exponents = [0] * 8
        for step in range(degree):
            exponents[(start + step) % 8] += 1
        terms.append((1, exponents))
    generators.append(ring.polynomial(terms))
exponents = [1] * 8
generators.append(ring.polynomial([(1, exponents), (-1, [0] * 8)]))
timer = threading.Timer(0.05, _thread.interrupt_main)
try:
    timer.start()
    ring.ideal(generators).groebner_basis(threads=1)
except KeyboardInterrupt:
    print("keyboard-interrupt")
else:
    raise SystemExit("the computation finished before the interrupt")
finally:
    timer.cancel()
    timer.join()
'''
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "keyboard-interrupt" in completed.stdout
