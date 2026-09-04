"""The certified path and the independent verifier."""

import pytest

import sylvester


def certified(ring, cyclic):
    return ring.ideal(cyclic(ring, 3)).groebner_basis_certified()


def test_the_certified_basis_is_the_computed_one(prime_ring, cyclic):
    ideal = prime_ring.ideal(cyclic(prime_ring, 3))
    accepted = ideal.groebner_basis_certified()
    assert [str(f) for f in accepted.basis] == [
        str(f) for f in ideal.groebner_basis()
    ]


def test_f4_writes_the_binary_schema_and_classic_the_json_one(prime_ring, cyclic):
    ideal = prime_ring.ideal(cyclic(prime_ring, 3))
    assert ideal.groebner_basis_certified().certificate[:4] == b"SYLV"
    assert ideal.groebner_basis_certified(backend="classic").certificate[:1] == b"{"


def test_verify_reads_the_bytes_back(prime_ring, cyclic):
    accepted = certified(prime_ring, cyclic)
    verified = sylvester.verify(accepted.certificate)
    assert verified.modulus == 32003
    assert verified.nvars == 3
    assert len(verified.input) == 3
    assert len(verified.basis) == len(accepted.basis)


def test_the_verified_basis_prints_under_synthetic_names(prime_ring, cyclic):
    verified = sylvester.verify(certified(prime_ring, cyclic).certificate)
    assert verified.basis[0].ring.variables == ["x1", "x2", "x3"]
    assert "x1" in str(verified.basis[-1])


def test_a_flipped_byte_is_rejected(prime_ring, cyclic):
    data = bytearray(certified(prime_ring, cyclic).certificate)
    data[len(data) // 2] ^= 0x01
    with pytest.raises(sylvester.CertificateInvalid):
        sylvester.verify(bytes(data))


def test_bytes_that_are_no_certificate_are_rejected():
    for data in [b"", b"not a certificate"]:
        with pytest.raises(sylvester.CertificateInvalid):
            sylvester.verify(data)


def test_a_rejection_is_bad_input_and_not_a_defect(prime_ring, cyclic):
    with pytest.raises(ValueError):
        sylvester.verify(b"{}")
