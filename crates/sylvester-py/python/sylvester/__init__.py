"""Gröbner bases over prime fields and over the rational numbers.

Every value comes from a ring, which fixes the coefficient domain, the
variable names, and the variable order:

    import sylvester

    ring = sylvester.PolynomialRing.prime_field(32003, ["x", "y", "z"])
    ideal = ring.ideal([ring.parse("x^2*y - 3*z + 1"), ring.parse("x*z - y")])
    basis = ideal.groebner_basis()

The monomial order is grevlex over the variable order, and certificates
name it `grevlex-v1`. There is no other order in this release.

`Ideal.groebner_basis` returns the value the engine computed.
`Ideal.groebner_basis_certified` returns one an independent verifier
accepted, and `sylvester.verify` checks certificate bytes on their own.
Over the rational numbers the engine is multimodular and there is no
certified path; `GroebnerBasis.lift` says what the run observed.

Every computation releases the GIL, so other Python threads keep running
while it runs.
"""

from sylvester._sylvester import (
    BasisError,
    BudgetExhausted,
    CertificateInvalid,
    CertifiedGroebnerBasis,
    ComputeReport,
    FiniteQuotient,
    GroebnerBasis,
    HilbertSeries,
    Ideal,
    InternalDefect,
    LimitExceeded,
    MemoryLimitExceeded,
    MultiplicationMatrix,
    ParseError,
    Polynomial,
    PolynomialRing,
    RationalEqualityCheck,
    ResultEnvelope,
    RingError,
    SylvesterError,
    Timeout,
    UnivariatePolynomial,
    VerifiedGroebnerBasis,
    __version__,
    verify,
)

__all__ = [
    "BasisError",
    "BudgetExhausted",
    "CertificateInvalid",
    "CertifiedGroebnerBasis",
    "ComputeReport",
    "FiniteQuotient",
    "GroebnerBasis",
    "HilbertSeries",
    "Ideal",
    "InternalDefect",
    "LimitExceeded",
    "MemoryLimitExceeded",
    "MultiplicationMatrix",
    "ParseError",
    "Polynomial",
    "PolynomialRing",
    "RationalEqualityCheck",
    "ResultEnvelope",
    "RingError",
    "SylvesterError",
    "Timeout",
    "UnivariatePolynomial",
    "VerifiedGroebnerBasis",
    "__version__",
    "verify",
]
