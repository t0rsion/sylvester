from collections.abc import Iterator, Sequence
from fractions import Fraction
from typing import Literal, overload

__version__: str

Coefficient = int | Fraction
TermInput = tuple[int | Fraction | tuple[int, int], Sequence[int]]
Backend = Literal["f4", "classic"]
Stop = Literal["unchanged", "contains_input"]

class SylvesterError(Exception): ...
class RingError(SylvesterError, ValueError): ...
class ParseError(SylvesterError, ValueError): ...

class BasisError(SylvesterError, ValueError):
    """A list the checked constructor rejected.

    Every instance carries all three attributes. `index` names the one
    element that failed a shape check, and `left` and `right` name the
    index pair whose S-polynomial did not reduce to zero. The attributes
    the failure does not name are None.
    """

    index: int | None
    left: int | None
    right: int | None

class CertificateInvalid(SylvesterError, ValueError): ...
class BudgetExhausted(SylvesterError, RuntimeError): ...
class Timeout(BudgetExhausted): ...
class MemoryLimitExceeded(BudgetExhausted): ...

class LimitExceeded(SylvesterError, RuntimeError):
    limit: int | None

class InternalDefect(SylvesterError, RuntimeError): ...

class PolynomialRing:
    @staticmethod
    def prime_field(p: int, variables: Sequence[str]) -> PolynomialRing: ...
    @staticmethod
    def rationals(variables: Sequence[str]) -> PolynomialRing: ...
    @property
    def modulus(self) -> int | None: ...
    @property
    def variables(self) -> list[str]: ...
    def parse(self, text: str) -> Polynomial: ...
    def polynomial(self, terms: Sequence[TermInput]) -> Polynomial: ...
    def ideal(self, generators: Sequence[Polynomial]) -> Ideal: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...

class Polynomial:
    def terms(self) -> list[tuple[Coefficient, list[int]]]: ...
    def degree(self) -> int | None: ...
    def is_zero(self) -> bool: ...
    @property
    def ring(self) -> PolynomialRing: ...
    def __str__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...

class Ideal:
    @property
    def generators(self) -> list[Polynomial]: ...
    @property
    def ring(self) -> PolynomialRing: ...
    def has_homogeneous_generators(self) -> bool: ...
    def groebner_basis(
        self,
        *,
        backend: Backend | None = None,
        timeout: float | None = None,
        memory_limit: int | None = None,
        threads: int | None = None,
        stop: Stop | None = None,
        extra_primes: int | None = None,
    ) -> GroebnerBasis: ...
    def groebner_basis_with_report(
        self,
        *,
        backend: Backend | None = None,
        timeout: float | None = None,
        memory_limit: int | None = None,
        threads: int | None = None,
        stop: Stop | None = None,
        extra_primes: int | None = None,
    ) -> tuple[GroebnerBasis, ComputeReport]: ...
    def groebner_basis_certified(
        self,
        *,
        backend: Backend | None = None,
        timeout: float | None = None,
        memory_limit: int | None = None,
        threads: int | None = None,
    ) -> CertifiedGroebnerBasis: ...

class GroebnerBasis:
    """A reduced Gröbner basis under grevlex.

    The value is a sequence: `len`, indexing with negative indices and
    slices, and iteration all work on it. There is no `__contains__`, so
    `f in basis` iterates and compares: it asks whether `f` is one of the
    listed polynomials. `basis.contains(f)` is the ideal membership test.
    """

    @staticmethod
    def from_polynomials(
        ring: PolynomialRing,
        polynomials: Sequence[Polynomial],
        *,
        timeout: float | None = None,
        memory_limit: int | None = None,
    ) -> GroebnerBasis: ...
    @property
    def ring(self) -> PolynomialRing: ...
    def normal_form(
        self,
        f: Polynomial,
        *,
        timeout: float | None = None,
        memory_limit: int | None = None,
    ) -> Polynomial: ...
    def contains(
        self,
        f: Polynomial,
        *,
        timeout: float | None = None,
        memory_limit: int | None = None,
    ) -> bool:
        """Report whether f belongs to the ideal this basis generates.

        The answer is about the ideal this basis generates, which is not
        always the ideal the caller started from. Over Q,
        lift()["established"] says what relates the two. Under
        "unchanged" neither answer is established for the input ideal.
        Under "contains_input" the ideal of this basis contains the input
        ideal, so False holds for the input ideal as well, and True does
        not.
        """
    def hilbert_series(
        self, *, timeout: float | None = None, memory_limit: int | None = None
    ) -> HilbertSeries:
        """The Hilbert series of the quotient by the leading monomial ideal
        of this basis.

        Over Q the series is of the ideal this basis generates. Under
        "unchanged" nothing relates that ideal to the input; under
        "contains_input" it contains the input ideal.
        """
    def krull_dimension(
        self, *, timeout: float | None = None, memory_limit: int | None = None
    ) -> int | None:
        """The Krull dimension of the quotient by the ideal this basis
        generates, or None for the unit ideal.

        Over Q, under "contains_input" the value is at most the Krull
        dimension of the quotient by the input ideal; under "unchanged"
        nothing follows.
        """
    def is_homogeneous(self) -> bool: ...
    def lift(self) -> dict[str, int | str] | None: ...
    def __len__(self) -> int: ...
    @overload
    def __getitem__(self, index: int) -> Polynomial: ...
    @overload
    def __getitem__(self, index: slice) -> list[Polynomial]: ...
    def __iter__(self) -> Iterator[Polynomial]: ...

class CertifiedGroebnerBasis:
    @property
    def basis(self) -> GroebnerBasis: ...
    @property
    def certificate(self) -> bytes: ...

class HilbertSeries:
    @property
    def numerator(self) -> list[int]: ...
    @property
    def denominator_power(self) -> int: ...
    def dimension(self) -> int | None: ...
    def multiplicity(self) -> int | None: ...
    def coefficient(self, degree: int) -> int: ...
    def __str__(self) -> str: ...

class ComputeReport:
    @property
    def backend(self) -> Backend: ...
    @property
    def counters(self) -> dict[str, int] | None: ...
    @property
    def elapsed(self) -> float: ...
    @property
    def threads_used(self) -> int: ...
    @property
    def modular(self) -> dict[str, int | str] | None: ...
    @property
    def modular_concurrency(self) -> int | None: ...

class VerifiedGroebnerBasis:
    @property
    def modulus(self) -> int: ...
    @property
    def nvars(self) -> int: ...
    @property
    def input(self) -> list[Polynomial]: ...
    @property
    def basis(self) -> list[Polynomial]: ...

def verify(data: bytes) -> VerifiedGroebnerBasis: ...
