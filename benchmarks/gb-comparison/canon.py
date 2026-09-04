"""Canonical form for a reduced Gröbner basis, shared by driver.py and
report.py.

Equal leading ideals do not imply equal bases. The cross-check compares
the whole basis, coefficients included, in this one module. A
tool-specific parse must not invent a second notion of equality.

A basis canonicalizes to a sorted tuple of monic polynomials over F_p.
Each polynomial is a tuple of (exponent tuple, coefficient) pairs sorted
grevlex descending. Two bases compare equal (`==`) exactly when this
form matches, term for term, coefficient for coefficient, mod p.

The rational functions below (the `_q` suffix) do the same job over `Q`,
with exact `fractions.Fraction` coefficients in place of a residue mod
`p`. The three reference tools disagree on normalization going in: msolve
content-clears (divides the integer coefficients, after clearing
denominators, by their gcd) rather than dividing by the leading
coefficient, and Singular needs `simplify(std(I), 1)` rather than
`option(redSB)` alone to reach a monic basis over `Q` (over `F_p`,
`option(redSB)` is already monic, because every nonzero leading
coefficient there has an inverse and Singular applies it). Groebner.jl's
`groebner` is monic already. `monic_q` below is the one place every
tool's basis is forced to the same normalization, the rational
counterpart of `monic`'s division by the leading coefficient mod `p`.
"""

import re
from fractions import Fraction
from functools import cmp_to_key

MONO_RE = re.compile(r"x(\d+)(?:\^(\d+))?")
Q_COEFF_RE = re.compile(r"^(\d+)(?:(//|/)(\d+))?")


def grevlex_cmp(a, b):
    """-1, 0, or 1: grevlex comparison of two exponent tuples, x1 > x2 > ... > xn."""
    da, db = sum(a), sum(b)
    if da != db:
        return -1 if da < db else 1
    for x, y in zip(reversed(a), reversed(b)):
        if x != y:
            return 1 if x < y else -1
    return 0


GREVLEX_KEY = cmp_to_key(grevlex_cmp)


def parse_mono(s, nvars):
    """'x2^2*x3' or 'x1' or '1' -> an exponent tuple of length nvars."""
    exps = [0] * nvars
    s = s.strip()
    if s in ("1", ""):
        return tuple(exps)
    for m in MONO_RE.finditer(s):
        idx = int(m.group(1)) - 1
        e = int(m.group(2)) if m.group(2) else 1
        exps[idx] += e
    return tuple(exps)


def split_terms(expr):
    """Split a flat sum of monomials at its top-level '+' and '-'.

    Every engine prints a polynomial as a flat sum: no parentheses, no
    unary minus inside an exponent. A plus/minus split is safe without a
    full expression parser.
    """
    expr = expr.replace(" ", "")
    return [t for t in re.findall(r"[+-]?[^+-]+", expr) if t]


def parse_term(chunk, nvars):
    """One signed term, '-3*x1^2*x2' or '7' or 'x1^2', to (exponents, coeff).

    Handles both conventions: an implicit coefficient of 1 and an implicit
    exponent of 1 (Singular, M2, Groebner.jl), and an explicit coefficient
    and exponent on every factor (msolve's own output, "1*x1^1").
    """
    sign = 1
    if chunk[0] == "-":
        sign, chunk = -1, chunk[1:]
    elif chunk[0] == "+":
        chunk = chunk[1:]
    m = re.match(r"(\d+)", chunk)
    if m:
        coeff = int(m.group(1))
        rest = chunk[m.end() :]
        rest = rest.removeprefix("*")
    else:
        coeff = 1
        rest = chunk
    return parse_mono(rest, nvars), sign * coeff


def parse_full_poly(expr, nvars, p):
    """A tool's printed polynomial expression to {exponent tuple: coeff mod p}."""
    expr = expr.strip()
    poly = {}
    if expr in ("", "0"):
        return poly
    for chunk in split_terms(expr):
        exps, coeff = parse_term(chunk, nvars)
        c = (poly.get(exps, 0) + coeff) % p
        if c == 0:
            poly.pop(exps, None)
        else:
            poly[exps] = c
    return poly


def poly_from_terms(term_list, p):
    """[(coeff, exponent tuple), ...] to {exponent tuple: coeff mod p}.

    For the runners that already return parsed (coefficient, exponent
    vector) pairs (sylvester, msolve-inproc) instead of a text expression.
    """
    poly = {}
    for coeff, exps in term_list:
        exps = tuple(exps)
        c = (poly.get(exps, 0) + coeff) % p
        if c == 0:
            poly.pop(exps, None)
        else:
            poly[exps] = c
    return poly


def monic(poly, p):
    """A poly dict to a monic term tuple, sorted grevlex descending.

    Raises ValueError on the zero polynomial: a reduced Gröbner basis
    never contains one, and a canonicalizer that silently dropped it
    would hide the defect that put it there.
    """
    if not poly:
        raise ValueError("zero polynomial in a basis")
    lt = max(poly, key=GREVLEX_KEY)
    inv = pow(poly[lt], p - 2, p)
    ordered = sorted(poly.items(), key=lambda kv: GREVLEX_KEY(kv[0]), reverse=True)
    return tuple((e, (c * inv) % p) for e, c in ordered)


def canon_basis(polys, p):
    """A list of poly dicts to the canonical, comparable form of the basis.

    Two engines whose bases agree exactly, coefficients included, produce
    equal values. The elements are sorted independently of the order an
    engine returned them in.
    """
    return tuple(sorted(monic(poly, p) for poly in polys))


def lms_of(polys):
    """A list of poly dicts to the sorted list of leading-monomial tuples.

    The weaker check: size and leading-monomial agreement without full
    coefficient agreement.
    """
    return sorted(max(poly, key=GREVLEX_KEY) for poly in polys)


def parse_term_q(chunk, nvars):
    """One signed rational term to (exponents, Fraction).

    A coefficient is a decimal integer, optionally followed by a
    denominator: Singular and M2 write "1/6" and msolve writes no
    fraction at all (it clears denominators, so every coefficient is an
    integer); Groebner.jl writes "1//6", Julia's own `Rational` format.
    `Q_COEFF_RE` reads both separators. A coefficient of 1 or -1 is
    implicit, as in the mod-`p` case: "x1^2" alone, not "1*x1^2".
    """
    sign = 1
    if chunk[0] == "-":
        sign, chunk = -1, chunk[1:]
    elif chunk[0] == "+":
        chunk = chunk[1:]
    m = Q_COEFF_RE.match(chunk)
    if m and m.group(1):
        numerator = int(m.group(1))
        denominator = int(m.group(3)) if m.group(3) else 1
        rest = chunk[m.end() :].removeprefix("*")
    else:
        numerator, denominator, rest = 1, 1, chunk
    return parse_mono(rest, nvars), Fraction(sign * numerator, denominator)


def parse_full_poly_q(expr, nvars):
    """A tool's printed rational polynomial expression to
    {exponent tuple: Fraction}."""
    expr = expr.strip()
    poly = {}
    if expr in ("", "0"):
        return poly
    for chunk in split_terms(expr):
        exps, coeff = parse_term_q(chunk, nvars)
        c = poly.get(exps, Fraction(0)) + coeff
        if c == 0:
            poly.pop(exps, None)
        else:
            poly[exps] = c
    return poly


def poly_from_terms_q(term_list):
    """[(Fraction, exponent tuple), ...] to {exponent tuple: Fraction}.

    For sylv-runner's `rational` mode, which hands back parsed
    (coefficient, exponent vector) pairs directly, the rational
    counterpart of `poly_from_terms`.
    """
    poly = {}
    for coeff, exps in term_list:
        exps = tuple(exps)
        c = poly.get(exps, Fraction(0)) + coeff
        if c == 0:
            poly.pop(exps, None)
        else:
            poly[exps] = c
    return poly


def monic_q(poly):
    """A poly dict of Fractions to a monic term tuple, sorted grevlex
    descending.

    Raises ValueError on the zero polynomial, as `monic` does.
    """
    if not poly:
        raise ValueError("zero polynomial in a basis")
    lt = max(poly, key=GREVLEX_KEY)
    inv = Fraction(1) / poly[lt]
    ordered = sorted(poly.items(), key=lambda kv: GREVLEX_KEY(kv[0]), reverse=True)
    return tuple((e, c * inv) for e, c in ordered)


def canon_basis_q(polys):
    """A list of poly dicts of Fractions to the canonical, comparable form
    of a rational basis.

    Four tools reach this: msolve content-cleared, Singular and
    Groebner.jl already close to monic (see the module docstring for
    each one's own normalization), and sylvester's basis, which the
    multimodular engine already returns monic. `monic_q` makes every one
    of them monic the same way, so agreement here means the whole basis
    agrees, coefficients included, not only its leading terms.
    """
    return tuple(sorted(monic_q(poly) for poly in polys))
