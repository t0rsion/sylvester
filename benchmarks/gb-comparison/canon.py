"""Canonical form for a reduced Gröbner basis, shared by driver.py and
report.py.

Equal leading ideals do not imply equal bases. The cross-check compares
the whole basis, coefficients included, in this one module. A
tool-specific parse must not invent a second notion of equality.

A basis canonicalizes to a sorted tuple of monic polynomials over F_p.
Each polynomial is a tuple of (exponent tuple, coefficient) pairs sorted
grevlex descending. Two bases compare equal (`==`) exactly when this
form matches, term for term, coefficient for coefficient, mod p.
"""

import re
from functools import cmp_to_key

MONO_RE = re.compile(r"x(\d+)(?:\^(\d+))?")


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
