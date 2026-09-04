"""Benchmark input generator.

Each system is built once as a list of polynomials. Each polynomial is a
dict from exponent tuple to integer coefficient. The same data is rendered
into every tool's input format.

Conventions:
- Variables are always x1 > x2 > ... > xN (grevlex everywhere).
- cyclic-n: n vars. f_d = sum_{i=0}^{n-1} prod_{j=0}^{d-1} x_{(i+j) mod n}
  for d = 1..n-1, and f_n = x1*...*xn - 1.
- katsura-n (classical POSSO/msolve definition): n+1 vars u_0..u_n mapped to
  x1..x_{n+1} (u_m = x_{m+1}).
  Equations: u_0 + 2*sum_{i=1}^{n} u_i - 1 = 0 and, for m = 0..n-1,
  sum_{i=-n}^{n} u_i * u_{m-i} - u_m = 0, where u_{-i} = u_i and u_i = 0
  for |i| > n.
- eco-n: n vars. f_k = (x_k + sum_{i=1}^{n-k-1} x_i * x_{i+k}) * x_n - k
  for k = 1..n-1, and f_n = x_1 + ... + x_{n-1} + 1.
- noon-n (Noonburg): n vars.
  f_i = 10 * x_i * sum_{j != i} x_j^2 - 11 * x_i + 10, for i = 1..n.

Field: F_p with p = 1073741827 for every tool.
"""

import os
from collections import defaultdict

P = 1073741827
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.path.join(HERE, "inputs")


def add_term(poly, exps, coeff):
    poly[tuple(exps)] += coeff


def cyclic(n):
    polys = []
    for d in range(1, n):
        poly = defaultdict(int)
        for i in range(n):
            exps = [0] * n
            for j in range(d):
                exps[(i + j) % n] += 1
            add_term(poly, exps, 1)
        polys.append(dict(poly))
    poly = defaultdict(int)
    add_term(poly, [1] * n, 1)
    add_term(poly, [0] * n, -1)
    polys.append(dict(poly))
    return n, polys


def katsura(n):
    nv = n + 1  # u_0..u_n -> x1..x_{n+1}
    polys = []
    # u_0 + 2*sum u_i - 1
    poly = defaultdict(int)
    for i in range(n + 1):
        exps = [0] * nv
        exps[i] = 1
        add_term(poly, exps, 1 if i == 0 else 2)
    add_term(poly, [0] * nv, -1)
    polys.append(dict(poly))
    # for m = 0..n-1: sum_{i=-n}^{n} u_i u_{m-i} - u_m
    for m in range(n):
        poly = defaultdict(int)
        for i in range(-n, n + 1):
            j = m - i
            if abs(j) > n:
                continue
            exps = [0] * nv
            exps[abs(i)] += 1
            exps[abs(j)] += 1
            add_term(poly, exps, 1)
        exps = [0] * nv
        exps[m] = 1
        add_term(poly, exps, -1)
        polys.append(dict(poly))
    return nv, polys


def eco(n):
    polys = []
    for k in range(1, n):
        poly = defaultdict(int)
        exps = [0] * n
        exps[k - 1] += 1
        exps[n - 1] += 1
        add_term(poly, exps, 1)
        # sum_{i=1}^{n-k-1} x_i x_{i+k} x_n
        for i in range(1, n - k):
            exps = [0] * n
            exps[i - 1] += 1
            exps[i + k - 1] += 1
            exps[n - 1] += 1
            add_term(poly, exps, 1)
        add_term(poly, [0] * n, -k)
        polys.append(dict(poly))
    poly = defaultdict(int)
    for i in range(1, n):
        exps = [0] * n
        exps[i - 1] = 1
        add_term(poly, exps, 1)
    add_term(poly, [0] * n, 1)
    polys.append(dict(poly))
    return n, polys


def noon(n):
    polys = []
    for i in range(n):
        poly = defaultdict(int)
        for j in range(n):
            if j == i:
                continue
            exps = [0] * n
            exps[i] += 1
            exps[j] += 2
            add_term(poly, exps, 10)
        exps = [0] * n
        exps[i] = 1
        add_term(poly, exps, -11)
        add_term(poly, [0] * n, 10)
        polys.append(dict(poly))
    return n, polys


def poly_str(poly, mul="*"):
    """Render canonical poly to a text form all four external tools accept."""
    parts = [signed_term(exps, coeff, mul) for exps, coeff in ordered_terms(poly)]
    if not parts:
        return "0"
    first_sign, first_term = parts[0]
    first = (first_sign if first_sign == "-" else "") + first_term
    return first + "".join(sign + term for sign, term in parts[1:])


def ordered_terms(poly):
    """Return nonzero terms in descending grevlex order."""
    key = lambda item: (-sum(item[0]), tuple(reversed(item[0])))
    return [(exps, coeff) for exps, coeff in sorted(poly.items(), key=key) if coeff]


def signed_term(exps, coeff, mul):
    """Render one nonzero term as a sign and unsigned body."""
    factors = [
        f"x{index + 1}^{exponent}" if exponent > 1 else f"x{index + 1}"
        for index, exponent in enumerate(exps)
        if exponent
    ]
    monomial = mul.join(factors)
    magnitude = abs(coeff)
    if not monomial:
        body = str(magnitude)
    elif magnitude == 1:
        body = monomial
    else:
        body = f"{magnitude}{mul}{monomial}"
    return ("-" if coeff < 0 else "+", body)


def render_sylvester(name, nv, polys):
    lines = [str(nv)]
    for poly in polys:
        terms = []
        for exps, coeff in sorted(poly.items()):
            if coeff == 0:
                continue
            terms.append(",".join([str(coeff)] + [str(e) for e in exps]))
        lines.append(";".join(terms))
    with open(os.path.join(INP, f"{name}.syl"), "w") as f:
        f.write("\n".join(lines) + "\n")


def render_singular(name, nv, polys):
    vars_ = ",".join(f"x{i + 1}" for i in range(nv))
    body = ",\n  ".join(poly_str(p) for p in polys)
    s = f"""system("--ticks-per-sec",1000);
ring r = {P},({vars_}),dp;
ideal I = {body};
option(redSB);
int t0 = rtimer;
ideal G = std(I);
int t1 = rtimer;
print("TIME_MS "+string(t1-t0));
print("SIZE "+string(size(G)));
int i;
for (i=1; i<=size(G); i++) {{ print("LM "+string(leadmonom(G[i]))); }}
for (i=1; i<=size(G); i++) {{ print("POLY "+string(G[i])); }}
exit;
"""
    with open(os.path.join(INP, f"{name}.sing"), "w") as f:
        f.write(s)


def render_msolve(name, nv, polys):
    vars_ = ",".join(f"x{i + 1}" for i in range(nv))
    body = ",\n".join(poly_str(p) for p in polys)
    with open(os.path.join(INP, f"{name}.ms"), "w") as f:
        f.write(f"{vars_}\n{P}\n{body}\n")


def render_m2(name, nv, polys):
    vars_ = ",".join(f"x{i + 1}" for i in range(nv))
    body = ",\n  ".join(poly_str(p) for p in polys)
    s = f"""R = ZZ/{P}[{vars_}, MonomialOrder => GRevLex];
I = ideal(
  {body}
);
t = elapsedTiming groebnerBasis I;
G = t#1;
<< "TIME_S " << t#0 << endl;
<< "SIZE " << numColumns G << endl;
scan(first entries leadTerm G, m -> << "LM " << toString m << endl);
scan(first entries G, m -> << "POLY " << toString m << endl);
exit 0;
"""
    with open(os.path.join(INP, f"{name}.m2"), "w") as f:
        f.write(s)


def render_julia(name, nv, polys):
    vars_list = ", ".join(f"x{i + 1}" for i in range(nv))
    vars_strs = ", ".join(f'"x{i + 1}"' for i in range(nv))
    body = ",\n  ".join(poly_str(p) for p in polys)
    s = f"""using Groebner, AbstractAlgebra
# warm-up on a tiny system (cyclic-3) so compilation is excluded from timing
let
    Rw, (a, b, c) = polynomial_ring(GF({P}), ["a", "b", "c"], internal_ordering=:degrevlex)
    groebner([a + b + c, a*b + b*c + c*a, a*b*c - 1], ordering=DegRevLex())
end
R, ({vars_list}) = polynomial_ring(GF({P}), [{vars_strs}], internal_ordering=:degrevlex)
sys = [
  {body}
];
G = nothing
for r in 1:3
    local t = @elapsed Gr = groebner(sys, ordering=DegRevLex())
    global G = Gr
    println("RUN ", t)
    flush(stdout)
    if t > 120.0
        println("TIMEOUT")
        exit(0)
    end
end
println("SIZE ", length(G))
for f in G
    println("LM ", leading_monomial(f))
end
for f in G
    println("POLY ", f)
end
"""
    with open(os.path.join(INP, f"{name}.jl"), "w") as f:
        f.write(s)


def render_sylvester_q(name, nv, polys):
    """The `.sylq` format: the same shape as `.syl` (see render_sylvester),
    under its own extension rather than a characteristic line, because the
    domain here comes from sylv-runner's own `rational` mode argument, not
    from the file. See runner/src/main.rs for the reader."""
    lines = [str(nv)]
    for poly in polys:
        terms = []
        for exps, coeff in sorted(poly.items()):
            if coeff == 0:
                continue
            terms.append(",".join([str(coeff)] + [str(e) for e in exps]))
        lines.append(";".join(terms))
    with open(os.path.join(INP, f"{name}.sylq"), "w") as f:
        f.write("\n".join(lines) + "\n")


def render_singular_q(name, nv, polys):
    """Like render_singular, over the rationals. `option(redSB)` alone
    content-clears rather than normalizing to monic over Q (unlike over
    F_p, where the leading coefficient is already 1 after reduction), so
    this adds `simplify(std(I), 1)`, confirmed against the installed
    Singular 4.4.1 (`system(--ticks-per-sec) ... option(redSB); ideal G =
    simplify(std(I), 1);` prints a monic reduced basis; plain
    `option(redSB); std(I)` does not)."""
    vars_ = ",".join(f"x{i + 1}" for i in range(nv))
    body = ",\n  ".join(poly_str(p) for p in polys)
    s = f"""system("--ticks-per-sec",1000);
ring r = 0,({vars_}),dp;
ideal I = {body};
option(redSB);
int t0 = rtimer;
ideal G = simplify(std(I), 1);
int t1 = rtimer;
print("TIME_MS "+string(t1-t0));
print("SIZE "+string(size(G)));
int i;
for (i=1; i<=size(G); i++) {{ print("LM "+string(leadmonom(G[i]))); }}
for (i=1; i<=size(G); i++) {{ print("POLY "+string(G[i])); }}
exit;
"""
    with open(os.path.join(INP, f"{name}.sing"), "w") as f:
        f.write(s)


def render_msolve_q(name, nv, polys):
    """Like render_msolve, with the characteristic line `0`. The installed
    msolve 0.10.1 then emits a rational basis with `-g 2`: content-cleared
    (the gcd of the cleared-denominator integer coefficients divided out),
    not normalized to monic. canon.py divides by the leading coefficient to
    compare it against the other three tools."""
    vars_ = ",".join(f"x{i + 1}" for i in range(nv))
    body = ",\n".join(poly_str(p) for p in polys)
    with open(os.path.join(INP, f"{name}.ms"), "w") as f:
        f.write(f"{vars_}\n0\n{body}\n")


def render_julia_q(name, nv, polys):
    """Like render_julia, over `QQ` in place of `GF(P)`. The warm-up system
    is over `QQ` too, so the rational specialization of Groebner.jl and
    AbstractAlgebra.jl compiles before the first timed run, the same
    reason the prime-field warm-up runs over `GF(P)`."""
    vars_list = ", ".join(f"x{i + 1}" for i in range(nv))
    vars_strs = ", ".join(f'"x{i + 1}"' for i in range(nv))
    body = ",\n  ".join(poly_str(p) for p in polys)
    s = f"""using Groebner, AbstractAlgebra
# warm-up on a tiny system (cyclic-3) so compilation is excluded from timing
let
    Rw, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"], internal_ordering=:degrevlex)
    groebner([a + b + c, a*b + b*c + c*a, a*b*c - 1], ordering=DegRevLex())
end
R, ({vars_list}) = polynomial_ring(QQ, [{vars_strs}], internal_ordering=:degrevlex)
sys = [
  {body}
];
G = nothing
for r in 1:3
    local t = @elapsed Gr = groebner(sys, ordering=DegRevLex())
    global G = Gr
    println("RUN ", t)
    flush(stdout)
    if t > 120.0
        println("TIMEOUT")
        exit(0)
    end
end
println("SIZE ", length(G))
for f in G
    println("LM ", leading_monomial(f))
end
for f in G
    println("POLY ", f)
end
"""
    with open(os.path.join(INP, f"{name}.jl"), "w") as f:
        f.write(s)


def main():
    os.makedirs(INP, exist_ok=True)
    instances = []
    for n in range(4, 9):
        instances.append((f"cyclic-{n}", *cyclic(n)))
    for n in range(4, 11):
        instances.append((f"katsura-{n}", *katsura(n)))
    for n in range(8, 12):
        instances.append((f"eco-{n}", *eco(n)))
    for n in range(3, 7):
        instances.append((f"noon-{n}", *noon(n)))
    names = []
    for name, nv, polys in instances:
        render_sylvester(name, nv, polys)
        render_singular(name, nv, polys)
        render_msolve(name, nv, polys)
        render_m2(name, nv, polys)
        render_julia(name, nv, polys)
        names.append(name)
        print(f"{name}: {nv} vars, {len(polys)} polys")
    with open(os.path.join(INP, "instances.txt"), "w") as f:
        f.write("\n".join(names) + "\n")

    # The rational cells: the same families, at sizes small enough that
    # Singular, msolve, and Groebner.jl finish in seconds over Q, where
    # coefficient growth, not monomial count, drives cost. Every name
    # carries a "-q" suffix so its four files (.ms, .sing, .jl, .sylq)
    # never collide with the prime-field files of the same family and
    # size, which already use the unsuffixed name under four of those five
    # extensions.
    q_instances = []
    for n in (4, 5, 6):
        q_instances.append((f"cyclic-{n}-q", *cyclic(n)))
    for n in (4, 5, 6, 7):
        q_instances.append((f"katsura-{n}-q", *katsura(n)))
    for n in (8, 9):
        q_instances.append((f"eco-{n}-q", *eco(n)))
    for n in (3, 4, 5):
        q_instances.append((f"noon-{n}-q", *noon(n)))
    q_names = []
    for name, nv, polys in q_instances:
        render_sylvester_q(name, nv, polys)
        render_singular_q(name, nv, polys)
        render_msolve_q(name, nv, polys)
        render_julia_q(name, nv, polys)
        q_names.append(name)
        print(f"{name}: {nv} vars, {len(polys)} polys (Q)")
    with open(os.path.join(INP, "instances-q.txt"), "w") as f:
        f.write("\n".join(q_names) + "\n")


if __name__ == "__main__":
    main()
