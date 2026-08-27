using Groebner, AbstractAlgebra
# warm-up on a tiny system (cyclic-3) so compilation is excluded from timing
let
    Rw, (a, b, c) = polynomial_ring(GF(1073741827), ["a", "b", "c"], internal_ordering=:degrevlex)
    groebner([a + b + c, a*b + b*c + c*a, a*b*c - 1], ordering=DegRevLex())
end
R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = polynomial_ring(GF(1073741827), ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"], internal_ordering=:degrevlex)
sys = [
  x1*x2*x10+x2*x3*x10+x3*x4*x10+x4*x5*x10+x5*x6*x10+x6*x7*x10+x7*x8*x10+x8*x9*x10+x1*x10-1,
  x1*x3*x10+x2*x4*x10+x3*x5*x10+x4*x6*x10+x5*x7*x10+x6*x8*x10+x7*x9*x10+x2*x10-2,
  x1*x4*x10+x2*x5*x10+x3*x6*x10+x4*x7*x10+x5*x8*x10+x6*x9*x10+x3*x10-3,
  x1*x5*x10+x2*x6*x10+x3*x7*x10+x4*x8*x10+x5*x9*x10+x4*x10-4,
  x1*x6*x10+x2*x7*x10+x3*x8*x10+x4*x9*x10+x5*x10-5,
  x1*x7*x10+x2*x8*x10+x3*x9*x10+x6*x10-6,
  x1*x8*x10+x2*x9*x10+x7*x10-7,
  x1*x9*x10+x8*x10-8,
  x9*x10-9,
  x1+x2+x3+x4+x5+x6+x7+x8+x9+1
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
