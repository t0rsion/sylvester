using Groebner, AbstractAlgebra
# warm-up on a tiny system (cyclic-3) so compilation is excluded from timing
let
    Rw, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"], internal_ordering=:degrevlex)
    groebner([a + b + c, a*b + b*c + c*a, a*b*c - 1], ordering=DegRevLex())
end
R, (x1, x2, x3, x4, x5, x6, x7, x8, x9) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9"], internal_ordering=:degrevlex)
sys = [
  x1*x2*x9+x2*x3*x9+x3*x4*x9+x4*x5*x9+x5*x6*x9+x6*x7*x9+x7*x8*x9+x1*x9-1,
  x1*x3*x9+x2*x4*x9+x3*x5*x9+x4*x6*x9+x5*x7*x9+x6*x8*x9+x2*x9-2,
  x1*x4*x9+x2*x5*x9+x3*x6*x9+x4*x7*x9+x5*x8*x9+x3*x9-3,
  x1*x5*x9+x2*x6*x9+x3*x7*x9+x4*x8*x9+x4*x9-4,
  x1*x6*x9+x2*x7*x9+x3*x8*x9+x5*x9-5,
  x1*x7*x9+x2*x8*x9+x6*x9-6,
  x1*x8*x9+x7*x9-7,
  x8*x9-8,
  x1+x2+x3+x4+x5+x6+x7+x8+1
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
