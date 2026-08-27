using Groebner, AbstractAlgebra
# warm-up on a tiny system (cyclic-3) so compilation is excluded from timing
let
    Rw, (a, b, c) = polynomial_ring(GF(1073741827), ["a", "b", "c"], internal_ordering=:degrevlex)
    groebner([a + b + c, a*b + b*c + c*a, a*b*c - 1], ordering=DegRevLex())
end
R, (x1, x2, x3, x4) = polynomial_ring(GF(1073741827), ["x1", "x2", "x3", "x4"], internal_ordering=:degrevlex)
sys = [
  x1+x2+x3+x4,
  x1*x2+x2*x3+x1*x4+x3*x4,
  x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4,
  x1*x2*x3*x4-1
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
