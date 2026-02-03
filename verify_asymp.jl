using IntU, Symbolics
@variables d
@symbolic_dimension U[1:d, 1:d]
measure = dU(U)
expr = abs(U[1,1])^4
res = integrate(expr, measure)
println("Exact result: ", res)
asymp = asymptotic(expr, measure, 4)
println("Asymptotic: ", asymp)
