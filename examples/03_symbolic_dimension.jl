using IntU
using Symbolics
using LinearAlgebra

println("Testing symbolic dimension d with unified clean interface...")

# 1. Integrate |U_11|^2
println("\n1. Integrating |U_11|^2 ...")
res1 = @integrate abs(U[1, 1])^2 dU(d)
println("Result: $res1 (Expected: 1/d)")

# 2. Integrate |U_11|^4
println("\n2. Integrating |U_11|^4 ...")
res2 = @integrate abs(U[1, 1])^4 dU(d)
# Expected: 2 / (d*(d+1))
println("Result: $res2")

# 3. Integrate a 2x2 minor
println("\n3. Integrating |U_11*U_22 - U_12*U_21|^2 ...")
res3 = @integrate abs(U[1, 1]*U[2, 2] - U[1, 2]*U[2, 1])^2 dU(d)
# Theoretical expectation: 2 / (d*(d-1))
println("Result: $res3")

# Simplify results for better readability
println("\nSimplified results:")
println("1. $(Symbolics.simplify(res1))")
println("2. $(Symbolics.simplify(res2))")
println("3. $(Symbolics.simplify(res3))")

# manual declaration using factory functions
println("\nDefining another symbolic matrix V using factory function:")
V = symbolic_unitary(:V, d)
res_v = integrate(abs(V[1, 1])^2, dU(d))
println("Result for V: $res_v")
