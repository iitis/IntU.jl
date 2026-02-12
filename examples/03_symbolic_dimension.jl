# examples/03_symbolic_dimension.jl
using IntU
using Symbolics
using LinearAlgebra

# Define symbolic dimension
@variables d
U = SymbolicMatrix(:U, :U)

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
minor = U[1, 1]*U[2, 2] - U[1, 2]*U[2, 1]
expr3 = abs(minor)^2
res3 = @integrate expr3 dU(d)
# Theoretical expectation: 2 / (d*(d-1))
println("Result: $res3")

# Simplify results for better readability
println("\nSimplified results:")
println("1. $(Symbolics.simplify(res1))")
println("2. $(Symbolics.simplify(res2))")
println("3. $(Symbolics.simplify(res3))")

# You can also use the @symbolic_dimension macro for backward compatibility,
# but the clean interface with SymbolicMatrix is now preferred.
println("\nUsing @symbolic_dimension (Legacy):")
@symbolic_dimension V[1:d, 1:d]
res_legacy = integrate(abs(V[1,1])^2, dU(V))
println("Legacy Result: $res_legacy")
