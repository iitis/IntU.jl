using IntU
using Symbolics
using LinearAlgebra

println("=== SU(d) Integration (Stable Range) ===\n")

@variables d
@symbolic_dimension U[1:d, 1:d]

# Define the Special Unitary Measure
# In the stable range (d large), SU(d) integrals coincide with U(d) for balanced moments.
println("--- 1. Define Measure ---")
measure = dSU(U, d)
println("Measure defined: ", measure)
println("")

# Example 1: Balanced Moment
# E[ |U_11|^2 ] = 1/d
println("--- 2. Balanced Moment ---")
println("Computing E[ |U_11|^2 ]...")
expr_balanced = U[1, 1] * conj(U[1, 1])
res_balanced = integrate(expr_balanced, measure)
println("Result: ", res_balanced)
println("Expected: 1/d")
println("")

# Example 2: Unbalanced Moment
# E[ U_11 ] = 0
println("--- 3. Unbalanced Moment ---")
println("Computing E[ U_11 ]...")
# Note: For SU(d), moments like E[ det(U) ] would be 1, but typical single entries are 0.
# The current implementation assumes stable range where unbalanced terms vanish.
expr_unbalanced = U[1, 1]
res_unbalanced = integrate(expr_unbalanced, measure)
println("Result: ", res_unbalanced)
println("Expected: 0")
println("")

# Example 3: Higher Order Moment
# E[ |U_11 U_22|^2 ] = 1 / (d^2 - 1)
println("--- 4. Higher Order Moment ---")
println("Computing E[ |U_11 U_22|^2 ]...")
expr_higher = abs2(U[1, 1]) * abs2(U[2, 2])
res_higher = Symbolics.simplify(integrate(expr_higher, measure); expand=true)
println("Result: ", res_higher)
println("Expected: 1 / (d^2 - 1)")

# Verify numerically (symbolic check)
expected = 1 / (d^2 - 1)
diff = Symbolics.simplify(res_higher - expected; expand=true)
println("Difference from expected: ", diff)

println("\nDone.")
