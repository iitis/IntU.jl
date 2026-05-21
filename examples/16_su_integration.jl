using IntegrateUnitary
using Symbolics
using LinearAlgebra

println("=== SU(d) Integration (Stable Range) ===\n")

@variables d
U = SymbolicMatrix(:U, :U)

println("--- 1. Define Measure ---")
println("Using measure dSU(d)")

# Example 1: Balanced Moment
# E[ |U_11|^2 ] = 1/d
println("\n--- 2. Balanced Moment ---")
println("Computing E[ |U_11|^2 ]...")
res_balanced = integrate(abs(U[1, 1])^2, dSU(d))
println("Result: ", res_balanced)
println("Expected: 1/d")

# Example 2: Matrix Integration (SU(d) -> Identity)
println("\n--- 3. Matrix Integration ---")
println("Integrating U * U' over SU(3) (should be Identity)")
res_mat = integrate(U * U', dSU(3))
display(res_mat)

# Example 3: Higher Order Moment
# E[ |U_11 U_22|^2 ] = 1 / (d^2 - 1)
println("\n--- 4. Higher Order Moment ---")
println("Computing E[ |U_11 U_22|^2 ]...")
res_higher = integrate(abs(U[1, 1])^2 * abs(U[2, 2])^2, dSU(d))
println("Result: ", Symbolics.simplify(res_higher))
println("Expected: 1 / (d^2 - 1)")

println("\nDone.")
