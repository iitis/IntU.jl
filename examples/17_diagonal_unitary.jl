using IntU
using IntU
using Symbolics
using LinearAlgebra

println("=== Diagonal Unitary Integration (Torus Group) ===\n")

@variables d
@variables V[1:10, 1:10]::Complex

# Define the Diagonal Unitary Measure
# This measure integrates over matrices of form V = diag(exp(iθ_1), ..., exp(iθ_d))
println("--- 1. Define Measure ---")
measure = dDiagUnitary(V, d)
println("Measure defined: ", measure)
println("")

# Example 1: Simple Phase Average
# E[ |V_11|^2 ] = 1
println("--- 2. Single Entry Moment ---")
println("Computing E[ |V_11|^2 ]...")
expr1 = abs2(V[1, 1])
res1 = integrate(expr1, measure)
println("Result: ", res1)
println("Expected: 1")
println("")

# Example 2: Independent Phases
# E[ V_11 * V_22^* ] = 0
println("--- 3. Independent Phases ---")
println("Computing E[ V_11 * conj(V_22) ]...")
expr2 = V[1, 1] * conj(V[2, 2])
res2 = integrate(expr2, measure)
println("Result: ", res2)
println("Expected: 0")
println("")

# Example 3: Non-diagonal entries
# Diagonal unitary matrices have V_ij = 0 for i != j.
# Integration over the torus correctly handles this.
println("--- 4. Non-diagonal entries ---")
println("Computing E[ |V_12|^2 ]...")
expr3 = abs2(V[1, 2])
res3 = integrate(expr3, measure)
println("Result: ", res3)
println("Expected: 0")
println("")

# Example 4: Higher Order Correlation
# E[ |V_11|^2 * |V_22|^2 ] = 1
println("--- 5. Higher Order Correlation ---")
println("Computing E[ |V_11|^2 * |V_22|^2 ]...")
expr4 = abs2(V[1, 1]) * abs2(V[2, 2])
res4 = integrate(expr4, measure)
println("Expected: 1")
println("")

# Example 5: Matrix Integration
println("--- 6. Matrix Integration ---")
println("Integrating V * V' (should be Identity)")
# Collect to ensure standard matrix structure
V_mat = collect(V)
res_V = integrate(V_mat * V_mat', measure)

println("Result[1,1]: ", res_V[1, 1])
println("Expected: 1")
println("Result[1,2]: ", res_V[1, 2])
println("Expected: 0")
println("")

println("\nDone.")
