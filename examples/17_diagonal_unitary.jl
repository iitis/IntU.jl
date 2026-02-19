using IntU
using Symbolics
using LinearAlgebra

println("=== Diagonal Unitary Integration (Torus Group) ===\n")

@variables d
V = SymbolicMatrix(:V, :DiagUnitary)

# Define the Diagonal Unitary Measure
# This measure integrates over matrices of form V = diag(exp(iθ_1), ..., exp(iθ_d))
println("--- 1. Define Measure ---")
# The new interface supports dDiagUnitary(d) directly
println("Using measure dDiagUnitary(d)")

# Example 1: Simple Phase Average
# E[ |V_11|^2 ] = 1
println("\n--- 2. Single Entry Moment ---")
println("Computing E[ |V_11|^2 ]...")
res1 = integrate(abs(V[1, 1])^2, dDiagUnitary(d))
println("Result: ", res1)
println("Expected: 1")

# Example 2: Independent Phases
# E[ V_11 * V_22^* ] = 0
println("\n--- 3. Independent Phases ---")
println("Computing E[ V_11 * conj(V_22) ]...")
res2 = integrate(V[1, 1] * conj(V[2, 2]), dDiagUnitary(d))
println("Result: ", res2)
println("Expected: 0")

# Example 3: Non-diagonal entries
# Diagonal unitary matrices have V_ij = 0 for i != j.
println("\n--- 4. Non-diagonal entries ---")
println("Computing E[ |V_12|^2 ]...")
res3 = integrate(abs(V[1, 2])^2, dDiagUnitary(d))
println("Result: ", res3)
println("Expected: 0")

# Example 4: Matrix Integration
println("\n--- 5. Matrix Integration over d=3 ---")
println("Integrating V * V' (should be Identity)")
res_V = integrate(V * V', dDiagUnitary(3))
display(res_V)

println("\nDone.")
