using IntU
using Symbolics
using LinearAlgebra

println("=== Diagonal Unitary Integration (Torus Group) ===\n")

@variables d
D = SymbolicMatrix(:D, :DiagUnitary)

println("--- 1. Define Measure ---")
println("Using measure dDiagUnitary(d)")

# Example 1: Simple Phase Average
# E[ |D_11|^2 ] = 1
println("\n--- 2. Single Entry Moment ---")
println("Computing E[ |D_11|^2 ]...")
res1 = integrate(abs(D[1, 1])^2, dDiagUnitary(d))
println("Result: ", res1)
println("Expected: 1")

# Example 2: Independent Phases
# E[ D_11 * D_22^* ] = 0
println("\n--- 3. Independent Phases ---")
println("Computing E[ D_11 * conj(D_22) ]...")
res2 = integrate(D[1, 1] * conj(D[2, 2]), dDiagUnitary(d))
println("Result: ", res2)
println("Expected: 0")

# Example 3: Non-diagonal entries
# Diagonal unitary matrices have D_ij = 0 for i != j.
println("\n--- 4. Non-diagonal entries ---")
println("Computing E[ |D_12|^2 ]...")
res3 = integrate(abs(D[1, 2])^2, dDiagUnitary(d))
println("Result: ", res3)
println("Expected: 0")

# Example 4: Matrix Integration
println("\n--- 5. Matrix Integration over d=3 ---")
println("Integrating D * D' (should be Identity)")
res_D = integrate(D * D', dDiagUnitary(3))
display(res_D)

println("\nDone.")
