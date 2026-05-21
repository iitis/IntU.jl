using IntU
using Symbolics

# Define variables
N = 3
U = SymbolicMatrix(:U, :U, N)

println("--- Unitary t-Design Example ---")

# Create a 2-design
k = 2
design2 = dDesign(N, k)
println("Created Unitary 2-Design for N=$N")

# Example 1: Integrating |u11|^2 (Degree 1)
println("\n1. Integrating |U[1,1]|^2 (Degree 1)")
res1 = integrate(abs(U[1, 1])^2, design2)
println("Result: $res1")
println("Expected (Haar): $(1/N)")

# Example 2: Integrating |u11 u22|^2 (Degree 2)
println("\n2. Integrating |U[1,1] U[2,2]|^2 (Degree 2)")
res2 = integrate(abs(U[1, 1] * U[2, 2])^2, design2)
println("Result: $res2")
println("Expected (Haar): $(1/(N^2 - 1))")

# Example 3: Attempting to integrate |u11|^6 (Degree 3)
println("\n3. Attempting to integrate |U[1,1]|^6 (Degree 3)")
# res3 = integrate(abs(U[1, 1])^6, design2)
# println("Result: $res3")
println("Note: Integrating degree 3 on a 2-design throws an error as expected.")

println("\n--- Comparison with full Haar Measure ---")
@variables d_sym
res3_haar = integrate(abs(U[1, 1])^6, dU(d_sym))
println("Haar Result for |U[1,1]|^6 over U(d): $res3_haar")
println("Simplified: ", Symbolics.simplify(res3_haar))
