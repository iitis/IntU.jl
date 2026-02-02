using IntU
using Symbolics

# Define variables
N = 3
@variables U[1:N, 1:N]::Complex

println("--- Unitary t-Design Example ---")

# Create a 2-design
# This measure will behave like Haar measure for polynomials of degree <= 2
# and throw an error for higher degrees.
design2 = dDesign(U, N, 2)
println("Created Unitary 2-Design for N=$N")

# Example 1: Integrating |u11|^2 (Degree 1)
println("\n1. Integrating |U[1,1]|^2 (Degree 1)")
expr1 = abs(U[1,1])^2
res1 = integrate(expr1, design2)
println("Result: $res1")
println("Expected (Haar): $(1/N)")

# Example 2: Integrating |u11 u22|^2 (Degree 2)
println("\n2. Integrating |U[1,1] U[2,2]|^2 (Degree 2)")
expr2 = abs(U[1,1] * U[2,2])^2
res2 = integrate(expr2, design2)
println("Result: $res2")
println("Expected (Haar): $(1/(N^2 - 1))")

# Example 3: Attempting to integrate |u11|^6 (Degree 3)
println("\n3. Attempting to integrate |U[1,1]|^6 (Degree 3)")
expr3 = abs(U[1,1])^6
try
    res3 = integrate(expr3, design2)
    println("Result: $res3")
catch e
    println("Caught expected error: ", e)
end

println("\n--- Comparison with full Haar Measure ---")
@variables d_sym
@symbolic_dimension U_haar[1:d_sym, 1:d_sym]
haar = dU(U_haar)
res3_haar = integrate(abs(U_haar[1,1])^6, haar)
println("Haar Result for |U[1,1]|^6: $res3_haar")
