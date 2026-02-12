using IntU
using Symbolics
using LinearAlgebra

# 1. Define the dimension (can be symbolic or numeric)
@variables d
U = SymbolicMatrix(:U)

println("Working with U(d) symbolically")

# 2. Example: Normalization of a matrix element
# Integral of |U_{11}|^2
println("\nIntegrating: abs(U[1, 1])^2")
result1 = @integrate abs(U[1, 1])^2 dU(d)
println("Result: ", result1)
println("Expected: 1/d")

# 3. Example: Correlated elements
# Integral of |U_{11} U_{22}|^2
println("\nIntegrating: abs(U[1, 1] * U[2, 2])^2")
result2 = @integrate abs(U[1, 1] * U[2, 2])^2 dU(d)
println("Result: ", result2)
println("Expected: 1/(d^2 - 1)")

# 5. Example: Matrix Integration
# You can integrate matrix-valued expressions directly.
println("\n5. Example: Matrix Integration (Fixed Dimension for demonstration)")
println("Integrating U * U' over U(2) (should be Identity)")
result_mat = @integrate U * U' dU(2)
println("Result:\n", result_mat)
println("Result is Identity? ", result_mat == I)

# 6. Example: Mixed Symbolic/Numeric
println("\n6. Example: Fixed dimension d=3")
res3 = @integrate abs(U[1,1])^2 dU(3)
println("<|U_11|^2> over U(3) = ", res3)
