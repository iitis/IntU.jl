using IntegrateUnitary
using Symbolics
using LinearAlgebra

# 1. Use the @integrate macro
# This automatically declares variables and matrices
println("Working with U(d) symbolically using @integrate")

# 2. Example: Normalization of a matrix element
# Integral of |U_{11}|^2
println("\nIntegrating: |U_11|^2")
result1 = @integrate abs2(U[1, 1]) dU(d)
println("Result: ", result1)
println("Expected: 1/d")

# 3. Example: Correlated elements
# Integral of |U_{11} U_{22}|^2
println("\nIntegrating: abs2(U[1, 1] * U[2, 2])")
result2 = @integrate abs2(U[1, 1] * U[2, 2]) dU(d)
println("Result: ", result2)
println("Expected: 1/(d^2 - 1)")

# 4. Example: Matrix Integration
println("\n4. Example: Matrix Integration using factory functions")
println("Integrating U * U' over U(2) (should be Identity)")
U = symbolic_unitary(:U, 2)
result_mat = integrate(U * U', dU(2))
println("Result:\n", result_mat)
println("Result is Identity? ", result_mat == I)

# 5. Example: Substituting values
println("\n5. Example: Using evaluate() to substitute d=5")
val_d5 = evaluate(result1, d => 5)
println("<|U_11|^2> over U(5) = ", val_d5)
