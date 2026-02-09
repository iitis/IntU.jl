using IntU
using Symbolics
using LinearAlgebra

# 1. Define the dimension and variables
d_val = 3
println("Working with U($d_val)")

@variables U[1:d_val, 1:d_val]::Complex
measure = dU(U, d_val)

# 2. Example: Normalization of a matrix element
# Integral of |U_{11}|^2
expr1 = abs(U[1, 1])^2
println("\nIntegrating: ", expr1)
result1 = integrate(expr1, measure)
println("Result: ", result1)
println("Expected: ", 1//d_val)

# 3. Example: Correlated elements
# Integral of |U_{11} U_{22}|^2
expr2 = abs(U[1, 1] * U[2, 2])^2
println("\nIntegrating: ", expr2)
result2 = integrate(expr2, measure)
println("Result: ", result2)
println("Expected: ", 1//(d_val^2 - 1))

# 5. Example: Matrix Integration
# You can integrate matrix-valued expressions directly.
# We use collect() to ensure we pass a standard Julia Matrix of symbolic numbers.
println("\n5. Example: Matrix Integration")
println("Integrating U * U' (should be Identity)")
expr_mat = collect(U * U')
result_mat = integrate(expr_mat, measure)
println("Result[1,1]: ", result_mat[1, 1])
println("Result is Identity? ", result_mat == I)

# 6. Example: General d formula (symbolic d is not fully supported yet,
# but we can show it matches the theoretical value for a specific d)
println("\nNote: Theoretical value involves Weingarten function Wg(1^2, d) = 1/(d^2-1)")
