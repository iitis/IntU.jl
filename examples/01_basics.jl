# examples/01_basics.jl
using IntU
using Symbolics

# 1. Define the dimension and variables
d_val = 3
println("Working with U($d_val)")

@variables U[1:d_val, 1:d_val]::Complex
measure = dU(U, d_val)

# 2. Example: Normalization of a matrix element
# Integral of |U_{11}|^2
expr1 = abs(U[1,1])^2
println("\nIntegrating: ", expr1)
result1 = integrate(expr1, measure)
println("Result: ", result1)
println("Expected: ", 1//d_val)

# 3. Example: Correlated elements
# Integral of |U_{11} U_{22}|^2
expr2 = abs(U[1,1] * U[2,2])^2
println("\nIntegrating: ", expr2)
result2 = integrate(expr2, measure)
println("Result: ", result2)
println("Expected: ", 1//(d_val^2 - 1))

# 4. Example: General d formula (symbolic d is not fully supported yet, 
# but we can show it matches the theoretical value for a specific d)
println("\nNote: Theoretical value involves Weingarten function Wg(1^2, d) = 1/(d^2-1)")
