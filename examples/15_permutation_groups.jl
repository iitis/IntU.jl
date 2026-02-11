using IntU
using Symbolics
using LinearAlgebra

# 1. Basic Integration over the Symmetric Group S_d
# We can use a fixed dimension and concrete variables
println("1. Integration over S_d (d=4)")
d_val = 4
@variables P[1:4, 1:4]
measure = dPerm(P, d_val)

# E[P_11] = 1/d
expr1 = P[1, 1]
println("Integrating: ", expr1)
result1 = integrate(expr1, measure)
println("Result: ", result1)
println("Expected: ", 1//d_val)

# E[P_11 * P_22] = 1 / (d * (d-1))
expr2 = P[1, 1] * P[2, 2]
println("\nIntegrating: ", expr2)
result2 = integrate(expr2, measure)
println("Result: ", result2)
println("Expected: ", 1//(d_val * (d_val - 1)))

# E[P_11 * P_12] = 0 (same row, different column is impossible for permutation matrix)
expr3 = P[1, 1] * P[1, 2]
println("\nIntegrating: ", expr3)
result3 = integrate(expr3, measure)
println("Result: ", result3)
println("Expected: 0")

# 2. Integration with Symbolic Dimension
println("\n2. Integration with Symbolic Dimension")
@variables d
# For symbolic arrays, we can use a small fixed size for indices
P_sym = [Symbolics.variable(:P, i, j) for i = 1:2, j = 1:2]
measure_sym = dPerm(P_sym, d)

expr4 = P_sym[1, 1] * P_sym[2, 2]
println("Integrating: ", expr4)
result4 = integrate(expr4, measure_sym)
println("Result: ", result4)
println("Simplified: ", Symbolics.simplify(result4))

# 3. Centered Permutation Group
# Matrices Y = P - J/d
println("\n3. Centered Permutation Group")
@variables Y[1:4, 1:4]
m_centered = dCPerm(Y, d_val)

# E[Y_11] = 0
expr5 = Y[1, 1]
println("Integrating: ", expr5)
result5 = integrate(expr5, m_centered)
println("Result: ", result5)

# E[Y_11^2] = (d-1)/d^2
expr6 = Y[1, 1]^2
println("\nIntegrating: ", expr6)
result6 = integrate(expr6, m_centered)
println("Result: ", result6)
println("Simplified: ", Symbolics.simplify(result6))
println("Expected (d=4): ", (d_val - 1) / d_val^2)

println("Expected (d=4): ", (d_val - 1) / d_val^2)

# 3b. Matrix Integration (Permutation Group)
println("\n3b. Matrix Integration")
println("Integrating P * P^T (should be Identity)")
# P is defined as simple Symbolic array [P_ij], so it's already a matrix.
# P * P' should integrate to I
# Convert to standard Matrix{Num} to ensure generic integration works
P_mat = collect(P)
res_P = integrate(P_mat * P_mat', measure)

println("Result[1,1]: ", res_P[1, 1])
println("Expected: 1")
println("Result[1,2]: ", res_P[1, 2])
println("Expected: 0")

# 4. Symbolic Traces with Symbolics.jl arrays
println("\n4. Symbolic Traces with Symbolics.jl arrays")
@variables A[1:2, 1:2]
# E[tr(P * A)] = sum(A) / d
expr7 = Symbolics.scalarize(IntU.tr(P_sym * A))
println("Integrating: tr(P * A)")
result7 = integrate(expr7, measure_sym)
println("Result: ", Symbolics.simplify(result7))
# Expected: (A[1,1] + A[1,2] + A[2,1] + A[2,2]) / d
