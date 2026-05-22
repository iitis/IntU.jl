using IntegrateUnitary
using Symbolics
using LinearAlgebra

# 1. Basic Integration over the Symmetric Group S_d
println("1. Integration over S_d (d symbolic) using @integrate")

# E[P_11] = 1/d
println("Integrating: P[1, 1]")
result1 = @integrate P[1, 1] dPerm(d)
println("Result: ", result1)
println("Expected: 1/d")

# E[P_11 * P_22] = 1 / (d * (d-1))
println("\nIntegrating: P[1, 1] * P[2, 2]")
result2 = @integrate P[1, 1] * P[2, 2] dPerm(d)
println("Result: ", result2)
println("Expected: 1/(d*(d-1))")

# E[P_11 * P_12] = 0
println("\nIntegrating: P[1, 1] * P[1, 2]")
result3 = @integrate P[1, 1] * P[1, 2] dPerm(d)
println("Result: ", result3)
println("Expected: 0")


# 2. Centered Permutation Group
# Matrices Y = P - J/d
println("\n2. Centered Permutation Group")

# E[Y_11] = 0
println("Integrating: Y[1, 1]")
result5 = @integrate Y[1, 1] dCPerm(d)
println("Result: ", result5)

# E[Y_11^2] = (d-1)/d^2
println("\nIntegrating: Y[1, 1]^2")
result6 = @integrate Y[1, 1]^2 dCPerm(d)
println("Result: ", result6)
println("Simplified: ", Symbolics.simplify(result6))


# 3. Matrix Integration
println("\n3. Matrix Integration over S_3")
println("Integrating P * P^T (should be Identity)")
P_mat = symbolic_permutation(:P, 3)
res_P = integrate(P_mat * P_mat', dPerm(3))
display(res_P)

# 4. Symbolic Traces
println("\n4. Symbolic Traces over S_d")
# E[tr(P * A)] = sum(A) / d
println("Integrating: tr(P * A)")
# Use @integrate, A is auto-declared as Constant, P as Perm
result7 = @integrate tr(P * A) dPerm(d)
println("Result: ", result7)
println("Expected: Sum(A_ij) / d")
