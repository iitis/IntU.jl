using IntU
using LinearAlgebra
using Symbolics

println("=== HCIZ Integral Examples ===\n")

# Example 1: Basic Numeric Integration for d=2
println("1. Basic Numeric Integration (d=2)")
A = diagm([1.0, 2.0])
B = diagm([0.5, 1.5])

# The HCIZ integral calculates ∫ dU exp(Tr(A U B U'))
res2 = hciz(A, B)
println("A = diag([1.0, 2.0]), B = diag([0.5, 1.5])")
println("hciz(A, B) = ", res2)
# Analytical check: (exp(1*0.5 + 2*1.5) - exp(1*1.5 + 2*0.5)) / ((1-2)*(0.5-1.5))
# = (exp(3.5) - exp(2.5)) / 1 ≈ 33.115 - 12.182 = 20.933
println()

# Example 2: Larger Dimension (d=5)
println("2. Larger Dimension (d=5)")
d = 5
A5 = diagm(collect(1.0:d))
B5 = diagm(reverse(collect(1.0:d)))
res5 = hciz(A5, B5)
println("d=$d, A=diag(1:5), B=diag(5:-1:1)")
println("hciz(A, B) = ", res5)
println()

# Example 3: Symbolic Integration
println("3. Symbolic Integration")
@variables a1 a2 b1 b2
# We can pass vectors of eigenvalues directly
sym_res = hciz([a1, a2], [b1, b2])
println("Symbolic d=2 result:")
display(sym_res)
println("\n")

# Example 4: Handling Degeneracies
println("4. Handling Degenerate Eigenvalues (Numerical Stability)")
# Values are exactly same, which would normally cause 0/0
Ad = diagm([1.0, 1.0])
Bd = diagm([0.5, 1.5])
res_deg = hciz(Ad, Bd)
println("Degenerate A=diag([1, 1]): ", res_deg)
# For a1=a2, the formula should limit to exp(a1*b1 + a1*b2) * (something)
# Our implementation uses a tiny perturbation to handle this numerically.
println()

# Example 5: Comparison with Haar Integration (Verification)
println("5. Verification against Weingarten integration (d=2)")
# exp(Tr(A U B U')) ≈ 1 + Tr(A U B U') + 0.5 * Tr(A U B U')^2 + ...
# We'll check the first non-trivial moment: ∫ dU Tr(A U B U')
# For d=2, ∫ dU U_ij U*_kl = δ_ik δ_jl / 2
# ∫ dU Tr(A U B U') = Σ_ijkl A_ii B_kk ∫ dU U_ik U*_ik = Σ_ik A_ii B_kk / d = Tr(A)Tr(B) / d

A_test = diagm([1.0, 0.5])
B_test = diagm([2.0, 1.0])
@variables U[1:2, 1:2]
measure = dU(U, 2)

# Avoiding matrix-symbolic multiplication ambiguity by performing multiplication in steps
expr_matmul = (A_test * U) * (B_test * adjoint(U))
term1 = integrate(IntU.tr(collect(expr_matmul)), measure)

# HCIZ should match if we take the linear part? No, HCIZ is the full exponential.
# But we can check small t: hciz(t*A, B) ≈ 1 + t*Tr(A)Tr(B)/d
t = 1e-4
res_hciz_small = hciz(t * A_test, B_test)
expected_linear = 1 + t * IntU.tr(A_test) * IntU.tr(B_test) / 2
println("hciz(ε*A, B): ", res_hciz_small)
println("1 + ε*Tr(A)Tr(B)/d: ", expected_linear)
println("Relative difference: ", (res_hciz_small - expected_linear) / expected_linear)
println()

# Example 6: Non-Diagonal Symbolic Matrices (d=2)
println("6. Non-Diagonal Symbolic Matrices (d=2)")
# Since LinearAlgebra.eigen doesn't support symbolic matrices directly,
# we compute the eigenvalues for the symbolic case and pass them to hciz.
@variables x y
A_non_diag = [0 x; x 0] # Eigenvalues: x, -x
B_non_diag = [0 y; y 0] # Eigenvalues: y, -y

# We pass the eigenvalue vectors
res_non_diag = hciz([x, -x], [y, -y])
println("A = [0 x; x 0], B = [0 y; y 0]")
println("hciz([x, -x], [y, -y]) = ")
display(res_non_diag)
println()
