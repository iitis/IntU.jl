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

# Example 2: Larger Dimension (d=5)
println("\n2. Larger Dimension (d=5)")
d_val = 5
A5 = diagm(collect(1.0:d_val))
B5 = diagm(reverse(collect(1.0:d_val)))
res5 = hciz(A5, B5)
println("d=$d_val, A=diag(1:5), B=diag(5:-1:1)")
println("hciz(A, B) = ", res5)

# Example 3: Symbolic Integration
println("\n3. Symbolic Integration")
@variables a1 a2 b1 b2
# We can pass vectors of eigenvalues directly
sym_res = hciz([a1, a2], [b1, b2])
println("Symbolic d=2 result:")
display(sym_res)

# Example 4: Comparison with Weingarten integration (d=2)
println("\n4. Verification against Weingarten integration (d=2)")
# exp(Tr(A U B U')) ≈ 1 + Tr(A U B U') + ...
# ∫ dU Tr(A U B U') = Tr(A)Tr(B) / d

A_test = diagm([1.0, 0.5])
B_test = diagm([2.0, 1.0])
U = SymbolicMatrix(:U, :U, 2)

println("Integrating Tr(A_test * U * B_test * U') over dU(2)...")
integrand = tr(A_test * U * B_test * U')
result_tr = @integrate integrand dU(2)

println("Result: ", result_tr)
println("Expected (Tr(A)Tr(B)/d): ", tr(A_test)*tr(B_test)/2)

# HCIZ for small perturbation
t = 1e-4
res_hciz_small = hciz(t * A_test, B_test)
expected_linear = 1 + t * tr(A_test) * tr(B_test) / 2
println("\nhciz(ε*A, B): ", res_hciz_small)
println("1 + ε*Tr(A)Tr(B)/d: ", expected_linear)

println("\nDone.")
