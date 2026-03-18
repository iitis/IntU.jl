using IntU
using Symbolics
using LinearAlgebra

# Use a symbolic dimension for general formulas
@variables d
U = SymbolicMatrix(:U, :U)

# --- 1. Trace Moments ---
println("\n--- Trace Moments: ∫ |Tr(U)|^(2k) dU ---")
println("For large d, Tr(U) converges to a complex Gaussian.")
println("Expected moments: k! (1, 2, 6, 24...) for d >= k")

for k = 1:3
    # |Tr(U)|^(2k) = (Tr(U)Tr(U'))^k
    expr = abs2(tr_lazy(U))^k

    print("k=$k (Moment $(2*k))... ")
    val = integrate(expr, dU(d))
    expected = factorial(k)
    println("Result: $val")

    # We can use symbolic substitution to check specific values
    val_at_d4 = Symbolics.substitute(val, Dict(d => 4))
    println("  At d=4: $val_at_d4 (Expected: $expected)")
end

# --- 2. Determinant of a Minor ---
println("\n--- Integral of a Minor ---")
# Minor M = U_11 U_22 - U_12 U_21 = det(U[1:2, 1:2])
minor_2x2 = U[1, 1]*U[2, 2] - U[1, 2]*U[2, 1]
expr_minor = minor_2x2 * conj(minor_2x2)

println("Integrating |U_11 U_22 - U_12 U_21|^2 over U(d) ...")
val_minor = integrate(expr_minor, dU(d))
println("Result: $val_minor")

# Theoretical expectation:
# For U in U(d), the mean square magnitude of a kxk minor is 1/binom(d, k).
# Here we have a 2x2 minor, so expected value is 1/binom(d, 2) = 2/(d*(d-1)).
expected_minor = 2 / (d * (d - 1))
println("Expected: $expected_minor")

# Check at d=4
val_minor_d4 = Symbolics.substitute(val_minor, Dict(d => 4))
println("Result at d=4: $val_minor_d4")
println("Expected at d=4: $(1//binomial(4, 2))")
