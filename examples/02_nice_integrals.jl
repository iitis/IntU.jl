using IntU
using Symbolics
using LinearAlgebra

# Use a symbolic dimension for general formulas
@variables d
U = SymbolicMatrix(:U, :U)

# --- 1. Trace Moments ---
println("\n--- Trace Moments: ∫ |Tr(U)|^(2k) dU ---")
println("Pure trace moments |tr(U)|^{2k} depend on d as a step function,")
println("so they require a concrete integer dimension (k > 1).")
println("For d >= k the result equals k!.")

# k=1 works with symbolic d
expr1 = abs2(tr_lazy(U))
val1 = integrate(expr1, dU(d))
println("k=1 (symbolic d): $val1")

# k >= 2 requires concrete d
U10 = SymbolicMatrix(:U, :U, 10)
for k = 1:3
    expr = abs2(tr_lazy(U10))^k
    print("k=$k (d=10)... ")
    val = integrate(expr, dU(10))
    println("Result: $val (expected: $(factorial(k)))")
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
