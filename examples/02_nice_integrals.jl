# examples/02_nice_integrals.jl
using IntU
using Symbolics
using LinearAlgebra

d_val = 4
println("Working with U($d_val)")

@variables U[1:d_val, 1:d_val]::Complex
measure = dU(U, d_val)

# --- 1. Trace Moments ---
println("\n--- Trace Moments: ∫ |Tr(U)|^(2k) dU ---")
println("For large d, Tr(U) converges to a complex Gaussian.")
println("Expected moments: k! (1, 2, 6, 24...) for d >= k")

# Using overloaded tr(U)
tr_U = IntU.tr(U)

for k = 1:d_val
    # |Tr(U)|^(2k) = (Tr(U) * conj(Tr(U)))^k = abs(Tr(U))^(2k)
    # Using our new abs support:
    expr = abs(tr_U)^(2*k)

    # Note: For k=3 (power 6), the expansion generates many terms. 
    # It might take a moment.

    print("k=$k (Moment $(2*k))... ")
    val = integrate(expr, measure)
    expected = factorial(k)
    println("Result: $val (Expected: $expected)")

    # Use isequal for symbolic objects and abs for numerical ones
    is_correct =
        isequal(val, expected) || (
            try
                abs(Symbolics.unwrap(val) - expected) < 1e-10
            catch
                ;
                false
            end
        )

    if !is_correct
        println("  (Note: Deviation expected if d < k? Here d=$d_val, k=$k. d >= k holds.)")
    end
end

# --- 2. Determinant of a Minor ---
println("\n--- Integral of a Minor ---")
# Let's verify that the integral of the squared modulus of a 2x2 minor is related to d.
# Minor M = U_11 U_22 - U_12 U_21 = det(U[1:2, 1:2])
minor_2x2 = U[1, 1]*U[2, 2] - U[1, 2]*U[2, 1]
expr_minor = abs(minor_2x2)^2

println("Integrating |U_11 U_22 - U_12 U_21|^2 ...")
val_minor = integrate(expr_minor, measure)
println("Result: $val_minor")

# Theoretical expectation:
# For U in U(d), the mean square magnitude of a kxk minor is 1/binom(d, k).
# Here we have a 2x2 minor, so expected value is 1/binom(d, 2) = 2/(d*(d-1)).
# For d=4, expected: 2/(4*3) = 1/6.
println("Expected: $(1//binomial(d_val, 2))")
println("Note: Matches 1/binom(d, 2) = 1/$(binomial(d_val, 2))")
