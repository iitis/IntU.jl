# Example 08: Symbolic Trace Logic
# ==================================
# This example demonstrates how to use the symbolic trace logic to integrate
# traces of random unitary matrices and constant matrices.

using IntegrateUnitary
using Symbolics
using LinearAlgebra

println("Integrating symbolic traces over Haar measure U(d)...")

# Case 1: Tr(U A U' B)
# --------------------
# We want to compute ∫ dU Tr(U A U' B).
# The @integrate macro handles LazyTrace creation automatically.
println("\nExpression 1: tr(U * A * U' * B)")

result = @integrate tr(U * A * U' * B) dU(d)
println("Result: ", result)
# Expected: tr(A) * tr(B) / d

# Case 2: Tr(U A U' B U C U')
# ---------------------------
# Integration of a longer product.
println("\nExpression 2: tr(U * A * U' * B * U * C * U')")

result2 = @integrate tr(U * A * U' * B * U * C * U') dU(d)
println("Result: ", Symbolics.simplify(result2))

# Case 3: Product of traces
# -------------------------
println("\nExpression 3: tr(U * A) * tr(U' * B)")

result_prod = @integrate tr(U * A) * tr(U' * B) dU(d)
println("Result: ", Symbolics.simplify(result_prod))

# Case 4: Power and evaluation
println("\nExpression 4: evaluate tr(U) * tr(U') for d=2")
expr4 = @integrate tr(U) * tr(U') dU(d)
val4 = evaluate(expr4, d => 2)
println("k=1 result for d=2: ", val4)

println("\nDone.")
