# Example 08: Symbolic Trace Logic
# ==================================
# This example demonstrates how to use the symbolic trace logic to integrate
# traces of random unitary matrices and constant matrices.

using IntU
using Symbolics
using LinearAlgebra

# 1. Define symbolic variables
@variables d
# Create a Haar-random SymbolicMatrix U of dimension d
U = SymbolicMatrix(:U, :U)

# 2. Define constant matrices using SymbolicMatrix
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)
C = SymbolicMatrix(:C)

println("Integrating symbolic traces over Haar measure U(d)...")

# Case 1: Tr(U A U' B)
# --------------------
# We want to compute ∫ dU Tr(U A U' B).
# The new interface handles traces of SymbolicMatrix products directly.
# tr(U*A*U'*B) creates a LazyTrace object.
println("\nExpression 1: tr(U * A * U' * B)")

result = integrate(tr(U * A * U' * B), dU(d))
println("Result: ", result)
# Expected: tr(A) * tr(B) / d

# Case 2: Tr(U A U' B U C U')
# ---------------------------
# Integration of a longer product.
println("\nExpression 2: tr(U * A * U' * B * U * C * U')")

result2 = integrate(tr(U * A * U' * B * U * C * U'), dU(d))
println("Result: ", Symbolics.simplify(result2))
# The result will involve Weingarten functions of 2nd order.


# Case 3: Product of traces
# -------------------------
# We can integrate products of traces, e.g. E[Tr(U A) Tr(U' B)].
println("\nExpression 3: tr(U * A) * tr(U' * B)")

# Multiplication of Symbolic Traces also produces a LazyTrace
expr_prod = tr(U * A) * tr(U' * B)

result_prod = integrate(expr_prod, dU(d))
println("Result: ", Symbolics.simplify(result_prod))
# Expected: tr(A B) / d

# Case 4: Powers of traces
println("\nExpression 4: tr(U)^k * tr(U')^k")
for k in 1:2
    expr_k = tr(U)^k * tr(U')^k
    res_k = integrate(expr_k, dU(d))
    println("k=$k result: ", res_k)
end

println("\nDone.")
