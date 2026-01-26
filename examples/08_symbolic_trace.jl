# Example 08: Symbolic Trace Logic
# ==================================
# This example demonstrates how to use the symbolic trace logic to integrate
# traces of random unitary matrices and constant matrices.

using IntU
using Symbolics
using LinearAlgebra

# 1. Define symbolic variables
@variables d
# Create a Haar measure for unitary matrix U of dimension d
# We create a symbolic matrix variable :U.
U_var = SymbolicMatrix(:U, false, :U)
# HaarMeasure expects an array for U, but for symbolic trace logic we only need the dimension d.
# We pass a dummy 1x1 symbolic array to satisfy the constructor.
@variables u_dummy[1:1, 1:1]
measure = dU(u_dummy, d)

# 2. Define constant matrices
# By default, SymbolicMatrix creates a constant matrix
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)
C = SymbolicMatrix(:C)

println("Integrating symbolic traces over Haar measure U(d)...")

# Case 1: Tr(U A U' B)
# --------------------
# We want to compute ∫ dU Tr(U A U' B).
# Construct the trace lazily:
trace_expr = tr_lazy(U_var * A * U_var' * B)
println("\nExpression 1: ", trace_expr)

result = integrate(trace_expr, measure)
println("Result: ", result)
# Expected: tr(A) * tr(B) / d

# Case 2: Tr(U A U' B U C U')
# ---------------------------
# Integration of a longer product.
trace_expr2 = tr_lazy(U_var * A * U_var' * B * U_var * C * U_var')
println("\nExpression 2: ", trace_expr2)

result2 = integrate(trace_expr2, measure)
println("Result: ", result2)
# The result will involve Weingarten functions of 2nd order (since U appears twice).

# Case 3: Pure State Integration
# ------------------------------
# We can also handle integration over pure states if we formulate them as
# projectors or traces involving U |0><0| U'.
# But currently symbolic_trace.jl focuses on full Unitary matrix U.

# Case 4: Higher moments
# ----------------------
# Let's verify E[Tr(U A) Tr(U' B)]?
# Currently `integrate` handles a SINGLE trace.
# Product of traces is not yet directly supported in `integrate` API for single calls,
# but can be handled if we expand the theory or use tensor product scaling.
# Supported: Single trace of product.

println("\nDone.")
