# Example 22: Streamlined UX Interface
# ===================================
# This example demonstrates the new, more ergonomic interface for 
# symbolic integration in IntegrateUnitary.jl.

using IntegrateUnitary
using Symbolics
using LinearAlgebra
@variables d

println("=== Streamlined UX Interface Demo ===")

# 1. The @integrate macro
# -----------------------
# The @integrate macro reduces boilerplate by automatically declaring 
# symbolic variables and matrices.

println("\n1. Using @integrate macro:")
# You don't need to declare d or U!
res1 = @integrate abs(U[1, 1])^2 dU(d)
println("Integral of |U_11|^2 over dU(d): ", res1)

# It also works with numerical dimensions:
res2 = @integrate abs2(U[1, 1]) dU(3)
println("Integral of |U_11|^2 over dU(3): ", res2)

# Multiple matrices? No problem. Unknown symbols are treated as Constant matrices.
res3 = @integrate tr(U * A * U' * B) dU(d)
println("Integral of tr(U A U' B): ", res3)


# 2. Evaluation with evaluate()
# -----------------------------
# You can easily substitute symbolic values in the result.
println("\n2. Substituting values with evaluate():")
println("Substituting d => 10 in res1:")
val_sub = evaluate(res1, d => 10)
println("Result: ", val_sub)

# It also works with multiple substitutions:
println("\nSubstituting d => 4, tr(A) => 2, tr(B) => 3 in res3:")
val_sub2 = evaluate(res3, [d => 4, tr(A) => 2, tr(B) => 3])
println("Result: ", val_sub2)


# 3. Factory functions
# --------------------
# If you prefer manual declaration, new factory functions make it cleaner.
println("\n3. Using factory functions:")

U = symbolic_unitary(:U, d)        # Shorthand for SymbolicMatrix(:U, :U, d)
println("U: ", U)

O = symbolic_orthogonal(:O, d)     # Shorthand for SymbolicMatrix(:O, :O, d)
println("O: ", O)

psi = symbolic_pure_state(:psi, d) # Shorthand for SymbolicMatrix(:psi, :psi, (d, 1))
println("psi: ", psi)

# These work seamlessly with the standard integrate() function:
res_o = integrate(O[1, 1]^4, dO(d))
println("Integral of O_11^4: ", res_o)

println("\nDemo completed.")
