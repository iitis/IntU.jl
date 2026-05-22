# Stiefel Manifold Integration Example
# ------------------------------------
# This example demonstrates integration over the Stiefel manifold V_k(C^d),
# which consists of d x k matrices with orthonormal columns (V'V = I_k).

using IntegrateUnitary
using Symbolics
using LinearAlgebra

println("=== Stiefel Manifold Integration Example ===")

@variables d
k = 2

println("\n--- Defining Variables ---")
println("Dimension d: symbolic")
println("Stiefel dimension k: $k")

V = SymbolicMatrix(:V, :V)

# 2. Normalization check
# E[V_{1,1} * conj(V_{1,1})] should be 1/d
println("\nCalculating E[|V_{1,1}|^2]...")
res1 = integrate(abs(V[1, 1])^2, dStiefel(d, k))
println("Result: ", Symbolics.simplify(res1))
println("Expected: 1/d")

# E[V_{1,1} * conj(V_{1,2})] should be 0 (orthogonality of columns)
println("\nCalculating E[V_{1,1} * conj(V_{1,2})]...")
res2 = integrate(V[1, 1] * conj(V[1, 2]), dStiefel(d, k))
println("Result: ", Symbolics.simplify(res2))
println("Expected: 0")

# 3. Higher moments
# E[|V_{1,1}|^2 * |V_{1,2}|^2]
println("\nCalculating E[|V_{1,1}|^2 * |V_{1,2}|^2]...")
res3 = integrate(abs(V[1, 1])^2 * abs(V[1, 2])^2, dStiefel(d, k))
println("Result: ", Symbolics.simplify(res3))
println("Expected: 1/(d*(d+1))")

# 4. Asymptotic expansion
println("\n--- Asymptotic Expansion ---")
asymp = asymptotic(abs(V[1, 1])^2 * abs(V[1, 2])^2, dStiefel(d, k), 2)
println("Asymptotic expansion of E[|V_{1,1}|^2 * |V_{1,2}|^2]:")
println(asymp)

# 5. Matrix Integration
println("\n--- Matrix Integration over V_2(C^3) ---")
println("Integrating V' * V (should be Identity I_k)")
# V is d x k. V' * V is k x k (here 2x2).
res_Id = integrate(V' * V, dStiefel(3, 2))
println("Result:")
display(res_Id)
println("Expected: I(2)")
