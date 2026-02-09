# Stiefel Manifold Integration Example
# ------------------------------------
# This example demonstrates integration over the Stiefel manifold V_k(C^d),
# which consists of d x k matrices with orthonormal columns (V'V = I_k).

using IntU
using Symbolics
using LinearAlgebra

println("=== Stiefel Manifold Integration Example ===")

@variables d

# 1. Define variables
# We deal with a d x k matrix V.
# Let's consider k = 2
k = 2

println("\n--- Defining Variables ---")
println("Dimension d: symbolic")
println("Stiefel dimension k: $k")

# To use IntU, we can define a symbolic matrix V
# We use explicit Symbolics.variable to ensure clear variable names
V = [Symbolics.variable(Symbol("V_$(i)_$(j)"), T=Complex{Num}) for i=1:3, j=1:k]
println("V (first 3 rows shown):")
display(V)

# Note: Even though we defined V with 3 rows for display, the integration
# parameter `dim` can be symbolic `d`. The indices of V are treated as abstract
# indices 1..d during integration if we use the measure correctly/abstractly.
# However, for concrete matrix multiplication V'V, we need dimensions to match.
# In IntU, we often define just the elements we need for the polynomial.

# Let's define specific elements for a moment calculation
v11 = Symbolics.variable(:v11, T=Complex{Num}) # V_{1,1}
v12 = Symbolics.variable(:v12, T=Complex{Num}) # V_{1,2}
v21 = Symbolics.variable(:v21, T=Complex{Num}) # V_{2,1}

# We manually mock the lookup by creating a measure with a "virtual" matrix structure
# or we use the specific elements in a constructed array.
# The cleanest way for arbitrary polynomials is to use a Symbolic Unitary 
# representation but limited to k columns.
# Currently StiefelMeasure takes a matrix V.
# Let's use a 2x2 block of V for demonstration (indices 1..2, 1..2)
V_sub = [Symbolics.variable(Symbol("V_$(i)_$(j)"), T=Complex{Num}) for i=1:2, j=1:2]
measure = dStiefel(V_sub, d, k)

# 2. Normalization check
# E[V_{1,1} * conj(V_{1,1})] should be 1/d
poly1 = V_sub[1,1] * conj(V_sub[1,1])
println("\nCalculating E[|V_{1,1}|^2]...")
res1 = integrate(poly1, measure)
println("Result: ", Symbolics.simplify(res1))
println("Expected: ", 1/d)

# E[V_{1,1} * conj(V_{1,2})] should be 0 (orthogonality of columns)
poly2 = V_sub[1,1] * conj(V_sub[1,2])
println("\nCalculating E[V_{1,1} * conj(V_{1,2})]...")
res2 = integrate(poly2, measure)
println("Result: ", Symbolics.simplify(res2))
println("Expected: ", 0)

# 3. Higher moments
# E[|V_{1,1}|^2 * |V_{1,2}|^2]
poly3 = abs2(V_sub[1,1]) * abs2(V_sub[1,2])
println("\nCalculating E[|V_{1,1}|^2 * |V_{1,2}|^2]...")
res3 = integrate(poly3, measure)
println("Result: ", Symbolics.simplify(res3))
println("Expected: ", 1/(d*(d+1)))

# 4. Asymptotic expansion
println("\n--- Asymptotic Expansion ---")
asymp = asymptotic(poly3, measure, 2)
println("Asymptotic expansion of E[|V_{1,1}|^2 * |V_{1,2}|^2]:")
println(asymp)

# 5. Matrix Integration
println("\n--- Matrix Integration ---")
println("Integrating V' * V (should be Identity I_k)")
# V is d x k. V' * V is k x k (here 2x2).
# We defined V using explicit loop so it is a Matrix{Num} already.
res_Id = integrate(V_sub' * V_sub, measure)

println("Result:")
display(res_Id)
println("Expected: I(2)")
