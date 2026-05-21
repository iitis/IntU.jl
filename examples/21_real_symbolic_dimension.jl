using IntegrateUnitary
using Symbolics
using LinearAlgebra

# Define symbolic dimension
@variables d

println("=== Symbolic Dimension Integration for Real Groups ===")

# --- Orthogonal Group O(d) ---
println("\n--- Orthogonal Group O(d) ---")
O = SymbolicMatrix(:O, :O)

# 1. Integrate O_11^2
println("1. Integrating O[1,1]^2 over O(d)")
res_o1 = integrate(O[1, 1]^2, dO(d))
println("Result: $res_o1 (Expected: 1/d)")

# 2. Integrate O_11^4
println("\n2. Integrating O[1,1]^4 over O(d)")
res_o2 = integrate(O[1, 1]^4, dO(d))
# Theoretical: 3 / (d(d+2))
println("Result: $(Symbolics.simplify(res_o2)) (Expected: 3 / (d*(d+2)))")

# 3. Integrate O_11^2 * O_22^2
println("\n3. Integrating O[1,1]^2 * O[2,2]^2 over O(d)")
res_o3 = integrate(O[1, 1]^2 * O[2, 2]^2, dO(d))
# Theoretical: (d+1) / (d(d-1)(d+2))
println("Result: $(Symbolics.simplify(res_o3)) (Expected: (d+1) / (d*(d-1)*(d+2)))")


# --- Symplectic Group Sp(d) ---
println("\n--- Symplectic Group Sp(d) ---")
println("1. Integrating abstract symplectic Weingarten functions")
p1 = [(1, 2)]
p2 = [(1, 2)]
wg_sp = weingarten_symplectic_val(p1, p2, d)
println("Wg^Sp(p1, p2, d) = $(Symbolics.simplify(wg_sp))")

println("\nDone.")
