# examples/21_real_symbolic_dimension.jl
using IntU
using Symbolics
using LinearAlgebra

# Define symbolic dimension
@variables d

println("=== Symbolic Dimension Integration for Real Groups ===")

# --- Orthogonal Group O(d) ---
println("\n--- Orthogonal Group O(d) ---")
@symbolic_dimension O[1:d, 1:d]
mO = dO(O, d)

# 1. Integrate O_11^2
println("\n1. Integrating O[1,1]^2...")
expr_o1 = O[1, 1]^2
res_o1 = integrate(expr_o1, mO)
println("Result: $res_o1 (Expected: 1/d)")

# 2. Integrate O_11^4
println("\n2. Integrating O[1,1]^4...")
expr_o2 = O[1, 1]^4
res_o2 = integrate(expr_o2, mO)
# Theoretical: 3 / (d(d+2))
println("Result: $(Symbolics.simplify(res_o2)) (Expected: 3 / (d*(d+2)))")

# 3. Integrate O_11^2 * O_22^2
println("\n3. Integrating O[1,1]^2 * O[2,2]^2...")
expr_o3 = O[1, 1]^2 * O[2, 2]^2
res_o3 = integrate(expr_o3, mO)
# Theoretical: (d+1) / (d(d-1)(d+2))
println("Result: $(Symbolics.simplify(res_o3)) (Expected: (d+1) / (d*(d-1)*(d+2)))")


# --- Symplectic Group Sp(d) ---
println("\n--- Symplectic Group Sp(d) ---")
@symbolic_dimension S[1:d, 1:d]
mSp = dSp(S, d)

# Note on Sp(d) and symbolic dimensions:
# Integration over the symplectic group involves the symplectic form J. 
# For entries like S[1,1], the result depends on the pairing of indices (e.g., 1 with d/2 + 1).
# When d is symbolic, these pairings are not integer-constant, so many fixed-index 
# integrals will evaluate to 0 unless the indices are explicitly paired.

println("\n1. Integrating |S[1,1]|^2 with symbolic d...")
# Since J(1,1) = 0, any term involving only S_11 and its conjugate (mapped to the same)
# might vanish depending on the exact index matching.
expr_s1 = abs(S[1, 1])^2
res_s1 = integrate(expr_s1, mSp)
println("Result: $res_s1")
println("(Note: For fixed indices and symbolic d, Sp(d) integrals often vanish due to J_{ij} = 0)")

# However, the Weingarten functions themselves support symbolic d:
println("\n2. Symplectic Weingarten value for symbolic d")
p1 = [(1, 2)]
p2 = [(1, 2)]
wg_sp = weingarten_symplectic_val(p1, p2, d)
println("Wg^Sp(p1, p2, d) = $(Symbolics.simplify(wg_sp))")
