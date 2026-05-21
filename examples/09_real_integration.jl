using IntegrateUnitary
using Symbolics
using LinearAlgebra

println("=== Orthogonal and Symplectic Group Integration ===")

@variables d

# --- Orthogonal Group O(d) ---
println("\n--- Orthogonal Group O(d) ---")
# O will be auto-declared by @integrate
# Note: for O(d), the measure is dO(d) and standard symbol is O

println("1. Integrating O[1,1]^2 over O(d)")
res1 = @integrate O[1, 1]^2 dO(d)
println("Result: ", res1, " (Expected: 1/d)")

println("\n2. Integrating O[1,1]^4 over O(d)")
res2 = @integrate O[1, 1]^4 dO(d)
println("Result: ", res2, " (Expected: 3/(d(d+2)))")

println("\n3. Orthogonality check: sum_k O[1,k]*O[1,k] over O(3)")
# Note: for fixed dimension demonstration
res3 = @integrate sum(O[1, k]^2 for k = 1:3) dO(3)
println("Result: ", res3, " (Expected: 1)")

println("\n4. Matrix Integration Example over O(2)")
println("Integrating O * O^T (should be Identity)")
# Use factory function for manual matrix creation if needed outside @integrate
O_mat = symbolic_orthogonal(:O, 2)
res_mat = integrate(O_mat * O_mat', dO(2))
println("Result matrix:\n", res_mat)


# --- Symplectic Group Sp(d) ---
println("\n--- Symplectic Group Sp(d) ---")
println("(Note: Sp(d) integration usually requires fixed even d for explicit contractions)")

# Measure dSp(Sp, d)
d_sp = 2

println("1. Integrating Sp[1,1]*Sp[2,2] over Sp(2)")
res_sp1 = @integrate Sp[1, 1] * Sp[2, 2] dSp(d_sp)
println("Result: ", res_sp1, " (Expected: 0.5)")

println("\n2. Integrating Sp[1,2]*Sp[2,1] over Sp(2)")
res_sp2 = @integrate Sp[1, 2] * Sp[2, 1] dSp(d_sp)
println("Result: ", res_sp2, " (Expected: -0.5)")

println("\n3. Integrating |Sp[1,1]|^2 over Sp(2)")
res_sp3 = @integrate abs(Sp[1, 1])^2 dSp(2)
println("Result: ", res_sp3, " (Expected: 0.5)")
