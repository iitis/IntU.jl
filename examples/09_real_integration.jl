using IntU
using Symbolics
using LinearAlgebra

println("=== Orthogonal and Symplectic Group Integration ===")

@variables d

# --- Orthogonal Group O(d) ---
println("\n--- Orthogonal Group O(d) ---")
@variables O[1:3, 1:3]
mO = dO(O, d)

println("1. Integrating |O[1,1]|^2 (should be 1/d)")
expr1 = O[1, 1]^2
res1 = integrate(expr1, mO)
println("Result: ", res1)

println("\n2. Integrating |O[1,1]|^4 (should be 3/(d(d+2)))")
expr2 = O[1, 1]^4
res2 = integrate(expr2, mO)
println("Result: ", Symbolics.simplify(res2))

println("\n3. Orthogonality check: sum_k O[1,k]*O[1,k] (should be 1)")
expr3 = sum(O[1, k]*O[1, k] for k = 1:3)
# Note: integrating over O(d) usually implies d matches the matrix size indices if we sum over all of them?
# If we treat d as symbolic but sum k=1..3, we get 3/d.
# If we set d=3, we get 1.
res3 = integrate(expr3, mO)
println("Result (symbolic sum k=1..3): ", res3)
println("Value at d=3: ", Symbolics.substitute(res3, Dict(d => 3)))


println("\n4. Matrix Integration Example")
println("Integrating O * O^T (should be Identity)")
# We need to collect symbolic array to standard Matrix{Num} to ensure element-wise integration
O_mat = collect(O)
res_mat = integrate(O_mat * O_mat', mO)

# Check first element
println("Result[1,1]: ", res_mat[1, 1])
println("Expected: 1")


# --- Symplectic Group Sp(d) ---
println("\n--- Symplectic Group Sp(d) ---")
println("(Note: Sp(d) integration requires d to be even)")

@variables S[1:2, 1:2]
# Note: For explicit index calculations involving conjugates, we currently require numeric d to evaluate J contractions.
d_val = 2
@variables S[1:2, 1:2]::Complex
mS = dSp(S, d_val)

println("1. Defining Measure dSp(S, $d_val) with complex variables")
println("Measure: ", mS)

println("\n2. Integrating S[1,1]*S[2,2] (Should be 0.5 for d=2)")
expr_sp1 = S[1, 1]*S[2, 2]
res_sp1 = integrate(expr_sp1, mS)
println("Result: ", res_sp1)

println("\n3. Integrating S[1,2]*S[2,1] (Should be -0.5)")
res_sp2 = integrate(S[1, 2]*S[2, 1], mS)
println("Result: ", res_sp2)

println("\n4. Integrating |S[1,1]|^2 (Should be 0.5 for d=2)")
expr_sp3 = abs(S[1, 1])^2
res_sp3 = integrate(expr_sp3, mS)
println("Result: ", res_sp3)
