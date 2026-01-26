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
expr1 = O[1,1]^2
res1 = integrate(expr1, mO)
println("Result: ", res1)

println("\n2. Integrating |O[1,1]|^4 (should be 3/(d(d+2)))")
expr2 = O[1,1]^4
res2 = integrate(expr2, mO)
println("Result: ", Symbolics.simplify(res2))

println("\n3. Orthogonality check: sum_k O[1,k]*O[1,k] (should be 1)")
expr3 = sum(O[1,k]*O[1,k] for k in 1:3)
# Note: integrating over O(d) usually implies d matches the matrix size indices if we sum over all of them?
# If we treat d as symbolic but sum k=1..3, we get 3/d.
# If we set d=3, we get 1.
res3 = integrate(expr3, mO)
println("Result (symbolic sum k=1..3): ", res3)
println("Value at d=3: ", Symbolics.substitute(res3, Dict(d => 3)))


# --- Symplectic Group Sp(d) ---
println("\n--- Symplectic Group Sp(d) ---")
println("(Note: Sp(d) integration requires d to be even)")

@variables S[1:2, 1:2]
mS = dSp(S, d)

println("1. Defining Measure dSp(S, d)")
println("Measure: ", mS)

try
    println("\n2. Attempting integration of |S[1,1]|^2 (Experimental)")
    expr_sp = S[1,1]^2
    res_sp = integrate(expr_sp, mS)
    println("Result: ", res_sp)
catch e
    println("Symplectic integration is currently a work in progress and throws an error as expected.")
    println("Error: ", e)
end
