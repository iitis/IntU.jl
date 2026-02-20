using IntU
using Symbolics
using LinearAlgebra

println("=== Circular Ensembles Integration ===\n")

@variables d

# --- 1. Circular Orthogonal Ensemble (COE) ---
println("--- 1. COE (Circular Orthogonal Ensemble) ---")
println("Matrix O is symmetric unitary: O = O^T.")
# The macro expects variable O for dCOE

# Moment: E[|O_11|^2]
# For COE, E[|O_ij|^2] = (1 + delta_ij) / (d + 1)
# So E[|O_11|^2] = 2/(d+1)
println("Integrating |O[1,1]|^2 over COE(d)...")
res_coe = @integrate abs(O[1, 1])^2 dCOE(d)
println("E[|O_11|^2] = $(Symbolics.simplify(res_coe)) (Expected: 2/(d+1))")

# Moment: E[|O_12|^2] = 1/(d+1)
println("Integrating |O[1,2]|^2 over COE(d)...")
res_coe_12 = @integrate abs(O[1, 2])^2 dCOE(d)
println("E[|O_12|^2] = $(Symbolics.simplify(res_coe_12)) (Expected: 1/(d+1))")

println("\n--- Matrix Integration (COE) over COE(2) ---")
println("Integrating O * O^T (should be Identity)")
res_S = @integrate O * O' dCOE(2)
display(res_S)


# --- 2. Circular Symplectic Ensemble (CSE) ---
println("\n--- 2. CSE (Circular Symplectic Ensemble) ---")
println("Matrix Sp is self-dual unitary: Sp = J Sp^T J^T.")
println("Note: CSE dimension must be even.")

# Moment: E[|Sp_11|^2]
# For CSE, E[|Sp_ii|^2] = 1/(d-1)
println("Integrating |Sp[1,1]|^2 over CSE(4)...")
res_cse = @integrate abs(Sp[1, 1])^2 dCSE(4)
println("E[|Sp_11|^2] = $(res_cse) (Expected: 1/(4-1) = 1/3)")


# --- 3. Circular Unitary Ensemble (CUE) ---
println("\n--- 3. CUE (Circular Unitary Ensemble) ---")
println("CUE is equivalent to Haar Unitary group.")
U = SymbolicMatrix(:U, :U)

# Moment: E[|U_11|^2] = 1/d
println("Integrating |U[1,1]|^2 over CUE(d)...")
res_cue = integrate(abs(U[1, 1])^2, dCUE(d))
println("E[|U_11|^2] = $(Symbolics.simplify(res_cue)) (Expected: 1/d)")

println("\nDone.")
