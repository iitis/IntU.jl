using IntU
using Symbolics
using LinearAlgebra

println("=== Circular Ensembles Integration ===\n")

@variables d

# --- 1. Circular Orthogonal Ensemble (COE) ---
println("--- 1. COE (Circular Orthogonal Ensemble) ---")
println("Matrix S is symmetric unitary: S = S^T.")
# Define a SymbolicMatrix for S
S_coe = SymbolicMatrix(:S, :COE)

# Moment: E[|S_11|^2]
# For COE, E[|S_ij|^2] = (1 + delta_ij) / (d + 1)
# So E[|S_11|^2] = 2/(d+1)
println("Integrating |S[1,1]|^2 over COE(d)...")
res_coe = @integrate abs(S_coe[1, 1])^2 dCOE(d)
println("E[|S_11|^2] = $(Symbolics.simplify(res_coe)) (Expected: 2/(d+1))")

# Moment: E[|S_12|^2] = 1/(d+1)
println("Integrating |S[1,2]|^2 over COE(d)...")
res_coe_12 = @integrate abs(S_coe[1, 2])^2 dCOE(d)
println("E[|S_12|^2] = $(Symbolics.simplify(res_coe_12)) (Expected: 1/(d+1))")

println("\n--- Matrix Integration (COE) over COE(2) ---")
println("Integrating S * S'\' (should be Identity)")
res_S = @integrate S_coe * S_coe' dCOE(2)
display(res_S)


# --- 2. Circular Symplectic Ensemble (CSE) ---
println("\n--- 2. CSE (Circular Symplectic Ensemble) ---")
println("Matrix S is self-dual unitary: S = J S^T J^T.")
println("Note: CSE dimension must be even.")

# Moment: E[|S_11|^2]
# For CSE, E[|S_ii|^2] = 1/(d-1)
S_cse = SymbolicMatrix(:S, :CSE)
println("Integrating |S[1,1]|^2 over CSE(4)...")
res_cse = @integrate abs(S_cse[1, 1])^2 dCSE(4)
println("E[|S_11|^2] = $(res_cse) (Expected: 1/(4-1) = 1/3)")


# --- 3. Circular Unitary Ensemble (CUE) ---
println("\n--- 3. CUE (Circular Unitary Ensemble) ---")
println("CUE is equivalent to Haar Unitary group.")
U = SymbolicMatrix(:U, :U)

# Moment: E[|U_11|^2] = 1/d
println("Integrating |U[1,1]|^2 over CUE(d)...")
res_cue = @integrate abs(U[1, 1])^2 dCUE(d)
println("E[|U_11|^2] = $(Symbolics.simplify(res_cue)) (Expected: 1/d)")

println("\nDone.")
