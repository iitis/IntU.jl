using IntU
using Symbolics
using LinearAlgebra

println("=== Circular Ensembles Integration ===\n")

@variables d

# --- 1. Circular Orthogonal Ensemble (COE) ---
println("--- 1. COE (Circular Orthogonal Ensemble) ---")
println("Matrix S is symmetric unitary: S = S^T.")
# Define a 2x2 symbolic matrix for S
N = 2
@variables S[1:N, 1:N]::Complex
# COE Measure
m_coe = dCOE(S, d)

# Moment: E[S_11 * S*_11]
# For COE, E[|S_ij|^2] = (1 + delta_ij) / (d + 1)
# So E[|S_11|^2] = 2/(d+1)
expr_coe = S[1, 1] * conj(S[1, 1])
res_coe = simplify(integrate(expr_coe, m_coe); expand = true)
println("E[|S_11|^2] = $(res_coe) (Expected: 2/(d+1))")

# Moment: E[|S_12|^2] = 1/(d+1)
expr_coe_12 = S[1, 2] * conj(S[1, 2])
res_coe_12 = simplify(integrate(expr_coe_12, m_coe); expand = true)
println("E[|S_12|^2] = $(res_coe_12) (Expected: 1/(d+1))")


# --- 2. Circular Symplectic Ensemble (CSE) ---
println("\n--- 2. CSE (Circular Symplectic Ensemble) ---")
println("Matrix S is self-dual unitary: S = J S^T J^T.")
println("Note: CSE dimension must be even. Here we use 2x2 (N=1).")
# Define a 2x2 symbolic matrix for S_cse
@variables S_cse[1:N, 1:N]::Complex
m_cse = dCSE(S_cse, d) # d corresponds to the full dimension 2N

# Moment: E[|S_11|^2]
# For CSE, E[|S_ii|^2] = 1/(d-1)
expr_cse = S_cse[1, 1] * conj(S_cse[1, 1])
res_cse = simplify(integrate(expr_cse, m_cse); expand = true)
println("E[|S_11|^2] = $(res_cse) (Expected: 1/(d-1))")

# Moment: E[|S_12|^2]
# For N=1 (d=2), |S_12|^2 is technically correlated with |S_11|^2 by unitarity.
# E[|S_{12}|^2] = (d-2)/(d-1)^2? For d=2, this is 0.
expr_cse_12 = S_cse[1, 2] * conj(S_cse[1, 2])
res_cse_12 = simplify(integrate(expr_cse_12, m_cse); expand = true)
println("E[|S_12|^2] = $(res_cse_12) (Expected: (d-2)/(d-1)^2 -> 0 for d=2)")

# --- 3. Circular Unitary Ensemble (CUE) ---
println("\n--- 3. CUE (Circular Unitary Ensemble) ---")
println("CUE is equivalent to Haar LoE.")
@variables U[1:N, 1:N]::Complex
m_cue = dCUE(U, d)

# Moment: E[|U_11|^2] = 1/d
expr_cue = U[1, 1] * conj(U[1, 1])
res_cue = simplify(integrate(expr_cue, m_cue); expand = true)
println("E[|U_11|^2] = $(res_cue) (Expected: 1/d)")

# Moment: E[|U_11 U_22|^2] = 1/(d^2 - 1)
expr_cue_4 = abs(U[1, 1] * U[2, 2])^2
res_cue_4 = simplify(integrate(expr_cue_4, m_cue); expand = true)
println("E[|U_11 U_22|^2] = $(res_cue_4) (Expected: 1/(d^2 - 1))")

println("\nDone.")
