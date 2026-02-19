using IntU
using Symbolics
using LinearAlgebra

println("=== Gaussian Ensembles Integration ===\n")

# --- 1. Explicit Matrix (Small N) ---
N = 2
println("1. Explicit Matrix (N=$N)")
# Use SymbolicMatrix to ensure correct metadata for integration engine
# GUE
H_gue = SymbolicMatrix(:H_gue, :GUE, N)
H_mat_gue = [H_gue[i,j] for i=1:N, j=1:N]

println("--- GUE ---")
res_gue = integrate(tr(H_mat_gue^2), dGUE(N))
println("<Tr(H^2)>_GUE = ", Symbolics.simplify(res_gue), " (Expected: $(N^2))")

println("--- GOE ---")
# GOE
H_goe_e = SymbolicMatrix(:H_goe, :GOE, N)
H_mat_goe = [H_goe_e[i,j] for i=1:N, j=1:N]
res_goe = integrate(tr(H_mat_goe^2), dGOE(N))
println("<Tr(H^2)>_GOE = ", Symbolics.simplify(res_goe), " (Expected: $(N^2 + N))")

println("--- GSE ---")
# GSE
H_gse_e = SymbolicMatrix(:H_gse, :GSE, N)
H_mat_gse = [H_gse_e[i,j] for i=1:N, j=1:N]
res_gse = integrate(tr(H_mat_gse^2), dGSE(N))
println("<Tr(H^2)>_GSE = ", Symbolics.simplify(res_gse), " (Expected: $(N^2 - N))")


println("\n--- Matrix Averages over GUE(2) ---")
res_sq = integrate(H_mat_gue^2, dGUE(N))
# Simplify result for display
res_sq_simp = map(x -> simplify(x), res_sq)
println("<H^2>_GUE (Should be Diagonal matrix N*I):")
display(res_sq_simp)

# --- 2. Symbolic Dimension and Traces ---
println("\n2. Symbolic Dimension and Traces")
@variables d
H = SymbolicMatrix(:H, :GUE) # Use standard tag :GUE for GUE

println("--- GUE ---")
# <Tr(H^4)> = 2d^3 + d
res_gue4 = integrate(tr(H^4), dGUE(d))
println("<Tr(H^4)>_GUE = ", simplify(res_gue4))

println("--- GOE ---")
# <Tr(H^4)> = 2d^3 + 5d^2 + 5d
H_goe = SymbolicMatrix(:H, :GOE)
res_goe4 = integrate(tr(H_goe^4), dGOE(d))
println("<Tr(H^4)>_GOE = ", simplify(res_goe4))

println("--- GSE ---")
# <Tr(H^4)> = 2d^3 - 5d^2 + 5d
H_gse = SymbolicMatrix(:H, :GSE)
res_gse4 = integrate(tr(H_gse^4), dGSE(d))
println("<Tr(H^4)>_GSE = ", simplify(res_gse4))

println("\nDone.")
