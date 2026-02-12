using IntU
using Symbolics
using LinearAlgebra

println("=== Gaussian Ensembles Integration ===\n")

# --- 1. Explicit Matrix (Small N) ---
N = 2
println("1. Explicit Matrix (N=$N)")
# Use SymbolicMatrix to ensure correct metadata for integration engine
H_sym = SymbolicMatrix(:H_mat, :GUE, N)
# Extract explicit matrix of variables (which now have metadata)
H_mat = [H_sym[i,j] for i=1:N, j=1:N]

println("--- GUE ---")
res_gue = @integrate tr(H_mat^2) dGUE(H_sym, N)
println("<Tr(H^2)>_GUE = ", Symbolics.simplify(res_gue), " (Expected: $(N^2))")

println("--- GOE ---")
res_goe = @integrate tr(H_mat^2) dGOE(H_sym, N)
println("<Tr(H^2)>_GOE = ", Symbolics.simplify(res_goe), " (Expected: $(N^2 + N))")

println("--- GSE ---")
res_gse = @integrate tr(H_mat^2) dGSE(H_mat, N)
println("<Tr(H^2)>_GSE = ", Symbolics.simplify(res_gse), " (Expected: $(N^2 - N))")


println("\n--- Matrix Averages over GUE(2) ---")
res_sq = @integrate H_mat^2 dGUE(H_mat, N)
# Simplify result for display
res_sq_simp = map(x -> simplify(x), res_sq)
println("<H^2>_GUE (Should be Diagonal matrix N*I):")
display(res_sq_simp)

# --- 2. Symbolic Dimension and Traces ---
println("\n2. Symbolic Dimension and Traces")
@variables d
H = SymbolicMatrix(:H)

println("--- GUE ---")
# <Tr(H^4)> = 2d^3 + d
res_gue4 = @integrate tr(H^4) dGUE(H, d)
println("<Tr(H^4)>_GUE = ", simplify(res_gue4))

println("--- GOE ---")
# <Tr(H^4)> = 2d^3 + 5d^2 + 5d
res_goe4 = @integrate tr(H^4) dGOE(H, d)
println("<Tr(H^4)>_GOE = ", simplify(res_goe4))

println("--- GSE ---")
# <Tr(H^4)> = 2d^3 - 5d^2 + 5d
res_gse4 = @integrate tr(H^4) dGSE(H, d)
println("<Tr(H^4)>_GSE = ", simplify(res_gse4))

println("\nDone.")
