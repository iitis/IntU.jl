using IntU
using Symbolics

println("=== Gaussian Ensembles Integration ===\n")

# --- 1. Explicit Matrix (Small N) ---
N = 2
println("1. Explicit Matrix (N=$N)")
H_explicit = [Symbolics.variable(:H, i, j) for i in 1:N, j in 1:N]

println("--- GUE ---")
res_gue2 = simplify(integrate(IntU.tr(H_explicit^2), dGUE(H_explicit, N)))
println("<Tr(H^2)>_GUE = ", res_gue2, " (Expected: $(N^2))")

println("--- GOE ---")
res_goe2 = simplify(integrate(IntU.tr(H_explicit^2), dGOE(H_explicit, N)))
println("<Tr(H^2)>_GOE = ", res_goe2, " (Expected: $(N^2 + N))")

println("--- GSE ---")
res_gse2 = simplify(integrate(IntU.tr(H_explicit^2), dGSE(H_explicit, N)))
println("<Tr(H^2)>_GSE = ", res_gse2, " (Expected: $(N^2 - N))")


# --- 2. Symbolic Dimension and Traces ---
println("\n2. Symbolic Dimension and Traces")
@variables d
H = SymbolicMatrix(:H)

println("--- GUE ---")
# <Tr(H^4)> = 2d^3 + d
res_gue4 = simplify(integrate(IntU.tr(H^4), dGUE(H, d)))
println("<Tr(H^4)>_GUE = ", res_gue4)

println("--- GOE ---")
# <Tr(H^4)> = 2d^3 + 5d^2 + 5d
res_goe4 = simplify(integrate(IntU.tr(H^4), dGOE(H, d)))
println("<Tr(H^4)>_GOE = ", res_goe4)

println("--- GSE ---")
# <Tr(H^4)> = 2d^3 - 5d^2 + 5d
res_gse4 = simplify(integrate(IntU.tr(H^4), dGSE(H, d)))
println("<Tr(H^4)>_GSE = ", res_gse4)

println("\nDone.")
