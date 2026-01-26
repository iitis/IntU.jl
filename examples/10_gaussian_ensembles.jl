using IntU
using Symbolics

println("=== Gaussian Ensembles Integration ===\n")

# --- 1. Explicit Matrix (Small N) ---
N = 2
println("1. Explicit Matrix (N=$N)")
H_explicit = [Symbolics.variable(:H, i, j) for i in 1:N, j in 1:N]

println("--- GUE ---")
res_gue2 = integrate(IntU.tr(H_explicit^2), dGUE(H_explicit, N))
println("<Tr(H^2)>_GUE = ", res_gue2, " (Expected: $(N^2))")

println("--- GOE ---")
res_goe2 = integrate(IntU.tr(H_explicit^2), dGOE(H_explicit, N))
println("<Tr(H^2)>_GOE = ", res_goe2, " (Expected: $(N^2 + N))")


# --- 2. Symbolic Dimension and Coordinate-free Trace ---
println("\n2. Symbolic Dimension and Traces")
@variables d
H = SymbolicMatrix(:H)

println("--- GUE ---")
# <Tr(H^4)> = 2d^3 + d
res_gue4 = integrate(IntU.tr(H^4), dGUE(H, d))
println("<Tr(H^4)>_GUE = ", res_gue4)

println("--- GOE ---")
# <Tr(H^4)> = 2d^3 + 5d^2 + 5d
res_goe4 = integrate(IntU.tr(H^4), dGOE(H, d))
println("<Tr(H^4)>_GOE = ", res_goe4)

println("\nDone.")
