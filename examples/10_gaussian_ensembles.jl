using IntU
using Symbolics

println("=== Gaussian Ensembles Integration ===")

N = 3
println("Using Matrix Size N = $N\n")

# Define symbolic matrix entries
H = [Symbolics.variable(:H, i, j) for i in 1:N, j in 1:N]

# --- GUE Section ---
println("--- Gaussian Unitary Ensemble (GUE) ---")
meas_gue = dGUE(H, N)

# < Tr(H^2) > = N^2
res_gue2 = integrate(IntU.tr(H^2), meas_gue)
println("< Tr(H^2) >_GUE = ", res_gue2, " (Expected: $(N^2))")

# < Tr(H^4) > = 2N^3 + N
res_gue4 = integrate(IntU.tr(H^4), meas_gue)
println("< Tr(H^4) >_GUE = ", res_gue4, " (Expected: $(2*N^3 + N))")

# --- GOE Section ---
println("\n--- Gaussian Orthogonal Ensemble (GOE) ---")
meas_goe = dGOE(H, N)

# < Tr(H^2) > = N^2 + N
res_goe2 = integrate(IntU.tr(H^2), meas_goe)
println("< Tr(H^2) >_GOE = ", res_goe2, " (Expected: $(N^2 + N))")

# < Tr(H^4) > = 2N^3 + 5N^2 + 5N
res_goe4 = integrate(IntU.tr(H^4), meas_goe)
println("< Tr(H^4) >_GOE = ", res_goe4, " (Expected: $(2*N^3 + 5*N^2 + 5*N))")

println("\nDone.")
