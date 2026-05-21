using IntegrateUnitary
using Symbolics
using LinearAlgebra

println("=== Ginibre Ensembles Integration Examples ===\n")

# --- 1. Complex Ginibre (GinUE) ---
println("1. Complex Ginibre Ensemble (GinUE)")
N = 2
G_ue = SymbolicMatrix(:G, :GinUE, N)
G_mat = [G_ue[i, j] for i = 1:N, j = 1:N]
meas_ginue = dGinUE(N)

println("--- Basic Moments ---")
# <Tr(G G')> = N^2
res_sq = integrate(tr(G_mat * G_mat'), meas_ginue)
println("<Tr(G G')> = ", Symbolics.simplify(res_sq), " (Expected: $(N^2))")

# --- 2. Real Ginibre (GinOE) ---
println("\n2. Real Ginibre Ensemble (GinOE)")
G_oe = SymbolicMatrix(:G, :GinOE, N)
Gr = [G_oe[i, j] for i = 1:N, j = 1:N]
meas_ginoe = dGinOE(N)

# <Tr(G G^T)> = N^2
res_oe_sq = integrate(tr(Gr * Gr'), meas_ginoe)
println("<Tr(G G^T)> = ", Symbolics.simplify(res_oe_sq), " (Expected: $(N^2))")


# --- 3. Graphical Calculus (Symbolic Dimension) ---
println("\n3. Graphical Calculus and Symbolic Dimension")
@variables d
G_sym = SymbolicMatrix(:G, :GinUE)
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)

# <Tr(G A G' B)> = Tr(A) * Tr(B)
println("Integrating tr(G * A * G' * B) over GinUE(d)...")
res_t2 = integrate(tr(G_sym * A * G_sym' * B), dGinUE(d))
println("<Tr(G A G' B)> = ", res_t2, " (Expected: tr(A)*tr(B))")


# --- 4. Matrix Integrals ---
println("\n4. Matrix Integrals over GinUE(2)")
A_const = [1 0; 0 2]
B_const = [1 1; 1 1]
# < G A G' B > = Tr(A) B
res_mat = integrate((G_mat * A_const) * (G_mat' * B_const), meas_ginue)
println("< G A G' B > = ")
display(map(simplify, res_mat))

println("\nDone.")
