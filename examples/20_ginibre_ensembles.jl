using IntU
using Symbolics
using LinearAlgebra
import LinearAlgebra: tr

println("=== Ginibre Ensembles Integration Examples ===\n")

# --- 1. Complex Ginibre (GinUE) ---
println("1. Complex Ginibre Ensemble (GinUE)")
N = 2
# Specification of T=Complex{Num} is important for GinUE to ensure conj(G) != G
G = [Symbolics.variable(:G, i, j, T = Complex{Num}) for i = 1:N, j = 1:N]
meas_ginue = dGinUE(G, N)

println("--- Basic Moments ---")
# <Tr(G G')> = N^2
res_sq = simplify(integrate(tr(G * G'), meas_ginue))
println("<Tr(G G')> = ", res_sq, " (Expected: $(N^2))")

# <Tr(G G' G G')> = N^4 + N^2 (for GinUE)
# Wait, for GinUE <G_ij conj(G_kl)> = delta_ik delta_jl. 
# Tr(G G' G G') = G_ij conj(G_kj) G_kl conj(G_il)
# Contractions:
# 1. (ij, kj) and (kl, il) -> delta_ik delta_jj * delta_ki delta_ll = N * N = N^2? No.
# Actually for GUE <Tr(H^4)> = 2N^3 + N. GinUE is different as it's not Hermitian.
res_fourth = simplify(integrate(tr(G * G' * G * G'), meas_ginue))
println("<Tr(G G' G G')> = ", res_fourth)

# --- 2. Real Ginibre (GinOE) ---
println("\n2. Real Ginibre Ensemble (GinOE)")
# entries are real, so G_ij is Hermitian to itself
G_real = [Symbolics.variable(:G_r, i, j) for i = 1:N, j = 1:N]
meas_ginoe = dGinOE(G_real, N)

# <Tr(G G^T)> = N^2
res_oe_sq = simplify(integrate(tr(G_real * transpose(G_real)), meas_ginoe))
println("<Tr(G G^T)> = ", res_oe_sq, " (Expected: $(N^2))")

# --- 3. Graphical Calculus (Symbolic Dimension) ---
println("\n3. Graphical Calculus and Symbolic Dimension")
@variables d
Gs = SymbolicMatrix(:G)
meas_s = dGinUE(Gs, d)

# <Tr(G G' G G')> with symbolic d
t = tr_lazy(Gs * Gs' * Gs * Gs')
res_s = simplify(integrate(t, meas_s))
println("<Tr(G G' G G')>_GinUE = ", res_s)

# Coordinate-free matrix integration using SymbolicMatrix for A and B
println("\n--- Coordinate-free Symbolic Matrix Integration ---")
As = SymbolicMatrix(:A)
Bs = SymbolicMatrix(:B)
# <Tr(G A G' B)> = Tr(A) * Tr(B)
t2 = tr_lazy(Gs * As * Gs' * Bs)
res_t2 = simplify(integrate(t2, meas_s))
println("<Tr(G A G' B)> = ", res_t2, " (Should be tr(A)*tr(B))")

# Matrix Integration
println("\n4. Matrix Integrals (Numerical A, B)")
A = [1 0; 0 2]
B = [1 1; 1 1]
# < G A G' B > = Tr(A) B
res_mat = integrate(G * A * G' * B, meas_ginue)
println("< G A G' B > = ")
display(map(simplify, res_mat))

# Symbolic A and B using Array Variables (avoids separate entries)
println("\n5. Matrix Integrals with Symbolic Array Variables")
# Use Matrix{Num} to ensure standard indexing works
As = [Symbolics.variable(:A, i, j) for i = 1:N, j = 1:N]
Bs = [Symbolics.variable(:B, i, j) for i = 1:N, j = 1:N]
res_sym = integrate(G * As * G' * Bs, meas_ginue)
println("< G As G' Bs > = ")
display(map(simplify, res_sym))

# Verification
expected_sym = tr(As) * Bs
is_correct = all(
    i -> IntU._symbolic_isequal(simplify(res_sym[i]), simplify(expected_sym[i])),
    eachindex(res_sym),
)
println("\nMatches Tr(A)*B: ", is_correct)

println("\nDone.")
