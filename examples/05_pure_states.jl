using IntU
using Symbolics
using LinearAlgebra

# Example: Average Purity of a Haar-random pure state
# For a pure state psi, rho = psi * psi'
# We can represent it as the first column of a Haar-random unitary U.

d1 = 2
d2 = 2
d = d1 * d2
U = SymbolicMatrix(:U, :U, d)

# Subsystem 1 density matrix from the first column of U
# rho1_ij = sum_k U_{ik} * conj(U_{jk}) where indices are mapped appropriately
# For simplicity in this symbolic example, we use the entry-wise definition:
println("Calculating average purity for a bipartite state of dimension d1*d2...")

# Pure state |psi> = U[:, 1]. We map the single index i to (a, b)
# psi_{a,b} = U[(a-1)*d2 + b, 1]
rho1(a, c) = sum(U[(a-1)*d2+k, 1] * conj(U[(c-1)*d2+k, 1]) for k = 1:d2)

# Subsystem Bipartite Trace Tr(rho1^2)
tr_rho1_2 = sum(rho1(a, c) * rho1(c, a) for a = 1:d1, c = 1:d1)

println("Integrating tr_rho1_2...")
avg_tr_rho1_2 = integrate(tr_rho1_2, dU(d))
println("Integration done.")

println("Average Bipartite Trace: ", avg_tr_rho1_2)
println("Expected (Analytical): (d1 + d2) / (d1 * d2 + 1)")

# Average Overlap between a random state and a fixed state |0>
# <psi|0> = U[1, 1]
overlap_expr = abs(U[1, 1])^2
println("Integrating overlap_expr...")
avg_overlap = integrate(overlap_expr, dU(d))
println("Integration done.")
println("\nAverage Overlap with |0>: ", avg_overlap)
println("Expected (1/d): 1/d")

# Convenience dPsi shorthand:
println("\ndPsi(d) can also be used for purely vector-based integration.")
println("Integrating with dPsi...")
psi = SymbolicMatrix(:psi, :psi, (d, 1))
res_psi = integrate(abs(psi[1, 1])^2, dPsi(d))
println("Integration done.")
println("Result using dPsi(d): ", res_psi)
