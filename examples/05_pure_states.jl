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
rho1(a, c) = sum(U[(a-1)*d2 + k, 1] * conj(U[(c-1)*d2 + k, 1]) for k=1:d2)

# Subsystem Purity Tr(rho1^2)
purity_expr = sum(rho1(a, c) * rho1(c, a) for a=1:d1, c=1:d1)

# Integrate over Haar measure dU(d) which induces the measure on pure states
# Integrate over Haar measure dU(d) which induces the measure on pure states
println("Integrating purity_expr...")
avg_purity = @integrate purity_expr dU(d)
println("Integration done.")

println("Average Purity: ", avg_purity)
println("Expected (Analytical): (d1 + d2) / (d1 * d2 + 1)")

# Average Fidelity between a random state and a fixed state |0>
# <psi|0> = U[1, 1]
fidelity_expr = abs(U[1, 1])^2
fidelity_expr = abs(U[1, 1])^2
println("Integrating fidelity_expr...")
avg_fidelity = @integrate fidelity_expr dU(d)
println("Integration done.")
println("\nAverage Fidelity with |0>: ", avg_fidelity)
println("Expected (1/d): 1/d")

# Convenience dPsi shorthand:
println("\ndPsi(d) can also be used for purely vector-based integration.")
println("\ndPsi(d) can also be used for purely vector-based integration.")
println("Integrating with dPsi...")
res_psi = @integrate abs(U[1,1])^2 dPsi(U[:, 1], d)
println("Integration done.")
println("Result using dPsi(d): ", res_psi)
