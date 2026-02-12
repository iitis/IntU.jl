# examples/06_quantum_information.jl
using IntU
using Symbolics
using LinearAlgebra

println("Quantum Information Example: Average Purity of a Subsystem")

# 1. Setup a bipartite system (d = d1 * d2)
d1 = 2
d2 = 2
d = d1 * d2
U = SymbolicMatrix(:U, :U)

# 2. Start with a product state |00>
# We represent the random state as the first column of U
# psi = U[:, 1]
# For a density matrix rho = |psi><psi|

# 3. Apply the random unitary (inherent in U definition)
# 4. Compute the reduced density matrix of subsystem 1 (trace out system 2)
# Since we are doing it symbolically, we define the RDM entry-wise:
# (rho1)_{ac} = sum_k U_{(a,k), 1} * conj(U_{(c,k), 1})
rho1(a, c) = sum(U[(a-1)*d2 + k, 1] * conj(U[(c-1)*d2 + k, 1]) for k=1:d2)

# 5. Measure the purity of the reduced density matrix: P = Tr(rho1^2)
# Here we just show the setup and use the automated integration
purity_expr = sum(rho1(a, c) * rho1(c, a) for a=1:d1, c=1:d1)

println("Integrating purity over the Haar measure dU(d)...")
avg_purity = @integrate purity_expr dU(d)

println("\nResults for symbolic d = d1*d2:")
println("Average Purity of Subsystem 1: ", Symbolics.simplify(avg_purity))

# Theoretical expectation for U(d1 * d2):
# <Tr(rho1^2)> = (d1 + d2) / (d1 * d2 + 1)
println("Theoretical expectation: (d1 + d2) / (d1 * d2 + 1)")

# Check for specific dimensions
