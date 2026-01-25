# examples/06_quantum_information.jl
using IntU
using Symbolics
using LinearAlgebra

println("Quantum Information Example: Average Purity of a Subsystem")

# 1. Setup a bipartite system (2 qubits, d=4)
# System A: qubit 1, System B: qubit 2
dims = (2, 2)
d = prod(dims)

@variables U[1:d, 1:d]::Complex
measure = dU(U, d)

# 2. Start with a product state |00>
# rho0 = |00><00|
rho0 = zeros(d, d)
rho0[1, 1] = 1.0

# 3. Apply a random unitary U
# rho = U * rho0 * U'
# Use collect to handle matrix multiplication with Symbolics.Arr
U_m = collect(U)
rho = U_m * rho0 * U_m'

# 4. Compute the reduced density matrix of subsystem A (trace out B)
println("Tracing out subsystem B (index 2)...")
rho_A = partial_trace(rho, dims, 2)

# 5. Measure the purity of the reduced density matrix
# purity = Tr(rho_A^2)
pur = purity(rho_A)
println("Integrating purity over the Haar measure...")
# Better: integrate the purity of the SUBSYSTEM
avg_pur_sub = integrate(pur, measure)

println("\nResults for d=$d:")
println("Average Purity of Subsystem A: ", Symbolics.simplify(avg_pur_sub))

# Theoretical expectation for U(d_A * d_B):
# <Tr(rho_A^2)> = (d_A + d_B) / (d_A * d_B + 1)
# For d_A=2, d_B=2: (2+2)/(4+1) = 4/5 = 0.8
expected = (dims[1] + dims[2]) / (prod(dims) + 1)
println("Theoretical expectation: ", expected)

# Results already match theoretical expectation of 0.8
