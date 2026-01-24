using IntU
using Symbolics
using LinearAlgebra

# Example: Average Purity of a Haar-random pure state
# For a pure state psi, rho = psi * psi'
# Purity Tr(rho^2) = Tr((psi * psi')^2) = (psi' * psi)^2 = 1 (always)
# But let's look at a subsystem purity.

d1 = 2
d2 = 2
d = d1 * d2
@variables psi_re[1:d] psi_im[1:d]
psi = [(psi_re[i] + im*psi_im[i]) for i in 1:d]
measure = dPsi(psi, d)

# Average Purity of subsystem 1: E[Tr(rho1^2)]
psi_mat = Array{Any}(undef, d1, d2)
let idx = 1
    for i in 1:d1
        for j in 1:d2
            psi_mat[i,j] = psi[idx]
            idx += 1
        end
    end
end

rho1 = [sum(psi_mat[i,k] * conj(psi_mat[j,k]) for k in 1:d2) for i in 1:d1, j in 1:d1]
purity_expr = sum(rho1[i,j] * conj(rho1[i,j]) for i in 1:d1, j in 1:d1)

println("Calculating average purity for d1=$d1, d2=$d2...")
avg_purity = integrate(purity_expr, measure)
println("Average Purity: ", avg_purity)
println("Expected (Analytical): ", (d1 + d2) / (d1 * d2 + 1))

# Average Fidelity between a random state and a fixed state |0>
fidelity_expr = conj(psi[1]) * psi[1]
avg_fidelity = integrate(fidelity_expr, measure)
println("Average Fidelity with |0>: ", avg_fidelity)
println("Expected (1/d): ", 1//d)
