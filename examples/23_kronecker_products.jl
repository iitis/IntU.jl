using IntU
using Symbolics
using LinearAlgebra

# Example 23: Kronecker Product Integration
# This example demonstrates how to integrate expressions involving 
# Kronecker products of symbolic matrices.

# 1. Integration with Concrete Dimensions
println("--- 1. Concrete Dimensions ---")
U = symbolic_unitary(:U, 3)
# Kronecker product remains lazy (SymbolicKron)
K = kron(U, U)
K_adj = adjoint(K)

# tr( (U ⊗ U) * (U† ⊗ U†) ) = tr(U*U†) * tr(U*U†) = 3 * 3 = 9
expr1 = tr(K * K_adj)
res1 = integrate(expr1, dU(3))
println("tr( (U ⊗ U) * (U† ⊗ U†) ) integrated over dU(3): ", res1)

# 2. Integration with Symbolic Dimensions
println("\n--- 2. Symbolic Dimensions ---")
@variables d
Ud = symbolic_unitary(:Ud, d)
Kd = kron(Ud, Ud)
Kd_adj = adjoint(Kd)

# tr( (Ud ⊗ Ud) * (Ud† ⊗ Ud†) ) = d^2
expr2 = tr(Kd * Kd_adj)
res2 = integrate(expr2, dU(d))
println("tr( (Ud ⊗ Ud) * (Ud† ⊗ Ud†) ) integrated over dU(d): ", simplify(res2))

# 3. Matrix-Valued Integration with @integrate
println("\n--- 3. Matrix-Valued Integration ---")
B = SymbolicMatrix(:B, :Constant, 9)

# @integrate handles :kron automatically
res3 = @integrate kron(U, U) * B * kron(U', U') dU(3)

println("Matrix-valued integral size: ", size(res3))
println("First element of result: ", res3[1, 1])

# But here we just show it computes successfully.
println("Result 3[1,1]: ", res3[1, 1])

# 4. Dimension-less Neighbors (Dimension Propagation)
println("\n--- 4. Dimension Propagation ---")
# B has no dimension specified
B_un = symbolic_unitary(:B_un, nothing)
# kron(U, U) is 4x4. B_un will correctly infer 4x4 size from its neighbors.
expr4 = kron(U, U) * B_un * adjoint(kron(U, U))
res4 = integrate(expr4, dU(2))

println("Integral with dimension-less B succeeded.")
println("Size of result 4: ", size(res4))

println("\nExample 23 completed successfully.")
