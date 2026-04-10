using IntU
using Symbolics
using LinearAlgebra

# Example 1: Basic scalar integrals
println("--- Example 1: Scalar Integrals (Symbolic d) ---")
@variables d
U = SymbolicMatrix(:U, :U)

println("1.1: Integrate |U[1,1]|^2")
res1_1 = integrate(abs(U[1, 1])^2, dU(d))
println("Result: ", res1_1, " (Expected: 1/d)")

println("\n1.2: Integrate |U[1,1]*U[2,2]|^2")
res1_2 = integrate(abs(U[1, 1] * U[2, 2])^2, dU(d))
println("Result: ", res1_2, " (Expected: 1/(d^2-1))")

println("\n1.3: Integrate U[1,1]*U[2,2]*conj(U[1,2]*U[2,1])")
res1_3 = integrate(U[1, 1] * U[2, 2] * conj(U[1, 2] * U[2, 1]), dU(d))
println("Result: ", res1_3, " (Expected: -1/(d(d^2-1)))")


# Example 2: Index-based integration (Low-level functionality)
println("\n--- Example 2: Index-based integration ---")
I1 = [1, 1, 1, 2, 2];
J1 = [2, 2, 1, 1, 1]
I2 = [1, 1, 1, 2, 2];
J2 = [2, 1, 1, 2, 1]
d2 = 6

u_idxs = collect(zip(I1, J1))
u_bar_idxs = collect(zip(I2, J2))
res2 = integrate_indices(u_idxs, u_bar_idxs, d2)
println("Result: ", res2, " (Expected: -1/16200)")


# Example 3: Multiple unitaries
println("\n--- Example 3: Multiple Unitaries ---")
@variables dU_dim, dV_dim
U = SymbolicMatrix(:U, :U, dU_dim)
V = SymbolicMatrix(:V, :V, dV_dim)
X = SymbolicMatrix(:X) # Constant matrix

# Complex expression with two independent unitaries
# ∫ dU ∫ dV (U ⊗ V) X (U ⊗ V)†

# Define elements for a small block demonstration
integrand3 = (U[1, 1] * V[1, 1]) * X[1, 1] * conj(U[1, 1] * V[1, 1])

println("Integrating over V...")
tmp = integrate(integrand3, IntU.HaarMeasure(dV_dim, IntU.MetadataMatcher(:V)))
println("Integrating over U...")
res3 = integrate(tmp, dU(dU_dim))
println("Result: ", res3)
println("Expected: X_1_1 / (dU_dim * dV_dim)")


# Example 4: Symbolic vector integration
println("\n--- Example 4: Vector moments ---")
@variables d
U = SymbolicMatrix(:U, :U)
X = SymbolicMatrix(:X)

expr4 = sum(abs(U[i, j])^2 for i = 1:2, j = 1:2)
res4 = integrate(expr4, dU(d))
println("Integrating sum of squares...")
println("Result: ", res4)
println("Expected: 4/d")

println("\nAll examples updated.")

println("\nAll examples implemented.")
