# examples/nb_examples.jl
using IntU
using Symbolics
using LinearAlgebra

# Example 1: Basic scalar integrals
println("--- Example 1 ---")
d1 = 3
@variables U1[1:d1, 1:d1]::Complex
m1 = dU(U1, d1)

println("1.1: Integrate |u[1,1]|^2")
res1_1 = integrate(abs(U1[1,1])^2, m1)
println("Result: ", res1_1, " (Expected: 1/3)")

println("\n1.2: Integrate |u[1,1]*u[2,2]|^2")
res1_2 = integrate(abs(U1[1,1]*U1[2,2])^2, m1)
println("Result: ", res1_2, " (Expected: 1/8)")

println("\n1.3: Integrate u[1,1]*u[2,2]*conj(u[1,2]*u[2,1])")
res1_3 = integrate(U1[1,1]*U1[2,2]*conj(U1[1,2]*U1[2,1]), m1)
println("Result: ", res1_3, " (Expected: -1/24)")


# Example 2: Index-based integration
println("\n--- Example 2 ---")
I1 = [1, 1, 1, 2, 2]; J1 = [2, 2, 1, 1, 1]
I2 = [1, 1, 1, 2, 2]; J2 = [2, 1, 1, 2, 1]
d2 = 6

u_idxs = collect(zip(I1, J1))
u_bar_idxs = collect(zip(I2, J2))
res2 = integrate_indices(u_idxs, u_bar_idxs, d2)
println("Result: ", res2, " (Expected: -1/16200)")


# Example 3: Kronecker product
println("\n--- Example 3 ---")
d3 = 2
@variables U3[1:d3, 1:d3]::Complex
m3 = dU(U3, d3)

U3_kron = kron(U3, U3)
integrand3 = kron(U3_kron, conj.(U3_kron))
res3 = integrate(integrand3, m3)
println("Result matrix (size ", size(res3), "):")
# Displaying a small part or just the fact it worked
println("Top-left element: ", res3[1,1])


# Example 4: Multiple unitaries
println("\n--- Example 4 ---")
d4 = 2
@variables U4[1:d4, 1:d4]::Complex
@variables V4[1:d4, 1:d4]::Complex
@variables X[1:d4^2, 1:d4^2]::Complex

mU = dU(U4, d4)
mV = dU(V4, d4)

UV = kron(collect(U4), collect(V4))
# Conjugate transpose of UV
UV_ct = collect(UV')

# integrand = UV * X * UV'
expr4 = UV * collect(X) * UV_ct

# Sequential integration
println("Integrating over V...")
tmp = integrate(expr4, mV)
println("Integrating over U...")
res4 = integrate(tmp, mU)
println("Result size: ", size(res4))
# Display one element
println("res[1,1]: ", res4[1,1])


# Example 5: Symbolic vector integration
println("\n--- Example 5 ---")
d5 = 2
@variables U5[1:d5, 1:d5]::Complex
@variables X5[1:d5^2, 1:d5^2]::Complex
m5 = dU(U5, d5)

# xi = 1/sqrt(d) * vec(U)
xi = (1 // 1) * vec(collect(U5)) # Use 1 instead of 1/sqrt(d) for now
z_mat = transpose(xi) * collect(X5) * conj.(xi)
z = z_mat[1]
zz = z * conj(z)

println("Integrating z*conj(z)...")
res5 = integrate(zz, m5)
println("Result (first few terms): ", first(string(res5), 100), "...")

println("\nAll examples implemented.")
