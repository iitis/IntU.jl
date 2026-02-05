# examples/03_symbolic_dimension.jl
using IntU
using Symbolics
using LinearAlgebra

# Define symbolic dimension
@variables d

println("Testing symbolic dimension d with @symbolic_dimension macro...")

# Use the new macro to define U
@symbolic_dimension U[1:d, 1:d]
measure = dU(U)

# 1. Integrate |U_11|^2
println("\n1. Integrating |U_11|^2 ...")
expr1 = abs(U[1, 1])^2
res1 = integrate(expr1, measure)
println("Result: $res1 (Expected: 1/d)")

# 2. Integrate |U_11|^4
println("\n2. Integrating |U_11|^4 ...")
expr2 = abs(U[1, 1])^4
res2 = integrate(expr2, measure)
# Expected Wg([1,1], d) * 2 = 2 / (d*(d+1))
# Wait, for |U_11|^4, n=2. Combinations are (1,1) (1,1) bar(1,1) bar(1,1).
# Weingarten states: <U_11 U_11 U_bar_11 U_bar_11> = sum_{sigma, tau in S_2} delta... Wg(sigma*tau^-1, d)
# For n=2, there are 2 permutations: id and (12).
# result = 2 * (Wg(id, d) + Wg((12), d)) 
# Wg(id, d) = 1/(d^2-1), Wg((12), d) = -1/(d(d^2-1))
# Sum = (d-1)/(d(d^2-1)) = 1/(d(d+1))
# Total = 2 / (d(d+1))
println("Result: $res2 (Expected: 2 / (d*(d+1)))")

# 3. Integrate a 2x2 minor
println("\n3. Integrating |U_11*U_22 - U_12*U_21|^2 ...")
# Convert to Matrix of expressions to force determinant expansion
# Since U is lazy infinite, we must explicitly slice it or construct the minor
# U[1:2, 1:2] usually works if getindex supports ranges, but our lazy implementation
# might only be optimized for scalar indexing implicitly or we can just build the matrix.
U_mat = [U[i, j] for i = 1:2, j = 1:2]
minor = det(U_mat)
expr3 = abs(minor)^2
res3 = integrate(expr3, measure)
# Theoretical expectation: 2 / (d*(d-1))
println("Result: $res3 (Expected: 2 / (d*(d-1)))")

# Simplify results for better readability
println("\nSimplified results:")
println("1. $(Symbolics.simplify(Symbolics.expand(res1)))")
println("2. $(Symbolics.simplify(Symbolics.expand(res2)))")
println("3. $(Symbolics.simplify(Symbolics.expand(res3)))")
