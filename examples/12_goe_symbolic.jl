using IntU
using Symbolics

println("=== GOE Integration with Symbolic Dimension ===")

@variables d
# For symbolic matrix entries with symbolic size, we use SymbolicMatrix wrapper
H = SymbolicMatrix(d, d, :H)
meas = dGOE(H, d)

println("Measure: dGOE(H, d)")

# 1. Expectation of Tr(H^2)
# Expected: d^2 + d
println("\n1. Integrating Tr(H^2)...")
expr1 = tr(H^2)
res1 = integrate(expr1, meas)
println("Result: ", res1)

# 2. Expectation of Tr(H^4)
# Expected: 2d^3 + 5d^2 + 5d
println("\n2. Integrating Tr(H^4)...")
expr2 = tr(H^4)
res2 = integrate(expr2, meas)
println("Result: ", res2)

# 3. Component-wise symbolic check
# Even if matrix is symbolic, we can look at component moments if we treat them atomicly
# or just use explicit matrices for low d.
# But here H is the SymbolicMatrix wrapper.

# 4. Asymptotic Expansion in 1/d
# Note: For Gaussian, expansions are usually in d, but we can expansion in 1/d 
# to see the leading order terms.
println("\n3. Leading order terms (order d^3) in Tr(H^4):")
# We don't have a direct 'leading_order' helper that is optimized for polynomials,
# but we can see it from the result.

# 5. Comparing GUE vs GOE
println("\n--- Comparison: GUE vs GOE ---")
meas_gue = dGUE(H, d)
res_gue4 = integrate(expr2, meas_gue)
println("GUE <Tr(H^4)>: ", res_gue4, " (Expected: 2d^3 + d)")
println("GOE <Tr(H^4)>: ", res2,     " (Expected: 2d^3 + 5d^2 + 5d)")

println("\nNote: GOE moments have extra sub-leading terms (d^2) due to the real symmetric symmetry.")
