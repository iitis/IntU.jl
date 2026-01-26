using IntU
using Symbolics

println("=== GOE Integration (Explicit Matrix) ===")

N = 2
println("Matrix Size N = $N")

# Define a real symmetric symbolic matrix
@variables h11 h12 h22
H_mat = [h11 h12;
         h12 h22]

# Measure for GOE
meas = dGOE(H_mat, N)

println("\n1. Verifying Var(H_ii) = 2 and Var(H_ij) = 1")
# < H_11^2 >
res11 = integrate(h11^2, meas)
println("< H_11^2 > = $res11 (Expected: 2)")

# < H_12^2 >
res12 = integrate(h12^2, meas)
println("< H_12^2 > = $res12 (Expected: 1)")

println("\n2. Trace Moments")
# < Tr(H^2) > = N^2 + N = 4 + 2 = 6
expr2 = IntU.tr(H_mat^2)
res2 = integrate(expr2, meas)
println("< Tr(H^2) > = $res2 (Expected: 6)")

# < Tr(H^4) > = 2N^3 + 5N^2 + 5N = 2(8) + 5(4) + 5(2) = 16 + 20 + 10 = 46
expr4 = IntU.tr(H_mat^4)
res4 = integrate(expr4, meas)
println("< Tr(H^4) > = $res4 (Expected: 46)")

println("\n3. Comparing with a non-symmetric matrix variable names")
# Note: we can use a generic array as long as we define the measure on it.
# The measure enforces the symmetry in the contraction logic.
@variables X[1:N, 1:N]
meas_X = dGOE(X, N)
# Even if X is not symmetric symbolically, < X_12 X_21 > will be 1 
# because H_atomic_lookup maps both to the same contraction rules.
res_sym = integrate(X[1,2]*X[2,1], meas_X)
println("< X_12 * X_21 >_GOE = $res_sym (Expected: 1)")

println("\nDone.")
