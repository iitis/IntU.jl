using IntU
using BenchmarkTools
using Symbolics

println("=== Diagonal Unitary Integration Benchmark ===\n")

@variables d
D = SymbolicMatrix(:D, :DiagUnitary, d)

measure = dDiagUnitary(d)

# Benchmark 1: Small moment
println("--- 1. Order 2 Moment E[|D_11|^2] ---")
expr1 = abs2(D[1, 1])
@btime integrate($expr1, $measure)

# Benchmark 2: Larger multiset
println("\n--- 2. Order 8 Moment E[|D_11|^2 * |D_22|^2 * |D_33|^2 * |D_44|^2] ---")
expr2 = abs2(D[1, 1]) * abs2(D[2, 2]) * abs2(D[3, 3]) * abs2(D[4, 4])
@btime integrate($expr2, $measure)

# Benchmark 3: Unbalanced (should be 0 and fast)
println("\n--- 3. Unbalanced Moment E[D_11^10 * conj(D_11)^9] ---")
expr3 = D[1, 1]^10 * conj(D[1, 1])^9
@btime integrate($expr3, $measure)

# Comparison with Full Haar U(d) for same indices (computationally much harder for Haar)
println("\n--- 4. Comparison: Order 4 Moment in Haar dU vs dDiag ---")
U = SymbolicMatrix(:U, :U, d)
measure_u = dU(d)
expr_comp = abs2(D[1, 1]) * abs2(D[2, 2])
expr_u = abs2(U[1, 1]) * abs2(U[2, 2])

println("Diagonal:")
@btime integrate($expr_comp, $measure)
println("Full Haar:")
@time integrate(expr_u, measure_u)

println("\nDone.")
