using IntU
using BenchmarkTools
using Symbolics

println("=== Diagonal Unitary Integration Benchmark ===\n")

@variables d
@variables V[1:20, 1:20]::Complex

measure = dDiagUnitary(V, d)

# Benchmark 1: Small moment
println("--- 1. Order 2 Moment E[|V_11|^2] ---")
expr1 = abs2(V[1, 1])
@btime integrate($expr1, $measure)

# Benchmark 2: Larger multiset
println("\n--- 2. Order 8 Moment E[|V_11|^2 * |V_22|^2 * |V_33|^2 * |V_44|^2] ---")
expr2 = abs2(V[1, 1]) * abs2(V[2, 2]) * abs2(V[3, 3]) * abs2(V[4, 4])
@btime integrate($expr2, $measure)

# Benchmark 3: Unbalanced (should be 0 and fast)
println("\n--- 3. Unbalanced Moment E[V_11^10 * conj(V_11)^9] ---")
expr3 = V[1, 1]^10 * conj(V[1, 1])^9
@btime integrate($expr3, $measure)

# Comparison with Full Haar U(d) for same indices (computationally much harder for Haar)
println("\n--- 4. Comparison: Order 4 Moment in Haar dU vs dDiag ---")
@variables U[1:20, 1:20]::Complex
measure_u = dU(U, d)
expr_comp = abs2(V[1,1]) * abs2(V[2,2])
expr_u = abs2(U[1,1]) * abs2(U[2,2])

println("Diagonal:")
@btime integrate($expr_comp, $measure)
println("Full Haar:")
@time integrate(expr_u, measure_u)

println("\nDone.")
