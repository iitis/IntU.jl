using IntU
using BenchmarkTools
using Symbolics

println("=== SU(d) Integration Benchmark ===\n")

@variables d
U = SymbolicMatrix(:U, :U, d)

measure = dSU(d)
measure_u = dU(d)

# Benchmark 1: Small balanced moment (order 2)
println("--- 1. Order 2 Moment E[|U_11|^2] ---")
expr2 = abs2(U[1, 1])
println("Integrating...")
@btime integrate($expr2, $measure)

# Benchmark 2: Order 4 Moment E[|U_11 U_22|^2]
println("\n--- 2. Order 4 Moment E[|U_11 U_22|^2] ---")
expr4 = abs2(U[1, 1]) * abs2(U[2, 2])
println("Integrating...")
@btime integrate($expr4, $measure)

# Benchmark 3: Order 6 Moment (on diagonal)
println("\n--- 3. Order 6 Moment E[|U_11|^6] ---")
expr6 = abs2(U[1, 1])^3
println("Integrating...")
@btime integrate($expr6, $measure)

println("\nDone.")
