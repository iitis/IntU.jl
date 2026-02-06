using IntU
using Symbolics
using BenchmarkTools

println("=== Permutation Group Integration Benchmark ===")

function benchmark_permutation(k, d_val)
    @variables P[1:d_val, 1:d_val]
    measure = dPerm(P, d_val)
    # Integral of P[1,1] * P[2,2] * ... * P[k,k]
    expr = prod(P[i, i] for i in 1:k)

    println("\nk = $k, d = $d_val (Integral of product of $k diagonal elements)")
    @btime integrate($expr, $measure)
end

function benchmark_centered_permutation(k, d_val)
    @variables Y[1:d_val, 1:d_val]
    measure = dCPerm(Y, d_val)
    # Integral of Y[1,1]^k
    expr = Y[1, 1]^k

    println("\nk = $k, d = $d_val (Integral of Y[1,1]^$k)")
    @btime integrate($expr, $measure)
end

function benchmark_trace(k, d_val)
    @variables P[1:d_val, 1:d_val]
    @variables A[1:d_val, 1:d_val]
    measure = dPerm(P, d_val)
    # Integral of tr(P * A)^k
    expr = Symbolics.scalarize(IntU.tr(P * A))^k

    println("\nk = $k, d = $d_val (Integral of tr(P * A)^$k)")
    @btime integrate($expr, $measure)
end

# Warmup
@variables P_warm[1:2, 1:2]
integrate(P_warm[1, 1], dPerm(P_warm, 2))
@variables Y_warm[1:2, 1:2]
integrate(Y_warm[1, 1], dCPerm(Y_warm, 2))
integrate(Symbolics.scalarize(IntU.tr(P_warm)), dPerm(P_warm, 2))

println("\n--- Permutation Group (S_d) ---")
benchmark_permutation(1, 100)
benchmark_permutation(2, 100)
benchmark_permutation(5, 100)
benchmark_permutation(10, 100)

println("\n--- Centered Permutation Group ---")
benchmark_centered_permutation(1, 10)
benchmark_centered_permutation(2, 10)
benchmark_centered_permutation(3, 10)
benchmark_centered_permutation(4, 10)

println("\n--- Symbolic Trace Benchmarks ---")
benchmark_trace(1, 4)
benchmark_trace(1, 8)
benchmark_trace(2, 4)
