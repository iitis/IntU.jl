using IntU
using Symbolics
using BenchmarkTools

function run_benchmarks()
    @variables d
    U = [Symbolics.variable(:u, i, j) for i in 1:2, j in 1:2]
    measure = dU(U, d)

    println("Benchmarking IntU performance...")

    # Basic integral
    expr1 = abs(U[1,1])^2
    println("\nIntegrating |u₁₁|²:")
    @btime integrate($expr1, $measure)

    # Higher order integral (more permutations)
    # |u11|^4 = u11 * u11 * u11_bar * u11_bar
    expr2 = abs(U[1,1])^4
    println("\nIntegrating |u₁₁|⁴:")
    @btime integrate($expr2, $measure)

    # Matrix integration
    U2 = [U[1,1] U[1,2]; U[2,1] U[2,2]]
    expr3 = U2 * U2'
    println("\nIntegrating U * U' (matrix):")
    @btime integrate($expr3, $measure)

    # Complex case with multiple entries
    expr4 = (U[1,1]*U[2,2]) * conj(U[1,2]*U[2,1])
    println("\nIntegrating u₁₁u₂₂ * conj(u₁₂u₂₁):")
    @btime integrate($expr4, $measure)

    println("\nBenchmarking complete.")
end

run_benchmarks()
