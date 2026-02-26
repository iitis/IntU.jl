using IntU
using Symbolics
using LinearAlgebra
using BenchmarkTools

function benchmark_matrix_integration(sizes)
    println("=== Matrix Integration Benchmark (Haar U * U') ===")
    @variables d

    for N in sizes
        println("\nBenchmarking Matrix Size N=$N")
        # Define N x N symbolic matrix
        # Note: integration dimension is symbolic 'd', but matrix size is fixed N
        # Define N x N symbolic matrix with metadata tag :U
        U_sym = SymbolicMatrix(:U, :U, d)
        U = U_sym[1:N, 1:N]
        measure = dU(d)

        # We need to collect to ensure we pass a Matrix{Num} to integrate
        expr = collect(U * U')
        println("  Expression size: $(size(expr))")

        # Warmup
        integrate(expr, measure)

        # Benchmark
        t = @benchmark integrate($expr, $measure)
        display(t)
    end
end

benchmark_matrix_integration([2, 3, 4])
