using IntU
using Symbolics
using LinearAlgebra
using BenchmarkTools

function benchmark_matrix_integration(sizes)
    println("=== Matrix Integration Benchmark (Haar U * U') ===")
    @variables d

    for N in sizes
        println("\nBenchmarking Matrix Size N=$N")
        U_sym = SymbolicMatrix(:U, :U, d)
        U = U_sym[1:N, 1:N]
        measure = dU(d)

        expr = collect(U * U')
        println("  Expression size: $(size(expr))")

        integrate(expr, measure)

        t = @benchmark integrate($expr, $measure)
        display(t)
    end
end

benchmark_matrix_integration([2, 3, 4])
