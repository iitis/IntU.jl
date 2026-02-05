using IntU
using Symbolics
using LinearAlgebra
using BenchmarkTools

function benchmark_trace_moments(d_vals, k_vals)
    for d in d_vals
        println("Benchmarking U($d)")
        @variables U[1:d, 1:d]::Complex
        measure = dU(U, d)
        tr_U = IntU.tr(U)

        for k in k_vals
            if k > d
                continue
            end
            expr = abs(tr_U)^(2*k)
            println("  k=$k (Moment $(2*k))")

            # Warmuo
            integrate(expr, measure)

            # Benchmark
            t = @benchmark integrate($expr, $measure)
            display(t)
            println()
        end
    end
end

println("=== Trace Moments Benchmark ===")
benchmark_trace_moments([2, 4], [1, 2, 3])
