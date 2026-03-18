using IntU
using Symbolics
using LinearAlgebra
using BenchmarkTools

function benchmark_trace_moments(d_vals, k_vals)
    for d in d_vals
        println("Benchmarking U($d)")
        U_sym = SymbolicMatrix(:U, :U, d)
        measure = dU(d)
        tr_U = tr(U_sym)

        for k in k_vals
            if k > d
                continue
            end
            expr = abs2(tr_U)^k
            println("  k=$k (Moment $(2*k))")

            # Warmuo
            integrate(expr, measure)

            t = @benchmark integrate($expr, $measure)
            display(t)
            println()
        end
    end
end

println("=== Trace Moments Benchmark ===")
benchmark_trace_moments([2, 4], [1, 2, 3])
