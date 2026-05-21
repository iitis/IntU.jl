using IntU
using Symbolics
using LinearAlgebra
using BenchmarkTools

function benchmark_minors(d_vals, sizes)
    for d in d_vals
        println("Benchmarking U($d)")
        U_sym = SymbolicMatrix(:U, :U, d)
        measure = dU(d)

        for s in sizes
            if s >= d
                continue
            end
            println("  $s x $s Minor")

            if s == 1
                expr = abs(U_sym[1, 1])^2
            elseif s == 2
                expr = abs(U_sym[1, 1]*U_sym[2, 2] - U_sym[1, 2]*U_sym[2, 1])^2
            else
                expr = abs(det(U_sym[1:s, 1:s]))^2
            end

            # Warmup
            integrate(expr, measure)

            t = @benchmark integrate($expr, $measure)
            display(t)
            println()
        end
    end
end

println("=== Minors Benchmark ===")
benchmark_minors([4, 6], [1, 2])
