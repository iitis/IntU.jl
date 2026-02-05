using IntU
using Symbolics
using LinearAlgebra
using BenchmarkTools

function benchmark_minors(d_vals, sizes)
    for d in d_vals
        println("Benchmarking U($d)")
        @variables U[1:d, 1:d]::Complex
        measure = dU(U, d)

        for s in sizes
            if s >= d
                continue
            end
            println("  $s x $s Minor")

            # Create a simple minor: det(U[1:s, 1:s])
            # For s=2: U[1,1]*U[2,2] - U[1,2]*U[2,1]
            if s == 1
                expr = abs(U[1, 1])^2
            elseif s == 2
                expr = abs(U[1, 1]*U[2, 2] - U[1, 2]*U[2, 1])^2
            else
                # Using det for larger minors
                expr = abs(det(U[1:s, 1:s]))^2
            end

            # Warmup
            integrate(expr, measure)

            # Benchmark
            t = @benchmark integrate($expr, $measure)
            display(t)
            println()
        end
    end
end

println("=== Minors Benchmark ===")
benchmark_minors([4, 6], [1, 2])
