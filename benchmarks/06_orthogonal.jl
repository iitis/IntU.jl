using IntegrateUnitary
using Symbolics
using BenchmarkTools

println("=== Orthogonal Group Integration Benchmark ===")

function benchmark_orthogonal(k, d_val)
    O_sym = SymbolicMatrix(:O, :O, d_val)
    measure = dO(d_val)
    # Integral of O[1,1]^(2k)
    expr = O_sym[1, 1]^(2*k)

    println("\nk = $k, d = $d_val (Integral of O[1,1]^$(2*k))")
    @btime integrate($expr, $measure)
end

# Warmup
O_warm_sym = SymbolicMatrix(:O, :O, 2)
integrate(O_warm_sym[1, 1]^2, dO(2))

benchmark_orthogonal(1, 2)
benchmark_orthogonal(1, 4)
benchmark_orthogonal(2, 2)
benchmark_orthogonal(2, 4)
benchmark_orthogonal(3, 4)
