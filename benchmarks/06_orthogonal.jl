using IntU
using Symbolics
using BenchmarkTools

println("=== Orthogonal Group Integration Benchmark ===")

function benchmark_orthogonal(k, d_val)
    @variables O[1:d_val, 1:d_val]
    measure = dO(O, d_val)
    # Integral of O[1,1]^(2k)
    expr = O[1,1]^(2*k)
    
    println("\nk = $k, d = $d_val (Integral of O[1,1]^$(2*k))")
    @btime integrate($expr, $measure)
end

# Warmup
@variables O_warm[1:2, 1:2]
integrate(O_warm[1,1]^2, dO(O_warm, 2))

benchmark_orthogonal(1, 2)
benchmark_orthogonal(1, 4)
benchmark_orthogonal(2, 2)
benchmark_orthogonal(2, 4)
benchmark_orthogonal(3, 4)
