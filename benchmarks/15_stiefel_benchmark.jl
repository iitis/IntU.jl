using BenchmarkTools
using IntU
using Symbolics
using LinearAlgebra

println("=== Stiefel Manifold Benchmark ===")

@variables d

k = 3
V_sym = SymbolicMatrix(:V, :V, (d, k))
measure = dStiefel(d, k)

function random_stiefel_poly(V, m)
    poly = 1
    for _ = 1:m
        i, j = rand(1:4), rand(1:3)
        poly *= V[i, j]

        p, q = rand(1:4), rand(1:3)
        poly *= conj(V[p, q])
    end
    return poly
end

SUITE = BenchmarkGroup()
SUITE["stiefel"] = BenchmarkGroup()

println("Benchmarking Stiefel integration (degree 2m, k=$k)...")

for m in [2, 3, 4]
    println("  Benchmark for degree 2m=$(2*m)...")
    poly = random_stiefel_poly(V_sym, m)

    integrate(poly, measure)

    SUITE["stiefel"]["degree_$(2*m)"] = @benchmarkable integrate($poly, $measure)
end

println("Benchmarking Stiefel integration (degree 4, variable k)...")
poly_fixed = V_sym[1, 1] * conj(V_sym[1, 1]) * V_sym[1, 2] * conj(V_sym[1, 2])

for k_val in [2, 4, 8]
    local V_k = SymbolicMatrix(Symbol("Vk_$k_val"), :V, (d, k_val))
    local poly_k = V_k[1, 1] * conj(V_k[1, 1]) * V_k[1, 2] * conj(V_k[1, 2])
    local measure_k = dStiefel(d, k_val)

    integrate(poly_k, measure_k)
    SUITE["stiefel"]["k_$(k_val)"] = @benchmarkable integrate($poly_k, $measure_k)
end
