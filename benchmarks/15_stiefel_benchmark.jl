using BenchmarkTools
using IntU
using Symbolics
using LinearAlgebra

println("=== Stiefel Manifold Benchmark ===")

@variables d

# Benchmark Setup
# We benchmark integration of polynomials of increasing degree over Stiefel manifold.
k = 3 # Fixed k
V_sym = SymbolicMatrix(:V, :U, (d, k))
measure = dStiefel(d, k)

# Helper to generate random polynomial term of degree 2m
function random_stiefel_poly(V, m)
    poly = 1
    # We use concrete ranges for random index selection since the integration is symbolic in d.
    # The actual values don't matter as long as they are distinct.
    for _ = 1:m
        i, j = rand(1:4), rand(1:3)
        # Add V_ij
        poly *= V[i, j]

        # Add conj(V_pq) to keep it balanced (usually required for non-zero result)
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

    # Pre-compile
    integrate(poly, measure)

    SUITE["stiefel"]["degree_$(2*m)"] = @benchmarkable integrate($poly, $measure)
end

# Also benchmark effect of increasing k (though Weingarten depends mostly on degree)
println("Benchmarking Stiefel integration (degree 4, variable k)...")
poly_fixed = V_sym[1, 1] * conj(V_sym[1, 1]) * V_sym[1, 2] * conj(V_sym[1, 2]) # Degree 4

for k_val in [2, 4, 8]
    # Re-create V and measure for new k
    # We need enough distinct indices
    local V_k = SymbolicMatrix(Symbol("Vk_$k_val"), :U, d)
    # Construct a valid poly for this V_k
    local poly_k = V_k[1, 1] * conj(V_k[1, 1]) * V_k[1, 2] * conj(V_k[1, 2])
    local measure_k = dStiefel(d, k_val)

    integrate(poly_k, measure_k)
    SUITE["stiefel"]["k_$(k_val)"] = @benchmarkable integrate($poly_k, $measure_k)
end
