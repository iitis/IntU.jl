using IntU
using Symbolics
using BenchmarkTools

println("=== Circular Ensembles Integration Benchmark ===")

@variables d

# --- Warmup ---
println("Warming up...")
@variables S_warm[1:2, 1:2]::Complex
integrate(S_warm[1,1] * conj(S_warm[1,1]), dCOE(S_warm, d))
integrate(S_warm[1,1] * conj(S_warm[1,1]), dCUE(S_warm, d))
integrate(S_warm[1,1] * conj(S_warm[1,1]), dCSE(S_warm, d))

# --- COE Benchmark ---
function benchmark_coe(k, N_matrix)
    println("\n--- COE: E[|S_{1,1}|^{2k}] (Matrix Size=$N_matrix) ---")
    @variables S[1:N_matrix, 1:N_matrix]::Complex
    measure = dCOE(S, d)
    expr = (S[1,1] * conj(S[1,1]))^k
    
    # We display time
    @btime integrate($expr, $measure)
end

# --- CSE Benchmark ---
function benchmark_cse(k, N_matrix)
    println("\n--- CSE: E[|S_{1,1}|^{2k}] (Matrix Size=$N_matrix) ---")
    @variables S[1:N_matrix, 1:N_matrix]::Complex
    measure = dCSE(S, d)
    expr = (S[1,1] * conj(S[1,1]))^k
    
    @btime integrate($expr, $measure)
end

# --- CUE Benchmark ---
function benchmark_cue(k, N_matrix)
    println("\n--- CUE: E[|U_{1,1}|^{2k}] (Matrix Size=$N_matrix) ---")
    @variables U[1:N_matrix, 1:N_matrix]::Complex
    measure = dCUE(U, d)
    expr = (U[1,1] * conj(U[1,1]))^k
    
    @btime integrate($expr, $measure)
end


# Run Benchmarks
N_val = 2

# k=1 (2nd moment)
benchmark_coe(1, N_val)
benchmark_cse(1, N_val)
benchmark_cue(1, N_val)

# k=2 (4th moment)
benchmark_coe(2, N_val)
benchmark_cse(2, N_val)
benchmark_cue(2, N_val)

# k=3 (6th moment) - Complexity grows
benchmark_coe(3, N_val)
# benchmark_cse(3, N_val) # CSE might be slower due to loop counting logic?
benchmark_cue(3, N_val)
