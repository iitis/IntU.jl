using IntU
using Symbolics
using BenchmarkTools

println("=== Symplectic Group Integration Benchmark ===")

function benchmark_symplectic(k, d_val)
    @variables S[1:d_val, 1:d_val]::Complex
    measure = dSp(S, d_val)
    # Integral of S[1,1]*S[2,2]... effectively k terms or so
    # For Sp, let's use a non-zero integral like S[1, m+1] * S[m+1, 1]
    m = div(d_val, 2)
    expr = S[1, m+1] * S[m+1, 1]
    
    if k > 1
        # Add more terms to increase complexity
        for i in 2:k
            expr *= S[i, m+i] * S[m+i, i]
        end
    end
    
    println("\nk = $k, d = $d_val (Product of $k Symplectic pairs)")
    @btime integrate($expr, $measure)
end

# Warmup
@variables S_warm[1:2, 1:2]::Complex
integrate(S_warm[1,2]*S_warm[2,1], dSp(S_warm, 2))

benchmark_symplectic(1, 2)
benchmark_symplectic(1, 4)
benchmark_symplectic(2, 4)
benchmark_symplectic(1, 6)
benchmark_symplectic(2, 6)
