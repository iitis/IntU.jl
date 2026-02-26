using IntU
using Symbolics
using BenchmarkTools

println("=== Symplectic Group Integration Benchmark ===")

function benchmark_symplectic(k, d_val)
    S_sym = SymbolicMatrix(:S, :Sp, d_val)
    measure = dSp(d_val)
    # Integral of S[1,1]*S[2,2]... effectively k terms or so
    # For Sp, let's use a non-zero integral like S[1, m+1] * S[m+1, 1]
    m = div(d_val, 2)
    expr = S_sym[1, m+1] * S_sym[m+1, 1]

    if k > 1
        # Add more terms to increase complexity
        for i = 2:k
            expr *= S_sym[i, m+i] * S_sym[m+i, i]
        end
    end

    println("\nk = $k, d = $d_val (Product of $k Symplectic pairs)")
    @btime integrate($expr, $measure)
end

# Warmup
S_warm_sym = SymbolicMatrix(:S, :Sp, 2)
integrate(S_warm_sym[1, 2]*S_warm_sym[2, 1], dSp(2))

benchmark_symplectic(1, 2)
benchmark_symplectic(1, 4)
benchmark_symplectic(2, 4)
benchmark_symplectic(1, 6)
benchmark_symplectic(2, 6)
