using IntU
using Symbolics
using LinearAlgebra
using BenchmarkTools
import LinearAlgebra: tr

function benchmark_ginibre(N_vals, powers)
    println("=== GinUE Benchmarks (Explicit Matrix) ===")
    for N in N_vals
        println("N = $N")
        G = [Symbolics.variable(:G, i, j, T=Complex{Num}) for i = 1:N, j = 1:N]
        meas = dGinUE(G, N)
        for p in powers
            expr = tr(G * G')^p
            println("  <Tr(G G')^$p>")
            t = @benchmark integrate($expr, $meas)
            display(t)
            println()
        end
    end

    println("\n=== GinUE Benchmarks (Symbolic Dimension) ===")
    @variables d
    Gs = SymbolicMatrix(:G)
    meas_s = dGinUE(Gs, d)
    for p in powers
        # For LazyTrace, it's tr(G G')^p
        expr = tr_lazy(Gs * Gs')^p
        println("  <Tr(G G')^$p> (Symbolic d)")
        t = @benchmark integrate($expr, $meas_s)
        display(t)
        println()
    end
end

println("=== Ginibre Ensembles Benchmark ===")
benchmark_ginibre([2, 4], [1, 2])
