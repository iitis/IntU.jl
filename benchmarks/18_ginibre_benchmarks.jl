using IntegrateUnitary
using Symbolics
using LinearAlgebra
using BenchmarkTools
import LinearAlgebra: tr

function benchmark_ginibre(N_vals, powers)
    println("=== GinUE Benchmarks (Explicit Matrix) ===")
    for N in N_vals
        println("N = $N")
        G_sym = SymbolicMatrix(:G, :GinUE, N)
        meas = dGinUE(N)
        for p in powers
            expr = tr(G_sym * G_sym')^p
            println("  <Tr(G G')^$p>")
            t = @benchmark integrate($expr, $meas)
            display(t)
            println()
        end
    end

    println("\n=== GinUE Benchmarks (Symbolic Dimension) ===")
    @variables d
    Gs = SymbolicMatrix(:G, :GinUE, d)
    meas_s = dGinUE(d)
    for p in powers
        expr = tr_lazy(Gs * Gs')^p
        println("  <Tr(G G')^$p> (Symbolic d)")
        t = @benchmark integrate($expr, $meas_s)
        display(t)
        println()
    end
end

println("=== Ginibre Ensembles Benchmark ===")
benchmark_ginibre([2, 4], [1, 2])
