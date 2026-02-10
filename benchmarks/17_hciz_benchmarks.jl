using IntU
using BenchmarkTools
using LinearAlgebra

println("=== HCIZ Integral Benchmarks ===\n")

function benchmark_hciz(dims)
    for d in dims
        println("Dimension d = $d")
        A = diagm(randn(d))
        B = diagm(randn(d))
        
        b = @benchmark hciz($A, $B)
        display(b)
        println("\n" * "-"^30)
    end
end

dims = [2, 5, 10, 15, 20]
benchmark_hciz(dims)

println("\nSymbolic HCIZ d=3 Comparison")
using Symbolics
@variables a[1:3] b[1:3]
b_sym = @benchmark hciz($a, $b)
display(b_sym)
println()
