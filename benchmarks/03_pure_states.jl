using IntU
using Symbolics
using BenchmarkTools
using LinearAlgebra

function benchmark_pure_vs_unitary()
    d = 4
    @variables psi[1:d]
    @variables U[1:d, 1:d]
    
    expr = abs(psi[1])^2 * abs(psi[2])^2
    unitary_expr = abs(U[1,1])^2 * abs(U[2,1])^2
    
    measure_psi = dPsi(psi, d)
    measure_u = dU(U, d)
    
    println("Benchmarking dPsi integration...")
    t1 = @benchmark integrate($expr, $measure_psi)
    display(t1)
    
    println("\nBenchmarking equivalent dU integration...")
    t2 = @benchmark integrate($unitary_expr, $measure_u)
    display(t2)
end

benchmark_pure_vs_unitary()
