using IntegrateUnitary
using Symbolics
using BenchmarkTools
using LinearAlgebra

function benchmark_pure_vs_unitary()
    d = 4
    psi_sym = SymbolicMatrix(:psi, :psi, (d, 1))
    U_sym = SymbolicMatrix(:U, :U, d)

    expr = abs(psi_sym[1, 1])^2 * abs(psi_sym[2, 1])^2
    unitary_expr = abs(U_sym[1, 1])^2 * abs(U_sym[2, 1])^2

    measure_psi = dPsi(d)
    measure_u = dU(d)

    println("Benchmarking dPsi integration...")
    t1 = @benchmark integrate($expr, $measure_psi)
    display(t1)

    println("\nBenchmarking equivalent dU integration...")
    t2 = @benchmark integrate($unitary_expr, $measure_u)
    display(t2)
end

benchmark_pure_vs_unitary()
