using IntU
using Symbolics
using BenchmarkTools

function benchmark_asymptotic()
    @variables d
    U = SymbolicMatrix(:U, :U, d)
    measure = dU(d)

    println("=== Asymptotic Expansion Benchmarks ===")

    # Case 1: Fourth moment |U_{11}|^4
    # Exact: 2/(d^2 - 1)
    println("\n1. Fourth Moment |U_{11}|^4 (Order 4)")
    expr1 = abs2(U[1, 1])^2
    t1 = @benchmark asymptotic($expr1, $measure, 4)
    display(t1)

    # Case 2: Eighth moment |U_{11}|^8
    # Exact: 24 / (d(d^2-1)(d^2-4)(d^2-9)) ... roughly
    println("\n2. Eighth Moment |U_{11}|^8 (Order 6)")
    expr2 = abs2(U[1, 1])^4
    t2 = @benchmark asymptotic($expr2, $measure, 6)
    display(t2)

    # Case 3: Pure State Fourth Moment
    println("\n3. Pure State Fourth Moment (Order 4)")
    # Use SymbolicMatrix for psi to get the :psi tag
    psi = SymbolicMatrix(:psi, :psi, d)
    m_psi = dPsi(d)
    expr3 = abs2(psi[1, 1])^2
    t3 = @benchmark asymptotic($expr3, $m_psi, 4)
    display(t3)
end

benchmark_asymptotic()
