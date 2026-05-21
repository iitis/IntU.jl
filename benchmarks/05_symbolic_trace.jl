using BenchmarkTools
using IntegrateUnitary
using Symbolics

function run_benchmarks()
    @variables d
    U = SymbolicMatrix(:U, false, :U)
    measure = dU(d)
    Ud = U'

    A = SymbolicMatrix(:A)
    B = SymbolicMatrix(:B)
    C = SymbolicMatrix(:C)
    D = SymbolicMatrix(:D)

    println("Benchmarking Symbolic Trace Integration...")

    # 1. First order: Tr(U A U' B)
    t1 = tr_lazy(U * A * Ud * B)
    println("1. Order 1 (2 matrices): Tr(U A U' B)")
    @btime integrate($t1, $measure)

    # 2. Second order: Tr(U A Ud B U C Ud D)
    t2 = tr_lazy(U * A * Ud * B * U * C * Ud * D)
    println("2. Order 2 (4 matrices): Tr(U A U' B U C U' D)")
    @btime integrate($t2, $measure)

    # 3. Third order
    E = SymbolicMatrix(:E)
    F = SymbolicMatrix(:F)
    t3 = tr_lazy(U * A * Ud * B * U * C * Ud * D * U * E * Ud * F)
    println("3. Order 3 (6 matrices)")
    @btime integrate($t3, $measure)

end

run_benchmarks()
