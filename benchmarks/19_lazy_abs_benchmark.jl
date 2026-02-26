using IntU
using Symbolics
using BenchmarkTools

println("Benchmarking @integrate abs(tr(U))^n dU(d)")

for n in [2, 4, 8, 12, 14]
    println("\nPower n = $n:")
    # Use concrete dimension for fastest path
    @btime @integrate abs(tr(U))^$n dU(10)
end

println("\nBenchmarking symbolic dimension:")
@btime @integrate abs(tr(U))^4 dU(d)
