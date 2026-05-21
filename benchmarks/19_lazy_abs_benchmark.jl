using IntegrateUnitary
using Symbolics
using BenchmarkTools

println("Benchmarking @integrate abs(tr(U))^n dU(d) with concrete d")

for n in [2, 4, 8, 12, 14]
    println("\nPower n = $n:")
    @btime @integrate abs(tr(U))^$n dU(10)
end

# Note: symbolic-d trace moments |tr(U)|^{2k} with k > 1 are not supported
# (result depends on d as a step function, not a polynomial).
# See src/library.jl for the ArgumentError.
