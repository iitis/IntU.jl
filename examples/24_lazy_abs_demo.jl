using IntU
using Symbolics
using LinearAlgebra

# 1. Basic usage and symbolic dimension
@variables d
println("Example 1: @integrate abs(tr(U))^2 dU(d)")
res1 = @integrate abs(tr(U))^2 dU(d)
println("Result: ", res1)
println()

# 2. Concrete dimension and high power
println("Example 2: @integrate abs(tr(U))^10 dU(10)")
@time res2 = @integrate abs(tr(U))^10 dU(10)
println("Result: ", res2)
println("Expected: 5! = 120 (for n=5, power=10, it is k! where k=power/2)")
# Actually for |tr(U)|^2k, the result is k!
println()

# 3. Square root of trace
println("Example 3: sqrt(tr(U) * tr(U'))")
expr3 = sqrt(tr(U) * conj(tr(U)))
println("Expression: ", expr3)
res3 = @integrate expr3^4 dU(d)
println("Result of (expr3^4) integration: ", res3)
