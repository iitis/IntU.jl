using IntegrateUnitary
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
println("Expected: 5! = 120 (for d >= k, |tr(U)|^{2k} gives k!; here k=5, d=10 >= 5)")
println()

# 3. Square root of trace — requires concrete d since the result
#    is |tr(U)|^4 which is a pure trace moment with k=2.
println("Example 3: sqrt(tr(U) * tr(U')) with concrete d")
expr3 = sqrt(tr(U) * conj(tr(U)))
println("Expression: ", expr3)
res3 = @integrate expr3^4 dU(10)
println("Result of (expr3^4) integration over dU(10): ", res3)
