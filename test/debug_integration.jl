using IntU
using Symbolics
using Test

@variables d
@variables U[1:1, 1:1]::Complex
measure = dU(U, d)

expr = abs(U[1,1])^2

function is_really_zero(x)
    res = Symbolics.expand(Symbolics.simplify(x))
    return string(Symbolics.unwrap(res)) == "0"
end

# 1. Check integration directly
println("--- Direct Integration ---")
res_int = integrate(expr, measure)
@show res_int
@test is_really_zero(res_int - 1/d)

# 2. Check |U11|^4 asymptotic expansion
println("\n--- Asymptotic Expansion |U11|^4 ---")
expr4 = abs(U[1,1])^4
res_asymp4 = asymptotic(expr4, measure, 4)
@show res_asymp4
# Verification by cross-multiplication to avoid division issues
diff = Symbolics.simplify(Symbolics.expand(res_asymp4 * d^4) - (2 - 2*d + 2*d^2))
@show diff
@test is_really_zero(diff)
