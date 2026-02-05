using IntU
using Test
using Symbolics
using LinearAlgebra

function is_really_zero(x)
    res = Symbolics.simplify(x)
    if iszero(res)
        return true
    end
    res = Symbolics.expand(res)
    if iszero(res)
        return true
    end
    # Robust check for rational zeros
    try
        num = Symbolics.numerator(res)
        if iszero(Symbolics.expand(num))
            return true
        end
    catch
    end
    return string(Symbolics.unwrap(res)) == "0"
end

@testset "Asymptotic Expansions" begin
    # 1. Test basic Weingarten asymptotic expansion Wg(1, d)
    # Wg(1, d) = 1/d
    # Expansion should match exactly 1/d term or series
    @variables d

    @testset "Weingarten Asymptotic" begin
        # Wg([1]) = 1/d
        wg1 = IntU.weingarten([1], d)
        res1 = asymptotic(wg1, d, 2)
        @test is_really_zero(res1 - 1/d)

        # Wg([2]) = -1/(d(d^2-1)) = -1/d^3 - 1/d^5 ...
        wg2 = IntU.weingarten([2], d)
        res2 = asymptotic(wg2, d, 5)
        # Should be -1/d^3 - 1/d^5
        @test is_really_zero(res2 - (-1/d^3 - 1/d^5))

        # Test generic expression
        term = d^2 + 2d + 1
        res3 = asymptotic(term, d, 2)
        @test is_really_zero(res3 - (d^2 + 2d + 1))

        # Test 1/(d-1) = 1/d + 1/d^2 + ...
        term_frac = 1/(d-1)
        res4 = asymptotic(term_frac, d, 2)
        @test is_really_zero(res4 - (1/d + 1/d^2))
    end

    @testset "Integrated Asymptotic" begin
        @variables d
        @variables U[1:1, 1:1]::Complex
        measure = dU(U, d)

        # Integrate |U11|^2
        expr = abs(U[1, 1])^2
        res = asymptotic(expr, measure, 2)
        @test is_really_zero(res - 1/d)

        # Integrate |U11|^4
        expr4 = abs(U[1, 1])^4
        res4_asymp = asymptotic(expr4, measure, 4)
        # ∫ |U11|^4 = 2 / (d(d+1)) = 2/d^2 - 2/d^3 + 2/d^4 ...
        # Verification by cross-multiplication
        diff = Symbolics.simplify(Symbolics.expand(res4_asymp * d^4) - (2 - 2*d + 2*d^2))
        @test is_really_zero(diff)
    end
end
