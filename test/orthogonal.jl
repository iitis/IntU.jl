using IntU
using Test
using Symbolics
using LinearAlgebra

@testset "Orthogonal Group Integration" begin
    # Helper to check if a symbolic expression is zero
    function is_like_zero(x)
        x_un = Symbolics.unwrap(x)
        if x_un isa Number
            return iszero(x_un)
        end
        # Try convert
        try
            return iszero(convert(Number, x_un))
        catch
        end
        # If it prints as 0, assume 0
        if string(x) == "0"
            return true
        end
        return false
    end

    # Define variables
    @variables d
    # Use SymbolicMatrix tagged as :O
    O = SymbolicMatrix(:O, :O, d)

    # Measure
    m = dO(d)

    @testset "2nd Moment" begin
        # int O_11^2 = 1/d
        expr = O[1, 1] * O[1, 1]
        res = integrate(expr, m)
        @test isequal(res, 1/d)

        # Off-diagonal
        expr2 = O[1, 1] * O[2, 2]
        res = integrate(expr2, m)
        @test isequal(res, 0)
    end

    @testset "Orthogonality Check" begin
        # sum_k O_ik O_jk = delta_ij
        d_val = 3
        O3 = SymbolicMatrix(:O, :O, d_val)
        m3 = dO(d_val)
        
        sum_val = 0
        for k = 1:d_val
            sum_val += integrate(O3[1, k] * O3[1, k], m3)
        end
        @test is_like_zero(Symbolics.simplify(sum_val - 1))

        # Symbolic sum check
        # We can't do symbolic sum over d directly without summation syntax, 
        # but we can test E[O_11^2 + O_12^2 + O_13^2] = 3/d
        O_sym = SymbolicMatrix(:O, :O, d)
        m_sym = dO(d)
        sum_3 = integrate(O_sym[1, 1]^2 + O_sym[1, 2]^2 + O_sym[1, 3]^2, m_sym)
        @test is_like_zero(Symbolics.simplify(sum_3 - 3/d))
    end

    @testset "4th Moment" begin
        expr = O[1, 1]^4
        res = integrate(expr, m)
        expected = 3 / (d * (d + 2))
        @test is_like_zero(Symbolics.simplify(res - expected))
    end

    @testset "Symplectic Integration" begin
        # Use SymbolicMatrix tagged as :Sp
        S = SymbolicMatrix(:S, :Sp, 2)
        mS = dSp(2)

        # S[1,1]^2 -> 0 per debug
        res1 = integrate(S[1, 1]^2, mS)
        @test is_like_zero(res1)

        # S[1,2]*S[2,1] -> -0.5
        res2 = integrate(S[1, 2]*S[2, 1], mS)
        @test to_numeric(real(res2)) ≈ -0.5

        # |S[1,1]|^2 -> 0.5
        res3 = integrate(abs(S[1, 1])^2, mS)
        @test to_numeric(real(res3)) ≈ 0.5
    end

    @testset "Matrix Integration Checks" begin
        O3 = SymbolicMatrix(:O, :O, 3)
        m3 = dO(3)
        res3 = integrate(O3 * O3', m3)
        @test map(to_numeric, res3) ≈ I(3)
    end
end
