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
        # If it prints as 0, assume 0 (hacky but works for simplified terms)
        if string(x) == "0"
            return true
        end
        return false
    end

    # Define variables
    @variables d
    @variables O[1:3, 1:3] # symbolic matrix

    # Measure
    m = dO(O, d)

    @testset "2nd Moment" begin
        # int O_ij O_kl = delta_ik delta_jl / d

        # Case 1: i=1, j=1, k=1, l=1 -> 1/d
        expr = O[1, 1] * O[1, 1]
        res = integrate(expr, m)
        @test isequal(res, 1/d)

        # Case 2: i=1, j=1, k=2, l=2 -> 0
        expr = O[1, 1] * O[2, 2]
        res = integrate(expr, m)
        @test isequal(res, 0)

        # Case 3: i=1, j=1, k=1, l=2 -> 0
        expr = O[1, 1] * O[1, 2]
        res = integrate(expr, m)
        @test isequal(res, 0)

        # Case 4: i=1, j=2, k=1, l=2 -> 1/d
        expr = O[1, 2] * O[1, 2]
        res = integrate(expr, m)
        @test isequal(res, 1/d)
    end

    @testset "Orthogonality Check" begin
        # int (O O^T)_ij = delta_ij
        # sum_k O_ik O_jk

        # Symbolic sum
        sum_val = 0
        for k = 1:3
            sum_val += integrate(O[1, k] * O[1, k], m)
        end
        # Check if simplified difference is zero
        diff = Symbolics.simplify(sum_val - 3/d)
        @test is_like_zero(diff)

        # If we set d=3
        m3 = dO(O, 3)
        res3 = 0
        for k = 1:3
            res3 += integrate(O[1, k] * O[1, k], m3)
        end
        @test res3 == 1
    end

    @testset "4th Moment" begin
        # int O_11^4
        expr = O[1, 1]^4
        res = integrate(expr, m)
        expected = 3 / (d * (d + 2))
        diff = Symbolics.simplify(res - expected)
        @test is_like_zero(diff)
    end

    @testset "Symplectic Integration" begin
        @variables S[1:2, 1:2]::Complex
        # Use d=2 (smallest valid Sp(2n))
        mS = dSp(S, 2)

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
        # O * O^T should integrate to I
        # Use collect to ensure Matrix{Num}
        O_conc = collect(O)
        res_mat = integrate(O_conc * O_conc', m)

        # Expect I(3)
        for i = 1:3, j = 1:3
            val = res_mat[i, j]
            expected = (i == j) ? 1.0 : 0.0

            # Simple numeric check for d
            val_sub = Symbolics.substitute(val, Dict(d => 3))
            num_val = to_numeric(val_sub)

            if num_val isa Number
                # Allow slight tolerance or exact match
                @test num_val ≈ expected atol=1e-12
            else
                # Fallback
                diff = Symbolics.simplify(val - expected)
                @test is_like_zero(diff)
            end
        end
    end
end
