using IntU
using Test
using Symbolics
using LinearAlgebra

@testset "HCIZ Integrals" begin
    @testset "Vandermonde Determinant" begin
        v = [1.0, 2.0, 3.0]
        # (1-2)(1-3)(2-3) = (-1)(-2)(-1) = -2
        @test IntU.vandermonde_det(v) ≈ -2.0
    end

    @testset "Base Case d=1" begin
        a = [1.2]
        b = [0.5]
        # ∫ dU exp(a*U*b*U') = exp(a*b)
        @test hciz(a, b) ≈ exp(1.2 * 0.5)
    end

    @testset "d=2 Case" begin
        a = [1.0, 2.0]
        b = [0.5, 1.5]
        res = hciz(a, b)
        # Analytical: (exp(1*0.5 + 2*1.5) - exp(1*1.5 + 2*0.5)) / ((1-2)(0.5-1.5))
        # = (exp(3.5) - exp(2.5)) / 1
        expected = exp(3.5) - exp(2.5)
        @test res ≈ expected
    end

    @testset "Degenerate Eigenvalues" begin
        a = [1.0, 1.0]
        b = [0.5, 1.5]
        res = hciz(a, b)
        # Should be finite
        @test !isnan(res)
        @test !isinf(res)
        
        # Check near degeneracy
        eps = 1e-8
        a_near = [1.0, 1.0 + eps]
        res_near = hciz(a_near, b)
        @test isapprox(res, res_near, atol=1e-1)
    end

    @testset "Symbolic Eigenvalues" begin
        @variables a1 a2 b1 b2
        res = hciz([a1, a2], [b1, b2])
        # Check that it's a symbolic expression
        @test res isa Num
    end

    @testset "Symbolic Matrix Input (Matrix{Num})" begin
        @variables a1 a2 b1 b2 x y
        @testset "Diagonal" begin
            A = [a1 0; 0 a2]
            B = [b1 0; 0 b2]
            res = hciz(A, B)
            @test res isa Num
        end
        @testset "2x2 Non-Diagonal" begin
            A = [0 x; x 0]
            B = [0 y; y 0]
            res = hciz(A, B)
            @test res isa Num
            # Comparison with manual eigenvalues [x, -x] and [y, -y]
            res_manual = hciz([x, -x], [y, -y])
            # Verify by substitution
            subs = Dict(x => 0.5, y => 0.3)
            val1_eval = eval(Symbolics.toexpr(Symbolics.substitute(res, subs)))
            val2_eval = eval(Symbolics.toexpr(Symbolics.substitute(res_manual, subs)))
            @test val1_eval ≈ val2_eval
        end
    end

    @testset "SymbolicMatrix Input" begin
        A = SymbolicMatrix(:A)
        B = SymbolicMatrix(:B)
        res = hciz(A, B, 2)
        @test res isa Num
        # Check that it contains A_1, A_2 etc
        s = string(res)
        @test occursin("A_1", s)
        @test occursin("B_1", s)
    end
end
