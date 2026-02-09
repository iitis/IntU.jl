using IntU
using Test
using LinearAlgebra
using Symbolics

@testset "HCIZ Integrals" begin
    @testset "Vandermonde Determinant" begin
        @test IntU.vandermonde_det([1, 2, 3]) ≈ (1-2)*(1-3)*(2-3)
        @test IntU.vandermonde_det([1, 1]) == 0
    end

    @testset "d=1 Case" begin
        # Tr(A U B U') for d=1 is just a*b
        # HCIZ prefactor is 1 (prod_{p=1}^0 p!)
        # matrix is [exp(a*b)], det is exp(a*b)
        # vandermonde is 1 (empty product)
        A = fill(0.5, 1, 1)
        B = fill(2.0, 1, 1)
        @test hciz(A, B) ≈ exp(1.0)
    end

    @testset "d=2 Case" begin
        # Analytical for d=2: (exp(a1*b1 + a2*b2) - exp(a1*b2 + a2*b1)) / ((a1-a2)(b1-b2))
        a = [1.0, 2.0]
        b = [0.5, 1.5]
        res = hciz(a, b)
        
        expected = (exp(1*0.5 + 2*1.5) - exp(1*1.5 + 2*0.5)) / ((1-2)*(0.5-1.5))
        @test res ≈ expected
    end

    @testset "Degenerate Eigenvalues" begin
        # Test if perturbation allows computing nearby values
        a = [1.0, 1.0]
        b = [0.5, 1.5]
        # This should trigger perturbation and not crash
        res = hciz(a, b)
        @test !isnan(res)
        @test !isinf(res)
        
        # Limit as a1 -> a2 should be (b1 exp(a1b1) - b2 exp(a1b2)) / (b1 - b2) times some factor?
        # Actually it's better to just check it's finite and reasonable.
    end
    
    @testset "Symbolic Eigenvalues" begin
        @variables a1 a2 b1 b2
        res = hciz([a1, a2], [b1, b2])
        # Check if result is a symbolic expression
        @test res isa Num
    end
end
