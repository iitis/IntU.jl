using IntU
using Test
using Symbolics

@testset "IntU Tests" begin
    # Define variables
    d_val = 3
    # Use literal dimensions for macro
    @variables U[1:3, 1:3]::Complex
    measure = dU(U, d_val)

    @testset "Example 1: |u11|^2" begin
        expr = abs(U[1,1])^2
        # Now Symbolics expanding to hypot(re, im)^2 is handled by IntU.
        
        res = integrate(expr, measure)
        @test isequal(res, 1//3)
    end

    @testset "Example 2: |u11 u22|^2" begin
        # |u11 u22|^2 
        expr = abs(U[1,1] * U[2,2])^2
        
        res = integrate(expr, measure)
        @test isequal(res, 1//8)
    end
    
    @testset "Example 3: u11 u22 conj(u12 u21)" begin
        expr = U[1,1] * U[2,2] * conj(U[1,2]) * conj(U[2,1])
        res = integrate(expr, measure)
        @test isequal(res, -1//24)
    end
    
    @testset "Unitarity Check: sum_k u_ik conj(u_jk) = delta_ij" begin
        # Integrate (sum_k u_ik conj(u_jk)) should be delta_ij * Volume?
        # Actually this is an identity U U^dag = I.
        # Integral of Constant 1 = 1.
        # Integral of constant 0 = 0.
        # Let's test term by term.
        # Int sum_k |u_1k|^2 = sum_k Int |u_1k|^2 = d * (1/d) = 1.
        
        sum_val = 0//1
        for k in 1:d_val
            sum_val += integrate(U[1,k] * conj(U[1,k]), measure)
        end
        @test isequal(sum_val, 1)
        
        # Off-diagonal
        sum_off = 0//1
        for k in 1:d_val
            sum_off += integrate(U[1,k] * conj(U[2,k]), measure)
        end
        @test isequal(sum_off, 0)
    end
    
    @testset "Weingarten Function Values" begin
        # Wg(1^2, d) = 1/(d^2-1)
        @test isequal(IntU.Weingarten.weingarten([1,1], 3), 1//8)
        
        # Wg(2, d) = -1/(d(d^2-1))
        @test isequal(IntU.Weingarten.weingarten([2], 3), -1//24)
    end

    @testset "Weingarten Unit Tests" begin
        # Import internal functions for testing if needed, or use IntU.Weingarten prefix
        using IntU.Weingarten: conjugate_partition, character_at_id, schur_polynomial_at_1, murnaghan_nakayama, weingarten

        @testset "conjugate_partition" begin
            @test isequal(conjugate_partition([1]), [1])
            @test isequal(conjugate_partition([2]), [1,1])
            @test isequal(conjugate_partition([1,1]), [2])
            @test isequal(conjugate_partition([2,1]), [2,1])
            @test isequal(conjugate_partition([3,1]), [2,1,1])
            @test isequal(conjugate_partition([4]), [1,1,1,1])
            @test isequal(conjugate_partition(Int[]), Int[])
        end

        @testset "schur_polynomial_at_1 (Dim of U(d) irrep)" begin
            # s_{1}(1^d) = d
            @test isequal(schur_polynomial_at_1([1], 3), 3//1)
            # s_{2}(1^d) = d(d+1)/2 symmetric tensor
            @test isequal(schur_polynomial_at_1([2], 3), 3*4//2) # 6
            # s_{1,1}(1^d) = d(d-1)/2 antisymmetric
            @test isequal(schur_polynomial_at_1([1,1], 3), 3*2//2) # 3
        end

        @testset "character_at_id (Dim of S_n irrep)" begin
            # S3
            # [3] -> 1 (trivial)
            @test isequal(character_at_id([3]), 1)
            # [1,1,1] -> 1 (sign)
            @test isequal(character_at_id([1,1,1]), 1)
            # [2,1] -> 2 (standard)
            @test isequal(character_at_id([2,1]), 2)
            
            # S4
            # [4] -> 1
            @test isequal(character_at_id([4]), 1)
            # [3,1] -> 3
            @test isequal(character_at_id([3,1]), 3)
            # [2,2] -> 2
            @test isequal(character_at_id([2,2]), 2)
            # [2,1,1] -> 3
            @test isequal(character_at_id([2,1,1]), 3)
        end

        @testset "murnaghan_nakayama (Character table values)" begin
            # S3 Character Table
            # Partitions: [3] (id), [2,1] (transposition), [1,1,1] (3-cycle) - Wait:
            # Cycle types correspond to classes.
            # Lambda [3] (Trivial): 1, 1, 1 everywhere.
            @test isequal(murnaghan_nakayama([3], [1,1,1]), 1)
            @test isequal(murnaghan_nakayama([3], [2,1]), 1)
            @test isequal(murnaghan_nakayama([3], [3]), 1)

            # Lambda [1,1,1] (Sign): 1, -1, 1
            @test isequal(murnaghan_nakayama([1,1,1], [1,1,1]), 1)
            @test isequal(murnaghan_nakayama([1,1,1], [2,1]), -1)
            @test isequal(murnaghan_nakayama([1,1,1], [3]), 1)

            # Lambda [2,1] (Standard): 2, 0, -1
            @test isequal(murnaghan_nakayama([2,1], [1,1,1]), 2)
            @test isequal(murnaghan_nakayama([2,1], [2,1]), 0)
            @test isequal(murnaghan_nakayama([2,1], [3]), -1)
        end
        
        @testset "Weingarten Function consistency" begin
             # Check basic property or redundancy
             # Wg([1,1], d) should be 1/(d^2-1)
             d = 3
             @test isequal(weingarten([1,1], d), 1//(d^2-1))
        end
    end
end
