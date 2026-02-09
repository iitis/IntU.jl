using IntU
using Test
using Symbolics
using LinearAlgebra

@testset verbose=true "Haar Measure Integration" begin
    # Define variables
    d_val = 3
    # Use literal dimensions for macro
    @variables U[1:3, 1:3]::Complex
    measure = dU(U, d_val)

    @testset verbose=true "Example 1: |u11|^2" begin
        expr = abs(U[1, 1])^2
        res = integrate(expr, measure)
        @test to_numeric(real(res)) ≈ 1/3
        @test to_numeric(imag(res)) ≈ 0
    end

    @testset verbose=true "Example 2: |u11 u22|^2" begin
        expr = abs(U[1, 1] * U[2, 2])^2
        res = integrate(expr, measure)
        @test to_numeric(real(res)) ≈ 1/8
        @test to_numeric(imag(res)) ≈ 0
    end

    @testset verbose=true "Example 3: u11 u22 conj(u12 u21)" begin
        expr = U[1, 1] * U[2, 2] * conj(U[1, 2]) * conj(U[2, 1])
        res = integrate(expr, measure)
        @test to_numeric(real(res)) ≈ -1/24
        @test to_numeric(imag(res)) ≈ 0
    end

    @testset verbose=true "Unitarity Check: sum_k u_ik conj(u_jk) = delta_ij" begin
        sum_val = 0//1
        for k = 1:d_val
            sum_val += integrate(U[1, k] * conj(U[1, k]), measure)
        end
        @test to_numeric(real(sum_val)) ≈ 1.0
        @test abs(to_numeric(imag(sum_val))) < 1e-12

        # Off-diagonal
        sum_off = 0//1
        for k = 1:d_val
            sum_off += integrate(U[1, k] * conj(U[2, k]), measure)
        end
        @test abs(to_numeric(real(sum_off))) < 1e-12
        @test abs(to_numeric(imag(sum_off))) < 1e-12
    end

    @testset "Matrix Integration Checks" begin
        # Example 4: Matrix integration U * U' -> I
        # We collect to ensure we are working with a Matrix{Num} and not a lazy Symbolic Array wrapper
        # which might cause dispatch or iteration issues in some versions of Symbolics.
        expr_mat = collect(U * U')
        res_matrix = integrate(expr_mat, measure)
        # res_matrix should be a 3x3 Matrix of numbers
        @test res_matrix isa AbstractMatrix
        @test size(res_matrix) == (3, 3)

        # Expected result is Identity
        I_mat = Matrix{Complex{Rational{Int}}}(I, 3, 3)
        # Simplify elementwise
        res_simp = map(x -> to_numeric(real(x)) + im*to_numeric(imag(x)), res_matrix)
        @test res_simp ≈ I_mat
        
        # Example 5: U' * U -> I
        expr_mat_2 = collect(U' * U)
        res_matrix_2 = integrate(expr_mat_2, measure)
        res_simp_2 = map(x -> to_numeric(real(x)) + im*to_numeric(imag(x)), res_matrix_2)
        @test res_simp_2 ≈ I_mat
    end
end
