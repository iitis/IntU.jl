using IntU
using Test
using Symbolics
using LinearAlgebra

@testset "Haar Measure Integration" begin
    d_val = 3
    @variables d
    U = SymbolicMatrix(:U, :U, d_val)
    measure = dU(d_val)

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
        res_matrix = integrate(U * U', measure)
        @test res_matrix isa AbstractMatrix
        @test size(res_matrix) == (d_val, d_val)

        # Expected result is Identity
        I_mat = Matrix(I, d_val, d_val)
        
        # Convert to numeric float for ≈ check
        res_num = zeros(Float64, d_val, d_val)
        for i=1:d_val, j=1:d_val
            val = to_numeric(real(res_matrix[i, j]))
            try
                res_num[i, j] = Float64(val)
            catch e
                println("ERROR converting at ($i, $j): val = $val, type = $(typeof(val))")
                rethrow(e)
            end
        end
        @test res_num ≈ I_mat

        res_matrix_2 = integrate(U' * U, measure)
        res_num_2 = map(x -> Float64(to_numeric(real(x))), res_matrix_2)
        @test res_num_2 ≈ I_mat
    end
end
