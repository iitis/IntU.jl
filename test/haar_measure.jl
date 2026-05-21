using IntegrateUnitary
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
        for i = 1:d_val, j = 1:d_val
            val = to_numeric(real(res_matrix[i, j]))
            res_num[i, j] = Float64(val)
        end
        @test res_num ≈ I_mat

        res_matrix_2 = integrate(U' * U, measure)
        res_num_2 = map(x -> Float64(to_numeric(real(x))), res_matrix_2)
        @test res_num_2 ≈ I_mat
    end

    @testset verbose=true "High-Order Single-Entry Moments" begin
        U10 = SymbolicMatrix(:U, :U, 10)
        res20 = integrate(abs(U10[1, 1])^20, dU(10))
        expected20 = factorial(big(10)) // prod(BigInt(10):BigInt(19))
        @test res20 == expected20

        U_sym = SymbolicMatrix(:U, :U, d)
        res20_sym = integrate(abs(U_sym[1, 1])^20, dU(d))
        expected20_sym = factorial(big(10)) / prod(d + k for k = 0:9)
        @test is_really_zero(Symbolics.simplify(res20_sym - expected20_sym))

        res_mixed = integrate(abs(U10[1, 1])^2 * abs(U10[1, 2])^20, dU(10))
        expected_mixed = factorial(big(1)) * factorial(big(10)) // prod(BigInt(10):BigInt(20))
        @test res_mixed == expected_mixed

        res_mixed_sym = integrate(abs(U_sym[1, 1])^2 * abs(U_sym[1, 2])^20, dU(d))
        expected_mixed_sym = factorial(big(1)) * factorial(big(10)) / prod(d + k for k = 0:10)
        @test is_really_zero(Symbolics.simplify(res_mixed_sym - expected_mixed_sym))

        res_mismatch = integrate(U10[1, 1]^10 * conj(U10[2, 2])^10, dU(10))
        @test res_mismatch == 0
    end

    @testset "_try_extract_int edge cases" begin
        @test IntegrateUnitary._try_extract_int(2) == 2
        @test IntegrateUnitary._try_extract_int(Num(2)) == 2
        @test IntegrateUnitary._try_extract_int(2 // 1) == 2
        @test IntegrateUnitary._try_extract_int(2.0) == 2
        @test IntegrateUnitary._try_extract_int(big(10)^15) == 10^15

        @variables d_tei
        @test IntegrateUnitary._try_extract_int(d_tei) === nothing
        @test IntegrateUnitary._try_extract_int(nothing) === nothing
        @test IntegrateUnitary._try_extract_int(3.5) === nothing
        @test IntegrateUnitary._try_extract_int(3 // 2) === nothing
        @test IntegrateUnitary._try_extract_int(big(10)^50) === nothing
    end

    @testset "BigInt Dimensions" begin
        d_big = big(10)^15
        U_big = symbolic_unitary(:U, d_big)
        res = integrate(abs(U_big[1, 1])^2, dU(d_big))
        @test res == 1 // big(10)^15
    end
end
