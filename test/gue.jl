using IntegrateUnitary
using Test
using Symbolics
using LinearAlgebra

@testset "GUE Integration" begin
    # Test for N=3
    N = 3
    H = SymbolicMatrix(:H, :GUE, N)
    meas = dGUE(N)

    @testset "Tr(H)" begin
        expr = IntegrateUnitary.tr(H)
        res = integrate(expr, meas)
        # Expected: 0
        @test to_numeric(res) == 0
    end

    @testset "Tr(H^2)" begin
        expr = IntegrateUnitary.tr(H^2)
        res = integrate(expr, meas)
        # Expected: N^2
        @test to_numeric(res) == N^2
    end

    @testset "Tr(H^4)" begin
        expr = IntegrateUnitary.tr(H^4)
        res = integrate(expr, meas)
        # Expected: 2N^3 + N
        expected = 2*N^3 + N
        @test to_numeric(res) == expected
    end

    @testset "Component Moments" begin
        # < H_11^2 > = 1
        res1 = integrate(H[1, 1]^2, meas)
        @test to_numeric(res1) == 1

        # < H_12 H_21 > = 1
        res2 = integrate(H[1, 2]*H[2, 1], meas)
        @test to_numeric(res2) == 1

        # < H_12^2 > = 0
        res3 = integrate(H[1, 2]^2, meas)
        @test to_numeric(res3) == 0
    end

    @testset "Matrix Integration Checks" begin
        # < H > = 0
        res_mat = integrate(H, meas)
        @test all(x -> to_numeric(x) == 0, res_mat)

        # < H^2 > should have diagonals N and off-diagonals 0
        res_mat2 = integrate(H^2, meas)
        expected_diag = N

        for i = 1:N, j = 1:N
            val = to_numeric(res_mat2[i, j])
            if i == j
                @test val == expected_diag
            else
                @test val == 0
            end
        end
    end
end
