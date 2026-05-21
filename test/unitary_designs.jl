using IntU
using Test
using Symbolics


@testset verbose=true "Unitary t-Designs" begin
    d_val = 3
    U = SymbolicMatrix(:U, :U, d_val)

    design2 = dDesign(d_val, 2)

    @testset verbose=true "Degree <= t (t=2)" begin
        # Degree 1 integrals (should match Haar)
        expr1 = abs(U[1, 1])^2
        res1 = integrate(expr1, design2)
        @test to_numeric(real(res1)) ≈ 1/3

        # Degree 2 integrals (should match Haar)
        expr2 = abs(U[1, 1] * U[2, 2])^2
        res2 = integrate(expr2, design2)
        @test to_numeric(real(res2)) ≈ 1/8

        sum_val = 0//1
        for k = 1:d_val
            sum_val += integrate(U[1, k] * conj(U[1, k]), design2)
        end
        @test to_numeric(real(sum_val)) ≈ 1.0
    end

    @testset verbose=true "Degree > t (t=2, degree=3)" begin
        # Degree 3 integral (should fail)
        expr3 = abs(U[1, 1])^6
        @test_throws ArgumentError integrate(expr3, design2)
    end

    design1 = dDesign(d_val, 1)

    @testset verbose=true "1-Design Constraints" begin
        # Degree 1 works
        expr1 = abs(U[1, 1])^2
        res1 = integrate(expr1, design1)
        @test to_numeric(real(res1)) ≈ 1/3

        # Degree 2 fails
        expr2 = abs(U[1, 1])^4
        @test_throws ArgumentError integrate(expr2, design1)
    end
end
