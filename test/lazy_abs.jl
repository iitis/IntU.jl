using IntegrateUnitary
using Test
using Symbolics

@testset "Lazy Abs and Sqrt" begin
    @variables d
    U = SymbolicMatrix(:U, :U, d)

    @testset "abs(tr(U))" begin
        expr = abs(tr(U))
        @test expr isa IntegrateUnitary.LazyPower
        @test expr.exponent == 0.5
        # |tr(U)| = (tr(U) * tr(U'))^0.5
    end

    @testset "abs(tr(U))^2 integration" begin
        res = @integrate abs(tr(U))^2 dU(d)
        @test isequal(res, 1)
    end

    @testset "abs(tr(U))^(2k) symbolic d raises error" begin
        @test_throws ArgumentError @integrate abs(tr(U))^4 dU(d)
        @test_throws ArgumentError @integrate abs(tr(U))^6 dU(d)
        @test_throws ArgumentError @integrate abs(tr(U))^8 dU(d)
    end
    
    @testset "sqrt(tr(U))" begin
        expr = sqrt(tr(U))
        @test expr isa IntegrateUnitary.LazyPower
        @test expr.exponent == 0.5
    end

    @testset "abs(tr(U))^4 integration (d=10)" begin
        res = @integrate abs(tr(U))^4 dU(10)
        @test res == 2
    end

    @testset "finite-d trace moments exactness" begin
        U1 = SymbolicMatrix(:U, :U, 1)
        @test integrate(abs(tr(U1))^2, dU(1)) == 1
        @test integrate(abs(tr(U1))^4, dU(1)) == 1
        @test integrate(abs(tr(U1))^6, dU(1)) == 1
        @test integrate(abs(tr(U1))^8, dU(1)) == 1

        U2 = SymbolicMatrix(:U, :U, 2)
        @test integrate(abs(tr(U2))^6, dU(2)) == 5
        @test integrate(abs(tr(U2))^8, dU(2)) == 14

        U3 = SymbolicMatrix(:U, :U, 3)
        @test integrate(abs(tr(U3))^8, dU(3)) == 23
    end
end
