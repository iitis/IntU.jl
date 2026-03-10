using IntU
using Test
using Symbolics

@testset "Lazy Abs and Sqrt" begin
    @variables d
    U = SymbolicMatrix(:U, :U, d)

    @testset "abs(tr(U))" begin
        expr = abs(tr(U))
        @test expr isa IntU.LazyPower
        @test expr.exponent == 0.5
        # |tr(U)| = (tr(U) * tr(U'))^0.5
    end

    @testset "abs(tr(U))^2 integration" begin
        res = @integrate abs(tr(U))^2 dU(d)
        @test isequal(res, 1)
    end

    @testset "sqrt(tr(U))" begin
        expr = sqrt(tr(U))
        @test expr isa IntU.LazyPower
        @test expr.exponent == 0.5
    end

    @testset "abs(tr(U))^4 integration (d=10)" begin
        res = @integrate abs(tr(U))^4 dU(10)
        @test res == 2
    end
end
