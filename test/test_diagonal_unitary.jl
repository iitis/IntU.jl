using IntegrateUnitary
using Test
using Symbolics

@testset "Diagonal Unitary Integration" begin
    @variables d
    D = SymbolicMatrix(:D, :DiagUnitary, d)

    measure = dDiagUnitary(d)

    @testset "Basic Moments" begin
        # E[|V_11|^2] = 1
        @test simplify(integrate(abs(D[1, 1])^2, measure)) == 1

        # E[V_11 * V_22^*] = 0
        @test simplify(integrate(D[1, 1] * conj(D[2, 2]), measure)) == 0

        # E[|V_11|^2 * |V_22|^2] = 1
        @test simplify(integrate(abs(D[1, 1])^2 * abs(D[2, 2])^2, measure)) == 1

        # E[V_11^2 * (V_11^*)^2] = 1
        @test simplify(integrate(D[1, 1]^2 * conj(D[1, 1])^2, measure)) == 1
    end

    @testset "Non-diagonal entries" begin
        @test simplify(integrate(D[1, 2] * conj(D[1, 2]), measure)) == 0
        @test simplify(integrate(D[1, 1] * conj(D[1, 2]), measure)) == 0
    end

    @testset "Symbolic Dimension" begin
        @test simplify(integrate(abs(D[1, 1])^2, measure)) == 1
    end

    @testset "Higher Order Moments" begin
        @test simplify(integrate(D[1, 1]^3 * conj(D[1, 1])^3, measure)) == 1
        @test simplify(integrate(D[1, 1]^3 * conj(D[1, 1])^2, measure)) == 0
    end
end
