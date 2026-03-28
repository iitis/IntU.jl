using IntU
using Test
using Symbolics
using LinearAlgebra

@testset verbose=true "Pure States Integration" begin
    @variables d
    psi = SymbolicMatrix(:psi, :psi, (d, 1))

    @testset "Diagonal Term" begin
        res = @integrate psi[1, 1] * conj(psi[1, 1]) dPsi(d)
        @test is_really_zero(res - 1/d)
    end

    @testset "Off-Diagonal Term" begin
        res = @integrate psi[1, 1] * conj(psi[2, 1]) dPsi(d)
        @test is_really_zero(res)
    end

    @testset "Fidelity Average" begin
        @variables r[1:2] i[1:2]
        phi = [r[1] + im*i[1], r[2] + im*i[2]]

        expr =
            (conj(psi[1, 1])*phi[1] + conj(psi[2, 1])*phi[2]) *
            (psi[1, 1]*conj(phi[1]) + psi[2, 1]*conj(phi[2]))
        res = @integrate expr dPsi(d)

        # Expected: sum_j |phi_j|^2 / d
        expected = (phi[1]*conj(phi[1]) + phi[2]*conj(phi[2])) / d
        @test is_really_zero(res - expected)
    end

    @testset "Dimension Enforcement" begin
        # psi should be (d, 1)
        @test_throws BoundsError begin
            @integrate psi[1, 2] dPsi(d)
        end

        psi_sq = SymbolicMatrix(:psi, :psi, d)
        res = integrate(psi_sq[1, 2], dPsi(d))
        @test is_really_zero(res)
    end

    @testset "Rectangular Matrix-Valued Integration" begin
        psi = SymbolicMatrix(:psi, :psi, (4, 1))
        res = integrate(psi, dPsi(4))
        @test size(res) == (4, 1)
    end
end
