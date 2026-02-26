using IntU
using Test
using Symbolics
using LinearAlgebra

@testset verbose=true "Pure States Integration" begin
    @variables d
    psi = SymbolicMatrix(:psi, :psi, d)

    @testset "Diagonal Term" begin
        res = integrate(psi[1, 1] * conj(psi[1, 1]), dPsi(d))
        @test is_really_zero(res - 1/d)
    end

    @testset "Off-Diagonal Term" begin
        res = integrate(psi[1, 1] * conj(psi[2, 1]), dPsi(d))
        @test is_really_zero(res)
    end

    @testset "Fidelity Average" begin
        # Use real components to avoid complex variable simplification issues
        @variables r[1:2] i[1:2]
        phi = [r[1] + im*i[1], r[2] + im*i[2]]

        inner_prod = conj(psi[1, 1])*phi[1] + conj(psi[2, 1])*phi[2]
        expr = inner_prod * conj(inner_prod)
        res = integrate(expr, dPsi(d))

        # Expected: sum_j |phi_j|^2 / d
        expected = (phi[1]*conj(phi[1]) + phi[2]*conj(phi[2])) / d
        @test is_really_zero(res - expected)
    end
end
