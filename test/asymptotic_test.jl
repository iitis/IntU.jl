using IntU
using Test
using Symbolics

@testset "Asymptotic Expansions" begin
    @variables d
    @variables U[1:2, 1:2]::Complex
    
    @testset "Haar Measure: |u11|^2" begin
        expr = U[1,1] * conj(U[1,1])
        res = asymptotic(expr, dU(U, d), 1)
        # 1/d
        @test iszero(Symbolics.simplify(res - 1/d))
    end

    @testset "Haar Measure: |u11|^4" begin
        expr = (U[1,1] * conj(U[1,1]))^2
        # Exact: 2 / (d(d + 1))
        # Asymptotic order 2: 2/d^2
        res2 = asymptotic(expr, dU(U, d), 2)
        @test iszero(Symbolics.simplify(res2 - 2/d^2))
        
        # Asymptotic order 4: 2/d^2 - 2/d^3 + 2/d^4
        res4 = asymptotic(expr, dU(U, d), 4)
        @test iszero(Symbolics.simplify(res4 - (2/d^2 - 2/d^3 + 2/d^4)))
    end

    @testset "Pure States: |psi1|^2" begin
        @variables psi[1:2]::Complex
        expr = psi[1] * conj(psi[1])
        # Exact: 1/d
        res = asymptotic(expr, dPsi(psi, d), 1)
        @test iszero(Symbolics.simplify(res - 1/d))
    end
    
    @testset "Numeric d" begin
        expr = U[1,1] * conj(U[1,1])
        # Even if measure has numeric d, asymptotic should return symbolic expression
        res = asymptotic(expr, dU(U, 3), 1)
        @test res isa Symbolics.Num
        @test !(Symbolics.unwrap(res) isa Number)
    end
end
