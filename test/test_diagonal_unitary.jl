using IntU
using Test
using Symbolics

@testset "Diagonal Unitary Integration" begin
    @variables d
    @variables V[1:10, 1:10]::Complex
    
    measure = dDiagUnitary(V, d)

    @testset "Basic Moments" begin
        # E[|V_11|^2] = 1
        @test simplify(integrate(abs(V[1,1])^2, measure)) == 1
        
        # E[V_11 * V_22^*] = 0 (independent phases)
        @test simplify(integrate(V[1,1] * conj(V[2,2]), measure)) == 0
        
        # E[|V_11|^2 * |V_22|^2] = 1
        @test simplify(integrate(abs(V[1,1])^2 * abs(V[2,2])^2, measure)) == 1
        
        # E[V_11^2 * (V_11^*)^2] = 1
        @test simplify(integrate(V[1,1]^2 * conj(V[1,1])^2, measure)) == 1
    end

    @testset "Non-diagonal entries" begin
        # Diagonal unitary matrix entries V_ij are zero for i != j
        # The measure should handle this by returning 0 if i != j is present in the integrand
        @test simplify(integrate(V[1,2] * conj(V[1,2]), measure)) == 0
        @test simplify(integrate(V[1,1] * conj(V[1,2]), measure)) == 0
    end

    @testset "Symbolic Dimension" begin
        # Result should be independent of d (as long as indices are within range)
        @test simplify(integrate(abs(V[1,1])^2, measure)) == 1
    end

    @testset "Higher Order Moments" begin
        # E[V_11^k * (V_11^*)^m] = delta_km
        @test simplify(integrate(V[1,1]^3 * conj(V[1,1])^3, measure)) == 1
        @test simplify(integrate(V[1,1]^3 * conj(V[1,1])^2, measure)) == 0
    end
end
