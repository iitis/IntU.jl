# test/qi_helpers.jl
using IntU
using Symbolics
using LinearAlgebra
using Test

@testset "QI Helpers" begin
    @testset "Purity and Fidelity" begin
        @variables rho[1:2, 1:2]::Complex
        @variables sigma[1:2, 1:2]::Complex
        
        @test isequal(purity(rho), tr(rho * rho))
        @test isequal(fidelity(rho, sigma), tr(rho * sigma))
        
        @variables U[1:2, 1:2]::Complex
        measure = dU(U, 2)
        rho_fixed = [1.0 0.0; 0.0 0.0]
        # Use collect to avoid Symbolics.Arr matrix multiplication ambiguity
        U_m = collect(U)
        rho_random = U_m * rho_fixed * U_m'
        
        avg_pur = average_purity(rho_random, measure)
        # Expected: 1 for any d if rho_random is pure
        @test to_numeric(real(avg_pur)) ≈ 1.0
        @test to_numeric(imag(avg_pur)) ≈ 0.0
    end
    
    @testset "Partial Trace" begin
        # Test on a Bell state |phi+> = 1/sqrt(2) (|00> + |11>)
        # rho = [1/2 0 0 1/2; 0 0 0 0; 0 0 0 0; 1/2 0 0 1/2]
        M = [0.5 0 0 0.5; 0 0 0 0; 0 0 0 0; 0.5 0 0 0.5]
        dims = (2, 2)
        
        # Trace out second subsystem
        rho_a = partial_trace(M, dims, 2)
        @test size(rho_a) == (2, 2)
        @test isequal(rho_a, [0.5 0; 0 0.5])
        
        # Trace out first subsystem
        rho_b = partial_trace(M, dims, 1)
        @test size(rho_b) == (2, 2)
        @test isequal(rho_b, [0.5 0; 0 0.5])
        
        # Test on a GHZ state |000> + |111>
        # rho = 1/2 (|000><000| + |000><111| + |111><000| + |111><111|)
        GHZ = zeros(8, 8)
        GHZ[1, 1] = 0.5
        GHZ[1, 8] = 0.5
        GHZ[8, 1] = 0.5
        GHZ[8, 8] = 0.5
        
        # Trace out 3rd qubit
        rho_12 = partial_trace(GHZ, (2, 2, 2), 3)
        @test size(rho_12) == (4, 4)
        # Result should be 1/2 (|00><00| + |11><11|)
        expected_12 = zeros(4, 4)
        expected_12[1, 1] = 0.5
        expected_12[4, 4] = 0.5
        @test isequal(rho_12, expected_12)
    end
end
