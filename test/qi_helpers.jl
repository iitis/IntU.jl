using IntU
using Symbolics
using LinearAlgebra
using Test

@testset verbose=true "QI Helpers" begin
    @testset verbose=true "Purity and Fidelity" begin
        @variables rho[1:2, 1:2]::Complex
        @variables sigma[1:2, 1:2]::Complex

        @test IntU._symbolic_isequal(purity(rho), IntU.tr(rho * rho))
        @test IntU._symbolic_isequal(fidelity(rho, sigma), IntU.tr(rho * sigma))

        U = SymbolicMatrix(:U, :U, 2)
        measure = dU(2)
        rho_fixed = [1.0 0.0; 0.0 0.0]
        rho_random = U * rho_fixed * U'

        avg_pur = average_purity(rho_random, measure)
        # Expected: 1 for any d if rho_random is pure
        @test to_numeric(real(avg_pur)) ≈ 1.0
        @test to_numeric(imag(avg_pur)) ≈ 0.0
    end

    @testset verbose=true "Partial Trace" begin
        # Test on a Bell state |phi+> = 1/sqrt(2) (|00> + |11>)
        # rho = [1/2 0 0 1/2; 0 0 0 0; 0 0 0 0; 1/2 0 0 1/2]
        M = [0.5 0 0 0.5; 0 0 0 0; 0 0 0 0; 0.5 0 0 0.5]
        dims = (2, 2)
        rho = M
        
        # Trace out the second qubit
        rho_A = partial_trace(rho, dims, 2)
        @test map(to_numeric, rho_A) ≈ [0.5 0; 0 0.5]
        
        # Trace out the first qubit
        rho_B = partial_trace(rho, dims, 1)
        @test map(to_numeric, rho_B) ≈ [0.5 0; 0 0.5]
    end
end
