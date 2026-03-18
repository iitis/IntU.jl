using IntU
using Symbolics
using LinearAlgebra
using Test

@testset verbose=true "QI Helpers" begin
    @testset verbose=true "Partial Trace" begin
        # Test on a Bell state |phi+> = 1/sqrt(2) (|00> + |11>)
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
