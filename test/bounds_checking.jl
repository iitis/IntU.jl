using IntU
using Symbolics
using Test

@testset "SymbolicMatrix Bounds Checking" begin
    # 1. Integer dimensions
    U10 = SymbolicMatrix(:U, :Haar, 10)
    @test_throws BoundsError U10[11, 1]
    @test_throws BoundsError U10[1, 11]
    @test_throws BoundsError U10[0, 1]
    @test_throws BoundsError U10[1, 0]
    @test !(U10[5, 5] isa AbstractMatrix)
    
    # 2. Symbolic dimensions
    @variables d
    Ud = SymbolicMatrix(:U, :Haar, d)
    
    # Ud[d+1, 1] should throw BoundsError because (d+1) - d = 1 > 0
    @test_throws BoundsError Ud[d+1, 1]
    @test_throws BoundsError Ud[1, d+1]
    
    # Ud[0, 1] should throw because simplify(0) < 1
    @test_throws BoundsError Ud[0, 1]
    
    # Valid symbolic access
    @test Ud[d, d] !== nothing

    # 3. Stiefel and Pure State specific bounds
    @variables k
    Vk = SymbolicMatrix(:V, :V, (d, k))
    @test_throws BoundsError Vk[1, k+1]
    @test_throws BoundsError Vk[d+1, 1]
    
    psi = SymbolicMatrix(:psi, :psi, (d, 1))
    @test_throws BoundsError psi[1, 2]
    @test_throws BoundsError psi[d+1, 1]

    # 4. Integration with out-of-bounds indices
    # This should fail during expression construction
    @test_throws BoundsError @integrate abs(U[11, 11])^2 dU(10)
    @test_throws BoundsError @integrate abs(U[d+1, d+1])^2 dU(d)
    
    # 5. Success case for integration
    @test (@integrate abs(U[1, 1])^2 dU(d)) !== nothing
    @test (@integrate abs(V[1, 1])^2 dStiefel(d, k)) !== nothing
    @test (@integrate abs(psi[1, 1])^2 dPsi(d)) !== nothing
end
