using IntU
using Symbolics
using Test

@testset "SymbolicMatrix Bounds Checking" begin
    # 1. Integer dimensions
    U10 = SymbolicMatrix(:U, :U, 10)
    @test_throws BoundsError U10[11, 1]
    @test_throws BoundsError U10[1, 11]
    @test_throws BoundsError U10[0, 1]
    @test_throws BoundsError U10[1, 0]
    @test !(U10[5, 5] isa AbstractMatrix)

    # 2. Symbolic dimensions
    @variables d
    Ud = SymbolicMatrix(:U, :U, d)

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

    # 6. Num-wrapped dimensions
    U_num = SymbolicMatrix(:U, :U, Num(2))
    @test_throws BoundsError U_num[1, 3]
    @test_throws BoundsError U_num[3, 1]
    @test U_num[1, 1] !== nothing

    # 7. Num-wrapped index vs measure dimension
    @test_throws BoundsError integrate(abs(U_num[1, 3])^2, dU(Num(2)))

    # 8. Num-wrapped dimension integration
    U_num2 = symbolic_unitary(:U, Num(2))
    res_num = integrate(U_num2 * U_num2', dU(Num(2)))
    @test size(res_num) == (2, 2)

    # 9. axes consistent with size for Num dimensions
    A_num = SymbolicMatrix(:A, :Constant, Num(3))
    @test size(A_num) == (3, 3)
    @test axes(A_num) == (Base.OneTo(3), Base.OneTo(3))
end
