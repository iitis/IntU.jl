using IntU
using Symbolics
using LinearAlgebra
using Test

@testset "Kronecker Enhancements" begin
    # 1. Symbolic Dimensions
    @variables d
    U = symbolic_unitary(:U, d)
    K = kron(U, U)
    @test K isa IntU.SymbolicKron
    @test size(K) == (d^2, d^2)

    # 2. Trace Distribution
    # tr(U ⊗ U) = tr(U) * tr(U)
    expr = tr(K)
    @test expr isa IntU.LazyTrace
    @test length(expr.cycles) == 2

    # 3. Multiplication Optimization
    # (U ⊗ U) * (U† ⊗ U†) = (U*U†) ⊗ (U*U†)
    K_adj = adjoint(K)
    prod_K = K * K_adj
    @test prod_K isa IntU.SymbolicKron
    @test prod_K.A isa IntU.SymbolicMatrixProduct

    # tr( (U ⊗ U) * (U† ⊗ U†) ) = tr(U*U†) * tr(U*U†) = d * d = d^2
    res = integrate(tr(prod_K), dU(d))
    @test Symbolics.value(Symbolics.simplify(res - d^2)) == 0

    # 4. Unknown Dimensions Fix (no BoundsError)
    U_un = symbolic_unitary(:U_un, nothing)
    @test size(U_un) == (nothing, nothing)
    A = SymbolicMatrix(:A, :Constant, 3)
    # This shouldn't crash with typemax
    P = U_un * A
    @test size(P) == (nothing, 3)

    # 5. DimensionMismatch detection
    U3 = symbolic_unitary(:U3, 3)
    U2 = symbolic_unitary(:U2, 2)
    # (3x3) * (2x2) should fail getindex if we force it
    P_bad = U3 * U2
    @test_throws DimensionMismatch P_bad[1, 1]
end
