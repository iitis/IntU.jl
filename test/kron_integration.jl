using IntegrateUnitary
using Symbolics
using LinearAlgebra
using Test

@testset "Kronecker Product Integration" begin
    U = symbolic_unitary(:U, 2)
    A = SymbolicMatrix(:A, :Constant, 4)

    # 1. Basic Kronecker product trace
    # tr( (U ⊗ U) A (U† ⊗ U†) ) = tr(A)
    expr = tr(kron(U, U) * A * kron(U', U'))
    res = integrate(expr, dU(2))

    expected = sum(i -> A[i, i], 1:4)
    @test Symbolics.value(Symbolics.simplify(res - expected)) == 0

    # 2. Mixed lazy/expanded products
    B = [1 0; 0 1]
    expr2 = tr(U * B * U')
    res2 = integrate(expr2, dU(2))
    @test isequal(Symbolics.simplify(res2), 2)

    U3 = symbolic_unitary(:U, 3)
    expr3 = tr(U3 * U3')
    res3 = integrate(expr3, dU(3))
    @test isequal(Symbolics.simplify(res3), 3)

    Sp = symbolic_symplectic(:Sp, 2)
    B = SymbolicMatrix(:B, :Constant, nothing)
    res_sp = integrate(kron(Sp, Sp) * B * kron(Sp', Sp'), dSp(2))
    @test res_sp isa Matrix
end
