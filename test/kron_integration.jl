using IntU
using Symbolics
using LinearAlgebra
using Test

@testset "Kronecker Product Integration" begin
    U = symbolic_unitary(:U, 2)
    A = SymbolicMatrix(:A, :Constant, 4)
    
    # 1. Basic Kronecker product trace
    # tr( (U ⊗ U) A (U† ⊗ U†) ) = tr(A)
    # The fix ensures tr() expands this "dirty" product so it can be integrated piece-wise
    expr = tr(kron(U, U) * A * kron(U', U'))
    res = integrate(expr, dU(2))
    
    # tr(A) = A_1_1 + A_2_2 + A_3_3 + A_4_4
    expected = sum(i -> A[i, i], 1:4)
    @test Symbolics.value(Symbolics.simplify(res - expected)) == 0
    
    # 2. Mixed lazy/expanded products
    # tr( U * B * U† ) where B is a "dirty" matrix should also work if it forms a SymbolicMatrixProduct
    B = [1 0; 0 1] # Identity as a dirty factor
    expr2 = tr(U * B * U')
    res2 = integrate(expr2, dU(2))
    @test isequal(Symbolics.simplify(res2), 2)
    
    # 3. Verify it still behaves for normal SymbolicMatrix
    U3 = symbolic_unitary(:U, 3)
    expr3 = tr(U3 * U3')
    res3 = integrate(expr3, dU(3))
    @test isequal(Symbolics.simplify(res3), 3)
end
