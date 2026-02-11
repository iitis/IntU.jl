using IntU
using Test
using Symbolics

@testset "SU(d) Integration" begin
    @variables d
    @symbolic_dimension U[1:d, 1:d]

    measure = dSU(U, d)

    # 1. Balanced polynomial: U_{11} conj(U_{11})
    # Should be 1/d, same as U(d)

    res1 = integrate(U[1, 1] * conj(U[1, 1]), measure)
    @test IntU._symbolic_isequal(res1, 1/d)

    # 2. Unbalanced polynomial: U_{11}
    # Should be 0
    res2 = integrate(U[1, 1], measure)
    @test IntU._symbolic_isequal(res2, 0)

    # 3. Two-body moment
    # E[ |U_{11}|^2 |U_{22}|^2 ]
    # = (1/d^2 - 1)

    # For U(d):
    mu = dU(U, d)
    res_haar = integrate(abs2(U[1, 1]) * abs2(U[2, 2]), mu)

    res_su = integrate(abs2(U[1, 1]) * abs2(U[2, 2]), measure)

    @test IntU._symbolic_isequal(res_su, res_haar)

    # 4. Check that it handles symbolic dimension correctly in types
    @test measure.dim === d
end
