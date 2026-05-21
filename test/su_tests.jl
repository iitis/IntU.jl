using IntegrateUnitary
using Test
using Symbolics

@testset "SU(d) Integration" begin
    @variables d
    U = SymbolicMatrix(:U, :U, d)

    measure = dSU(d)

    # Should be 1/d
    res1 = integrate(U[1, 1] * conj(U[1, 1]), measure)
    @test IntegrateUnitary._symbolic_isequal(res1, 1/d)

    # Should be 0
    res2 = integrate(U[1, 1], measure)
    @test IntegrateUnitary._symbolic_isequal(res2, 0)

    # E[ |U_{11}|^2 |U_{22}|^2 ]
    mu = dU(d)
    res_haar = integrate(abs2(U[1, 1]) * abs2(U[2, 2]), mu)
    res_su = integrate(abs2(U[1, 1]) * abs2(U[2, 2]), measure)

    @test IntegrateUnitary._symbolic_isequal(res_su, res_haar)

    @test measure.dim === d
end
