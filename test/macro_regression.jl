using IntU
using Test
using Symbolics

@testset "@integrate Macro Regression" begin
    @variables d
    
    @testset "Unitary Family" begin
        @test isequal(@integrate(abs(U[1,1])^2, dU(d)), 1/d)
        @test isequal(@integrate(abs(U[1,1])^2, dCUE(d)), 1/d)
        @test isequal(@integrate(abs(U[1,1])^2, dSU(d)), 1/d)
    end

    @testset "Orthogonal Family" begin
        @test isequal(@integrate(O[1,1]^2, dO(d)), 1/d)
        # dCOE
        res_coe = @integrate(O[1,1], dCOE(d))
        @test is_really_zero(res_coe)
        res_coe2 = @integrate(O[1,1]^2, dCOE(d))
        @test is_really_zero(res_coe2)
    end

    @testset "Symplectic Family" begin
        # dSp
        res_sp = @integrate(Sp[1,1]^2, dSp(d))
        @test is_really_zero(res_sp)
        # Testing symbolic index fix: Sp[1, d+1] in 2d-dim
        # Note: d must be even for dSp(d) to be formally correct, 2d is always even.
        res_sp2 = @integrate(Sp[1, d+1], dSp(2d))
        @test is_really_zero(res_sp2)
        
        # dCSE
        # Just check it integrates
        res_cse = @integrate(Sp[1,1]^2, dCSE(d))
        @test !isequal(res_cse, Sp[1,1]^2)
    end

    @testset "Pure States" begin
        @test isequal(@integrate(abs(psi[1,1])^2, dPsi(d)), 1/d)
    end

    @testset "Permutations" begin
        # dPerm
        res_p = @integrate(P[1,1], dPerm(d))
        @test isequal(res_p, 1/d)
        # dCPerm
        res_y = @integrate(Y[1,1], dCPerm(d))
        @test is_really_zero(res_y)
    end

    @testset "Diagonal Unitary" begin
        @test isequal(@integrate(abs(V[1,1])^2, dDiagUnitary(d)), 1)
    end

    @testset "Stiefel" begin
        @test isequal(@integrate(abs(V[1,1])^2, dStiefel(d, 1)), 1/d)
    end

    @testset "Gaussian Ensembles" begin
        @test isequal(@integrate(H[1,1]^2, dGUE(d)), 1)
        @test isequal(@integrate(H[1,1]^2, dGOE(d)), 2)
        @test isequal(@integrate(H[1,1]^2, dGSE(2d)), 1)
    end

    @testset "Ginibre Ensembles" begin
        @test isequal(@integrate(abs(G[1,1])^2, dGinUE(d)), 1)
        @test isequal(@integrate(G[1,1]^2, dGinOE(d)), 1)
        @test isequal(@integrate(G[1,1]^2, dGinSE(2d)), 1)
    end
end
