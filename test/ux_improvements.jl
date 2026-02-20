
using IntU
using Test
using Symbolics

@testset "UX Improvements" begin
    @testset "Factory Functions" begin
        U = symbolic_unitary(:U, :d)
        @test U isa SymbolicMatrix
        @test U.special_type == :U
        @test U.dim == :d
        
        O = symbolic_orthogonal(:O, :d)
        @test O.special_type == :O
        
        Sp = symbolic_symplectic(:Sp, :d)
        @test Sp.special_type == :Sp
        
        psi = symbolic_pure_state(:psi, :d)
        @test psi.special_type == :psi
        
        P = symbolic_permutation(:P, :d)
        @test P.special_type == :Perm
    end

    @testset "@integrate Macro" begin
        # Test auto-definition of d and U
        @variables d
        res = @integrate abs(U[1,1])^2 dU(d)
        @test isequal(res, 1/d)
        
        # Test numerical dimension
        res_num = @integrate abs(U[1,1])^2 dU(2)
        @test res_num == 1//2
        
        # Test auto-definition of multiple matrices
        res_tr = @integrate tr(U * A * U' * B) dU(d)
        @test IntU.is_number(evaluate(res_tr, [d => 2, tr(A) => 1, tr(B) => 1]))
    end

    @testset "Evaluate Helper" begin
        @variables d
        expr = 1/d
        @test evaluate(expr, d => 2) == 1//2
        @test evaluate(expr, Dict(d => 3)) == 1//3
    end

    @testset "Show Method" begin
        U = symbolic_unitary(:U, :d)
        buf = IOBuffer()
        show(buf, U)
        @test occursin("U (U)", String(take!(buf)))
        
        A = SymbolicMatrix(:A) # Constant
        show(buf, A)
        @test String(take!(buf)) == "A"
    end
end
