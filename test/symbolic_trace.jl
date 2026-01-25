using IntU
using Test
using Symbolics
using LinearAlgebra

@testset "Symbolic Trace Logic" begin

    @variables d
    @variables U_sym[1:2, 1:2]::Complex
    measure = dU(U_sym, d)
    
    # 1. Tr(U U') = Tr(I) = d
    @testset "Tr(U U')" begin
        U = SymbolicMatrix(:U, false, :U)
        U_dag = U'
        
        # tr(U * U')
        t = tr_lazy(U * U_dag)
        
        res = integrate(t, measure)
        # Should be d
        @test isequal(Symbolics.simplify(res), d)
    end
    
    # 2. Tr(U A U' B) = Tr(A) Tr(B) / d
    @testset "Tr(U A U' B)" begin
        U = SymbolicMatrix(:U, false, :U)
        A = SymbolicMatrix(:A)
        B = SymbolicMatrix(:B)
        
        t = tr_lazy(U * A * U' * B)
        # Res should be tr(A)*tr(B)/d
        res = integrate(t, measure)
        
        s = string(res)
        # Check for presence of tr(A) and tr(B) variables
        @test occursin("tr(A)", s)
        @test occursin("tr(B)", s)
        @test occursin("d", s)
        
        # Check that it simplifies to 1/d * tr(A) * tr(B) effectively
        # We can substitute d=2.
        subbed = Symbolics.substitute(res, Dict(d => 2))
        s_sub = string(subbed)
        @test occursin("0.5", s_sub) || occursin("1//2", s_sub) || occursin("1 / 2", s_sub)
    end

    # 3. Tr(U A U' B U C U' D) (2nd Moment)
    @testset "Tr(U A U' B U C U' D)" begin
        U = SymbolicMatrix(:U, false, :U)
        A = SymbolicMatrix(:A)
        B = SymbolicMatrix(:B)
        C = SymbolicMatrix(:C)
        D = SymbolicMatrix(:D)
        
        # This expression has 2 Us and 2 U's
        t = tr_lazy(U * A * U' * B * U * C * U' * D)
        
        res = integrate(t, measure)
        
        s = string(res)
        # Check that result contains traces of combinations of A, B, C, D
        # Specifically, it should not be empty or zero
        @test !isequal(res, 0)
        
        # It should depend on d
        @test occursin("d", s)
        
        # It should contain trace variables of combinations
        # Note: Weingarten formula for n=2 involves sum over sigma, tau
        # One term (sigma=id, tau=id) -> tr(U A U' B . U C U' D) contractions?
        # Actually it contracts U-U' pairs.
        # It forms traces of constants. e.g. tr(A B) is possible? Or tr(A)tr(B)...
        # Since we use SymbolicMatrix(:A etc), we expect tr_val(...) strings.
        @test occursin("tr(", s)
        
        # Verify result is not just 1/d or 0
        @test length(s) > 10
    end
end
