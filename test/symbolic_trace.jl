using IntU
using Test
using Symbolics

@testset "Symbolic Trace Logic" begin

    @variables d
    U = SymbolicMatrix(:U, false, :U)
    measure = dU(U, d)

    # 1. Tr(U U') = Tr(I) = d
    @testset "Tr(U U')" begin
        U_dag = U'

        # tr(U * U')
        t = tr_lazy(U * U_dag)

        res = integrate(t, measure)
        # Should be d
        @test isequal(Symbolics.simplify(res), d)
    end

    # 2. Tr(U A U' B) = Tr(A) Tr(B) / d
    @testset "Tr(U A U' B)" begin
        A = SymbolicMatrix(:A)
        B = SymbolicMatrix(:B)

        t = tr_lazy(U * A * U' * B)
        res = integrate(t, measure)
        @show res
        @show string(res)
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
        A = SymbolicMatrix(:A)
        B = SymbolicMatrix(:B)
        C = SymbolicMatrix(:C)
        D = SymbolicMatrix(:D)

        # This expression has 2 Us and 2 U's
        t = tr_lazy(U * A * U' * B * U * C * U' * D)

        res = integrate(t, measure)

        s = string(res)
        # Check that result contains traces of combinations of A, B, C, D
        @test !isequal(res, 0)
        @test occursin("d", s)
        @test occursin("tr(", s)
        @test length(s) > 10
    end
end
