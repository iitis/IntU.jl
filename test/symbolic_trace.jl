using IntU
using Test
using Symbolics

@testset "Symbolic Trace Logic" begin

    @variables d
    U = SymbolicMatrix(:U, :U, d)
    measure = dU(d)

    # 1. Tr(U U') = Tr(I) = d
    @testset "Tr(U U')" begin
        U_dag = U'

        t = tr_lazy(U * U_dag)

        res = integrate(t, measure)
        @test isequal(Symbolics.simplify(res), d)
    end

    # 2. Tr(U A U' B) = Tr(A) Tr(B) / d
    @testset "Tr(U A U' B)" begin
        A = SymbolicMatrix(:A)
        B = SymbolicMatrix(:B)

        t = tr_lazy(U * A * U' * B)
        res = integrate(t, measure)
        s = string(res)
        @test occursin("tr(A)", s)
        @test occursin("tr(B)", s)
        @test occursin("d", s)

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

        t = tr_lazy(U * A * U' * B * U * C * U' * D)

        res = integrate(t, measure)

        s = string(res)
        @test !isequal(res, 0)
        @test occursin("d", s)
        @test occursin("tr(", s)
        @test length(s) > 10
    end
end
