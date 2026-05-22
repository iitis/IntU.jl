using IntegrateUnitary
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

    @testset "tr() rejects non-square matrices" begin
        psi = SymbolicMatrix(:psi, :psi, (4, 1))
        @test_throws ArgumentError tr(psi)

        rect = SymbolicMatrix(:R, :Constant, (2, 3))
        @test_throws ArgumentError tr(rect)

        # Product yielding non-square
        A_ns = SymbolicMatrix(:A, :Constant, (2, 3))
        B_ns = SymbolicMatrix(:B, :Constant, (3, 4))
        @test_throws ArgumentError tr(A_ns * B_ns)
    end

    @testset "show(::SymbolicMatrix) conjugation marker" begin
        U2 = SymbolicMatrix(:U, :U, 2)
        s_U = sprint(show, U2)
        s_adj = sprint(show, U2')
        @test s_U != s_adj
        @test occursin("ᴴ", s_adj) || occursin("'", s_adj)
        @test string(U2) != string(conj(U2))
    end

    @testset "Product indexing expands with Num dimensions" begin
        A_num = SymbolicMatrix(:A, :Constant, Num(2))
        B_num = SymbolicMatrix(:B, :Constant, Num(2))
        val = (A_num * B_num)[1, 1]
        @test !occursin("sum_", string(val))
    end

    @testset "Product dimension mismatch detection" begin
        A_bad = SymbolicMatrix(:A, :Constant, (2, 3))
        B_bad = SymbolicMatrix(:B, :Constant, (4, 2))
        @test_throws DimensionMismatch size(A_bad * B_bad)
    end
end
