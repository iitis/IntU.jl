using IntU
using Test
using Symbolics

@testset "Integral Library" begin
    @variables d

    @testset "Haar Unitary Trace" begin
        U = SymbolicMatrix(:U, :U)
        A = SymbolicMatrix(:A)
        B = SymbolicMatrix(:B)

        expr = IntU.tr(U * A * U' * B)
        res = integrate(expr, dU(d))

        # Expected: (tr(A) * tr(B)) / d
        expected = (IntU.tr_val([A]) * IntU.tr_val([B])) / d
        @test isequal(res, expected)
    end

    @testset "GUE Moments" begin
        H = SymbolicMatrix(:H, :GUE)

        @test isequal(integrate(IntU.tr(H^2), dGUE(d)), d^2)
        @test isequal(integrate(IntU.tr(H^4), dGUE(d)), 2d^3 + d)
        @test isequal(integrate(IntU.tr(H^6), dGUE(d)), 5d^4 + 10d^2)
    end

    @testset "GOE Moments" begin
        H = SymbolicMatrix(:H, :GOE)
        @test isequal(integrate(IntU.tr(H^2), dGOE(d)), d^2 + d)
        @test isequal(integrate(IntU.tr(H^4), dGOE(d)), 2d^3 + 5d^2 + 5d)
    end

    @testset "GSE Moments" begin
        H = SymbolicMatrix(:H, :GSE)
        @test isequal(integrate(IntU.tr(H^2), dGSE(d)), d^2 - d)
        @test isequal(integrate(IntU.tr(H^4), dGSE(d)), 2d^3 - 5d^2 + 5d)
    end

    @testset "Fallback Check" begin
        # Something not in library
        U_sym = SymbolicMatrix(:U, :U)
        expr = U_sym[1, 1] * conj(U_sym[1, 1])
        res = integrate(expr, dU(d))
        @test isequal(res, 1/d)
    end

    @testset "Prefactor Handling" begin
        # Verify that check_gaussian_library correctly applies the prefactor
        H = SymbolicMatrix(:H, :GUE)
        # 3 * tr(H^2) should yield 3 * d^2, not d^2
        expr = IntU.LazyTrace(Vector{AbstractMatrix}[[H, H]], Num(3))
        res = IntU.check_gaussian_library(expr, dGUE(d), :GUE)
        @test res !== nothing
        @test isequal(res, 3 * d^2)

        # tr(H^2) with prefactor=1 should yield d^2
        expr1 = IntU.LazyTrace(Vector{AbstractMatrix}[[H, H]], Num(1))
        res1 = IntU.check_gaussian_library(expr1, dGUE(d), :GUE)
        @test isequal(res1, d^2)
    end
end
