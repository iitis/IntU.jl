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
        H = SymbolicMatrix(:H, :H)

        @test isequal(integrate(IntU.tr(H^2), dGUE(d)), d^2)
        @test isequal(integrate(IntU.tr(H^4), dGUE(d)), 2d^3 + d)
        @test isequal(integrate(IntU.tr(H^6), dGUE(d)), 5d^4 + 10d^2)
    end

    @testset "GOE Moments" begin
        H = SymbolicMatrix(:H, :H)
        @test isequal(integrate(IntU.tr(H^2), dGOE(d)), d^2 + d)
        @test isequal(integrate(IntU.tr(H^4), dGOE(d)), 2d^3 + 5d^2 + 5d)
    end

    @testset "GSE Moments" begin
        H = SymbolicMatrix(:H, :H)
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
end
