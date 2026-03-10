using IntU
using Test
using Symbolics

@testset "Gaussian Symbolic Trace Integration" begin
    @variables d

    @testset "GUE Symbolic Trace" begin
        H = SymbolicMatrix(:H, :GUE)
        # tr(H^2)
        expr = tr(H^2)
        res = integrate(expr, dGUE(d))
        @test isequal(simplify(res), simplify(d^2))

        # Product of traces
        expr2 = tr(H) * tr(H)
        res2 = integrate(expr2, dGUE(d))
        @test isequal(simplify(res2), simplify(d))
    end

    @testset "GOE Symbolic Trace" begin
        H = SymbolicMatrix(:H, :GOE)
        # tr(H^2)
        expr = tr(H^2)
        res = integrate(expr, dGOE(d))
        @test isequal(simplify(res), simplify(d^2 + d))
    end

    @testset "GSE Symbolic Trace" begin
        H = SymbolicMatrix(:H, :GSE)
        # tr(H^2)
        expr = tr(H^2)
        res = integrate(expr, dGSE(d))
        # <tr(H^2)>_GSE(d) = d^2 - d
        @test isequal(simplify(res), simplify(d^2 - d))
    end

    @testset "GinUE Symbolic Trace" begin
        G = SymbolicMatrix(:G, :GinUE)
        # tr(G * G')
        expr = tr(G * G')
        res = integrate(expr, dGinUE(d))
        @test isequal(simplify(res), simplify(d^2))

        # tr(G) * tr(G')
        expr2 = tr(G) * tr(G')
        res2 = integrate(expr2, dGinUE(d))
        @test isequal(simplify(res2), simplify(d))
    end

    @testset "GinOE Symbolic Trace" begin
        G = SymbolicMatrix(:G, :GinOE)
        # tr(G * G')
        expr = tr(G * G')
        res = integrate(expr, dGinOE(d))
        # tr(G G^T)
        @test isequal(simplify(res), simplify(d^2))

        # tr(G^2)
        expr2 = tr(G^2)
        res2 = integrate(expr2, dGinOE(d))
        @test isequal(simplify(res2), simplify(d))
    end
end
