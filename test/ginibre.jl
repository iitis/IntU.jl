using IntU
using Test
using Symbolics
using LinearAlgebra
import LinearAlgebra: tr

@testset "Ginibre Ensembles" begin
    N = 2
    G = SymbolicMatrix(:G, :GinUE, N)

    @testset "Complex Ginibre (GinUE)" begin
        meas = dGinUE(N)

        @testset "Basic Moments" begin
            # < Tr(G) > = 0
            @test to_numeric(integrate(tr(G), meas)) == 0

            # < Tr(G*G') > = N^2
            expr = tr(G * G')
            @test to_numeric(integrate(expr, meas)) == N^2

            # < G_11 conj(G_11) > = 1
            @test to_numeric(integrate(G[1, 1] * conj(G[1, 1]), meas)) == 1

            # < G_11 conj(G_22) > = 0
            @test to_numeric(integrate(G[1, 1] * conj(G[2, 2]), meas)) == 0

            # < G_11 G_22 > = 0
            @test to_numeric(integrate(G[1, 1] * G[2, 2], meas)) == 0
        end

        @testset "Matrix Integrals" begin
            # < G * A * G' * B > = Tr(A) * B
            A = [1 2; 3 4]
            B = [5 6; 7 8]
            res = integrate(G * A * G' * B, meas)
            expected = tr(A) * B
            @test all(to_numeric.(res) .== expected)
        end

        @testset "Graphical Calculus (LazyTrace)" begin
            Gs = SymbolicMatrix(:G, :GinUE, N)
            meas_s = dGinUE(N)

            # < Tr(Gs * Gs') >
            t = tr_lazy(Gs * Gs')
            res = integrate(t, meas_s)
            @test to_numeric(res) == N^2

            A = [1 0; 0 1]
            B = [1 0; 0 1]
            t2 = tr_lazy(Gs * A * Gs' * B)
            res2 = integrate(t2, meas_s)
            @test to_numeric(res2) == tr(A) * tr(B)
        end
    end

    @testset "Real Ginibre (GinOE)" begin
        G_oe = SymbolicMatrix(:G, :GinOE, N)
        meas = dGinOE(N)
        # < Tr(G * G^T) > = N^2
        @test to_numeric(integrate(IntU.tr(G_oe * transpose(G_oe)), meas)) == N^2

        # < G_11^2 > = 1
        @test to_numeric(integrate(G_oe[1, 1]^2, meas)) == 1
    end

    @testset "Symplectic Ginibre (GinSE)" begin
        G_se = SymbolicMatrix(:G, :GinSE, N)
        meas = dGinSE(N)
        @test to_numeric(integrate(IntU.tr(G_se * G_se'), meas)) !== nothing
    end

    @testset "Asymptotic Expansion" begin
        d_val = Symbolics.variable(:d)
        Gd = SymbolicMatrix(:G, :GinUE, d_val)
        meas = dGinUE(d_val)
        expr = tr(Gd * Gd')

        asymp = asymptotic(expr, meas, 1)
        val = Symbolics.substitute(asymp, Dict(d_val => 10))
        @test IntU._symbolic_isequal(val, 100)
    end
end
