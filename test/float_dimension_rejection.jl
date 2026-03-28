using IntU
using Test
using Symbolics

@testset "Float Dimension Rejection" begin
    @testset "Public Constructors Reject Raw Floats" begin
        for ctor in (
            dU,
            dSU,
            dO,
            dCUE,
            dCOE,
            dPsi,
            dDiagUnitary,
            dPerm,
            dCPerm,
            dGUE,
            dGOE,
            dGinUE,
            dGinOE,
        )
            @test_throws ArgumentError ctor(2.0)
        end

        for ctor in (dSp, dCSE, dGSE, dGinSE)
            @test_throws ArgumentError ctor(2.0)
            @test_throws ArgumentError ctor(3.0)
        end

        @test_throws ArgumentError dStiefel(4.0, 2)
        @test_throws ArgumentError dStiefel(4, 2.0)

        @test_throws ArgumentError dDesign(4.0, 2)
        @test_throws ArgumentError dDesign(4, 2.0)
    end

    @testset "Public Constructors Reject Wrapped Float Constants" begin
        d2f = Num(2.0)
        d3f = Num(3.0)

        for ctor in (
            dU,
            dSU,
            dO,
            dCUE,
            dCOE,
            dPsi,
            dDiagUnitary,
            dPerm,
            dCPerm,
            dGUE,
            dGOE,
            dGinUE,
            dGinOE,
        )
            @test_throws ArgumentError ctor(d2f)
        end

        for ctor in (dSp, dCSE, dGSE, dGinSE)
            @test_throws ArgumentError ctor(d2f)
            @test_throws ArgumentError ctor(d3f)
        end

        @test_throws ArgumentError dStiefel(d2f, 2)
        @test_throws ArgumentError dStiefel(4, d2f)

        @test_throws ArgumentError dDesign(d2f, 2)
        @test_throws ArgumentError dDesign(4, d2f)
    end

    @testset "Symbolic Inputs Remain Valid" begin
        @variables d k t

        @test dU(d) isa IntU.HaarMeasure
        @test dSU(d) isa IntU.SpecialUnitary
        @test dO(d) isa IntU.OrthogonalMeasure
        @test dSp(d) isa IntU.SymplecticMeasure
        @test dCUE(d) isa IntU.HaarMeasure
        @test dCOE(d) isa IntU.COEMeasure
        @test dCSE(d) isa IntU.CSEMeasure
        @test dPsi(d) isa IntU.PureStateMeasure
        @test dDiagUnitary(d) isa IntU.DiagonalUnitaryMeasure
        @test dPerm(d) isa IntU.PermutationMeasure
        @test dCPerm(d) isa IntU.CenteredPermutationMeasure
        @test dGUE(d) isa IntU.GUEMeasure
        @test dGOE(d) isa IntU.GOEMeasure
        @test dGSE(d) isa IntU.GSEMeasure
        @test dGinUE(d) isa IntU.GinUEMeasure
        @test dGinOE(d) isa IntU.GinOEMeasure
        @test dGinSE(d) isa IntU.GinSEMeasure
        @test dStiefel(d, k) isa IntU.StiefelMeasure
        @test dDesign(d, t) isa IntU.UnitaryDesign
    end

    @testset "Regression: Large Float dU Rejected Early" begin
        @test_throws ArgumentError dU(1.0e10)
    end

    @testset "Integration Guards for Manually Constructed Measures" begin
        U = SymbolicMatrix(:U, :U)
        H = SymbolicMatrix(:H, :GUE)
        V = SymbolicMatrix(:V, :V, (4, 2))

        @test_throws ArgumentError integrate(abs(U[1, 1])^2, IntU.HaarMeasure(2.0, nothing))
        @test_throws ArgumentError integrate(H[1, 1]^2, IntU.GUEMeasure(2.0, nothing))
        @test_throws ArgumentError integrate(
            abs(V[1, 1])^2,
            IntU.StiefelMeasure(4, 2.0, nothing),
        )
        @test_throws ArgumentError integrate(
            abs(U[1, 1])^2,
            IntU.UnitaryDesign(4, 2.0, nothing),
        )
    end

    @testset "Integration Guards Cover LazyTrace and Specialized Paths" begin
        @variables d

        U_lazy = symbolic_unitary(:U_lazy, d)
        H_lazy = SymbolicMatrix(:H_lazy, :GUE, d)
        V_lazy = SymbolicMatrix(:V_lazy, :V, (4, 2))

        @test_throws ArgumentError integrate(
            tr(U_lazy * U_lazy'),
            IntU.HaarMeasure(2.0, nothing),
        )

        @test_throws ArgumentError integrate(
            tr(H_lazy * H_lazy),
            IntU.GUEMeasure(2.0, nothing),
        )

        @test_throws ArgumentError integrate(
            V_lazy' * V_lazy,
            IntU.StiefelMeasure(4, 2.0, nothing),
        )
    end

    @testset "Num-Wrapped Constraint Validation" begin
        # Odd Num dimension rejected for Sp and CSE
        @test_throws ArgumentError dSp(Num(3))
        @test_throws ArgumentError dCSE(Num(3))
        @test dSp(Num(4)) isa IntU.SymplecticMeasure
        @test dCSE(Num(4)) isa IntU.CSEMeasure

        # Stiefel k > d rejected for Num
        @test_throws ArgumentError dStiefel(Num(3), Num(4))
        @test dStiefel(Num(4), Num(2)) isa IntU.StiefelMeasure
    end
end
