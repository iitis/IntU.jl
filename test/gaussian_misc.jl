using IntU
using Symbolics
using Test

@testset "Gaussian Miscellaneous" begin
    @testset "GSE Dimension Parity" begin
        # Should work for even
        @test dGSE(2) isa IntU.GSEMeasure
        @test dGSE(4) isa IntU.GSEMeasure

        # Should throw for odd
        @test_throws ArgumentError dGSE(1)
        @test_throws ArgumentError dGSE(3)

        # Should work for symbolic
        @variables d
        @test dGSE(d) isa IntU.GSEMeasure
    end

    @testset "GSE Asymptotic" begin
        @variables d_sym
        H = SymbolicMatrix(:H, :H, d_sym)
        meas = dGSE(d_sym)

        # <Tr(H^2)>_GSE = d^2 - d
        res2 = asymptotic(IntU.tr(H^2), meas)
        # Expected: d^2 - d
        @test simplify(res2 - (d_sym^2 - d_sym)) == 0
    end

    @testset "Scalar Consistency" begin
        # Verify that GUEMeasure and GOEMeasure still work after refactoring
        N = 2
        H = SymbolicMatrix(:H, :H, N)

        # GUE
        m_gue = dGUE(N)
        @test integrate(IntU.tr(H^2), m_gue) == N^2

        # GOE
        m_goe = dGOE(N)
        @test integrate(IntU.tr(H^2), m_goe) == N^2 + N

        # GSE
        # For d=2, Tr(H^2) = 2^2 - 2 = 2
        m_gse = dGSE(N)
        @test integrate(IntU.tr(H^2), m_gse) == 2
    end
end
