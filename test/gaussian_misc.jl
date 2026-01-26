using IntU
using Symbolics
using Test

@testset "Gaussian Miscellaneous" begin
    @testset "GSE Dimension Parity" begin
        H = SymbolicMatrix(:H)
        # Should work for even
        @test dGSE(H, 2) isa IntU.GSEMeasure
        @test dGSE(H, 4) isa IntU.GSEMeasure
        
        # Should throw for odd
        @test_throws ArgumentError dGSE(H, 1)
        @test_throws ArgumentError dGSE(H, 3)
        
        # Should work for symbolic
        @variables d
        @test dGSE(H, d) isa IntU.GSEMeasure
    end
    
    @testset "GSE Asymptotic" begin
        @variables d_sym
        H = SymbolicMatrix(:H)
        meas = dGSE(H, d_sym)
        
        # <Tr(H^2)>_GSE = d^2 - d
        res2 = asymptotic(IntU.tr(H^2), meas)
        # Expected: d^2 - d
        @test simplify(res2 - (d_sym^2 - d_sym)) == 0
    end
    
    @testset "Scalar Consistency" begin
        # Verify that GUEMeasure and GOEMeasure still work after refactoring
        N = 2
        H_explicit = [Symbolics.variable(:H, i, j) for i in 1:N, j in 1:N]
        
        # GUE
        m_gue = dGUE(H_explicit, N)
        @test integrate(IntU.tr(H_explicit^2), m_gue) == N^2
        
        # GOE
        m_goe = dGOE(H_explicit, N)
        @test integrate(IntU.tr(H_explicit^2), m_goe) == N^2 + N
        
        # GSE
        # For d=2, Tr(H^2) = 2^2 - 2 = 2
        m_gse = dGSE(H_explicit, N)
        @test integrate(IntU.tr(H_explicit^2), m_gse) == 2
    end
end
