using Test
using IntU
using Symbolics
using LinearAlgebra

@testset "Stiefel Manifold V_k(C^d)" begin
    @variables d


    @testset "Normalization" begin
        k_fixed = 2
        d_fixed = 4

        V = SymbolicMatrix(:V, :V, (d_fixed, k_fixed))

        measure = dStiefel(d_fixed, k_fixed)

        # Check E[V' * V] = I_k
        # Check individual entries to avoid collect symtype issues in some Symbolics versions
        for i = 1:k_fixed, j = 1:k_fixed
            val = integrate((V'*V)[i, j], measure)
            @test val ≈ (i == j ? 1 : 0)
        end

        # Check E[V V'] = (k/d) * I_d
        res_outer = integrate(V * V', measure)
        # res_outer is a Matrix{Num}, comparison with Diagonal works if we iterate or use ≈
        for i = 1:d_fixed, j = 1:d_fixed
            @test res_outer[i, j] ≈ (i == j ? k_fixed // d_fixed : 0)
        end
    end

    @testset "Symbolic d Normalization" begin
        d_sym = d # Using the variable d defined above
        k_fixed = 2
        V = SymbolicMatrix(:V, :V, (d_sym, k_fixed))
        measure = dStiefel(d_sym, k_fixed)

        # 1. Scalar integration with symbolic d
        # E[|V_11|^2] = 1/d
        res11 = integrate(abs2(V[1, 1]), measure)
        @test IntU._symbolic_isequal(Symbolics.simplify(res11), 1/d_sym)

        # E[|V_11|^2 * |V_12|^2] = 1/(d(d+1)) for complex
        # Note: Stiefel integration logic currently maps to unitary
        res1112 = integrate(abs2(V[1, 1]) * abs2(V[1, 2]), measure)
        @test IntU._symbolic_isequal(
            Symbolics.simplify(res1112 - 1/(d_sym * (d_sym + 1))),
            0,
        )

        # 2. Verify that direct SymbolicMatrixProduct integration with symbolic d throws ArgumentError
        # (It requires expansion which fails for non-numeric matrix sizes)
        @test_throws ArgumentError integrate(V' * V, measure)
    end

    @testset "Consistency with dPsi (k=1)" begin
        d_sym = d

        # Stiefel version
        V = SymbolicMatrix(:V, :V, (d_sym, 1))
        m_stiefel = dStiefel(d_sym, 1)

        # Pure state version (tag :psi)
        psi = SymbolicMatrix(:psi, :psi, (d_sym, 1))
        m_psi = dPsi(d_sym)

        # If we use V[1,1] with m_stiefel it should match psi[1,1] with m_psi
        res_stiefel = integrate(V[1, 1] * conj(V[1, 1]), m_stiefel)
        res_psi = integrate(psi[1, 1] * conj(psi[1, 1]), m_psi)

        @test isequal(res_stiefel, res_psi)
        @test isequal(res_stiefel, 1/d_sym)
    end

    @testset "Bounds Enforcement" begin
        @variables d k

        # V should be (d, k)
        # Accessing V[1, k+1] should throw BoundsError
        @test_throws BoundsError begin
            @integrate V[1, k+1] dStiefel(d, k)
        end

        # Test that integrate with invalid indices returns 0 (manual bypass of bounds check)
        V_large = SymbolicMatrix(:V, :V, (d, k+1))
        res = integrate(V_large[1, k+1], dStiefel(d, k))
        @test isequal(res, 0)
    end
end
