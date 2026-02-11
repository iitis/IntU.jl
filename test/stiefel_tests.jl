using Test
using IntU
using Symbolics
using LinearAlgebra

@testset "Stiefel Manifold V_k(C^d)" begin
    @variables d

    # Cases to test:
    # 1. Normalization E[V'V] = I_k
    # 2. Covariance of entries
    # 3. Connection to Haar states (k=1)
    # 4. Connection to Haar unitaries (k=d)

    @testset "Normalization" begin
        k = 3
        # Create a symbolic matrix V of size d x k
        # We can simulate this by creating a discrete matrix of variables
        @variables V_discrete[1:4, 1:k]
        # But we need integration with symbolic d.
        # IntU supports direct integration of array expressions.

        # Let's define V elements symbolically
        # We need a way to represent V explicitly for now as IntU needs explicit arrays or SymbolicUnitary
        # StiefelMeasure currently mirrors fallback_integrate which takes Array{Num}

        # We'll use a small fixed d for explicit matrix construction test
        d_fixed = 4
        k_fixed = 2
        # Use explicit symbol construction to ensure simple Num entries
        V = [
            Symbolics.variable(Symbol("V_$(i)_$(j)"), T = Complex{Num}) for
            i = 1:d_fixed, j = 1:k_fixed
        ]

        measure = dStiefel(V, d_fixed, k_fixed)

        # Check E[V' * V] = I_k
        res = integrate(V' * V, measure)
        @test res == I(k_fixed)

        # Check trace(V * V') is NOT I_d, but typically related to projection P
        # E[V V'] = (k/d) * I_d
        res_outer = integrate(V * V', measure)
        @test res_outer == (k_fixed // d_fixed) * I(d_fixed)
    end

    @testset "Symbolic d Normalization" begin
        # Symbolic d test
        d_sym = d # Using the variable d defined above
        k_fixed = 2

        # We cannot create a matrix of size d_sym x k_fixed easily in Julia Symbolics
        # But we can verify individual indices.

        @variables V_11 V_12 V_21 V_22 # Just some elements
        # We cheat and say V is just these relevant elements for the test query
        # But for StiefelMeasure we need to map V_ij -> (i,j)
        # So let's construct a "dummy" V that matches indices we care about

        # Test: E[|V_11|^2] = 1/d
        # To do this, we need to pass a structure that integrate() can parse.
        # The measure takes 'V'.

        # Let's create a proxy struct or just use a small symbolic array and pretend it's part of a larger one?
        # No, the logic loops over the input V array to build the map.
        # So V must contain the symbols we integrate.

        V_small = [
            Symbolics.variable(Symbol("Vs_$(i)_$(j)"), T = Complex{Num}) for
            i = 1:2, j = 1:2
        ]
        # Make sure d is treated as integer or symbolic appropriately
        measure = dStiefel(V_small, d_sym, 2)

        # E[V_11 * conjugate(V_11)] should be 1/d
        res = integrate(V_small[1, 1] * conj(V_small[1, 1]), measure)
        @test isequal(res, 1/d)

        # E[V_11 * conjugate(V_12)] should be 0 (orthogonality of columns)
        res = integrate(V_small[1, 1] * conj(V_small[1, 2]), measure)
        @test isequal(res, 0)

        # E[V_11 * conjugate(V_21)] should be 0 (orthogonality of random vector entries essentially?)
        # Wait, for U: E[U_11 U*_21] = 0. Yes.
        res = integrate(V_small[1, 1] * conj(V_small[2, 1]), measure)
        @test isequal(res, 0)
    end

    @testset "Consistency with dPsi (k=1)" begin
        # Stiefel(d, 1) should be equivalent to dPsi(psi, d)
        d_sym = d
        # Define complex psi explicitly
        psi = [Symbolics.variable(Symbol("psi_$i"), T = Complex{Num}) for i = 1:3]

        # Stiefel version
        V = reshape(psi, 3, 1) # Treat psi as 3x1 matrix
        m_stiefel = dStiefel(V, d_sym, 1)

        # Pure state version
        m_psi = dPsi(psi, d_sym)

        test_poly = abs2(psi[1])

        res_stiefel = integrate(test_poly, m_stiefel)
        res_psi = integrate(test_poly, m_psi)

        @test isequal(res_stiefel, res_psi)
        @test isequal(res_stiefel, 1/d_sym)
    end

    @testset "Consistency with dU (k=d)" begin
        # Stiefel(d, d) should be equivalent to dU(U, d)
        d_fixed = 3
        # U must be complex
        U = [
            Symbolics.variable(Symbol("U_$(i)_$(j)"), T = Complex{Num}) for
            i = 1:d_fixed, j = 1:d_fixed
        ]

        m_stiefel = dStiefel(U, d_fixed, d_fixed)
        m_U = dU(U, d_fixed)

        test_poly = U[1, 1] * conj(U[1, 1]) * U[2, 2] * conj(U[2, 2])

        res_stiefel = integrate(test_poly, m_stiefel)
        res_U = integrate(test_poly, m_U)

        @test isequal(res_stiefel, res_U)
    end

end
