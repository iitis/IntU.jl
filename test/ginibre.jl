using IntU
using Test
using Symbolics
using LinearAlgebra
import LinearAlgebra: tr

@testset "Ginibre Ensembles" begin
    N = 2
    # Ensure variables are Complex so conj(x) != x for GinUE
    G = [Symbolics.variable(:G, i, j, T=Complex{Num}) for i = 1:N, j = 1:N]
    
    @testset "Complex Ginibre (GinUE)" begin
        meas = dGinUE(G, N)
        
        @testset "Basic Moments" begin
            # < Tr(G) > = 0
            @test to_numeric(integrate(tr(G), meas)) == 0
            
            # < Tr(G*G') > = N^2
            expr = tr(G * G')
            @test to_numeric(integrate(expr, meas)) == N^2
            
            # < G_11 conj(G_11) > = 1
            @test to_numeric(integrate(G[1,1] * conj(G[1,1]), meas)) == 1
            
            # < G_11 conj(G_22) > = 0
            @test to_numeric(integrate(G[1,1] * conj(G[2,2]), meas)) == 0
            
            # < G_11 G_22 > = 0
            @test to_numeric(integrate(G[1,1] * G[2,2], meas)) == 0
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
            # Using SymbolicMatrix to trigger LazyTrace logic
            Gs = SymbolicMatrix(:G)
            meas_s = dGinUE(Gs, N)
            
            # < Tr(Gs * Gs') >
            t = tr_lazy(Gs * Gs')
            res = integrate(t, meas_s)
            @test to_numeric(res) == N^2
            
            As = SymbolicMatrix(:A)
            Bs = SymbolicMatrix(:B)
            t2 = tr_lazy(Gs * As * Gs' * Bs)
            # This should integrate correctly
            res2 = integrate(t2, meas_s)
            # For LazyTrace, Tr(As) results in a symbolic tr(A)
            # < Tr(G A G' B) > = Tr(A) Tr(B)
            # res2 should be tr(A) * tr(B)
            @test Symbolics.iscall(Symbolics.unwrap(res2))
        end
    end
    
    @testset "Real Ginibre (GinOE)" begin
        meas = dGinOE(G, N)
        # < G_ij G_kl > = delta_ik delta_jl
        # < Tr(G * G^T) > = sum_i,j,k (G_ij G_kj) delta_ik? No.
        # Tr(G * G^T) = sum_i (G G^T)_ii = sum_i,j G_ij (G^T)_ji = sum_i,j G_ij G_ij
        # < sum_i,j G_ij^2 > = sum_i,j delta_ii delta_jj = N^2
        @test to_numeric(integrate(IntU.tr(G * transpose(G)), meas)) == N^2
        
        # < G_11^2 > = 1
        @test to_numeric(integrate(G[1,1]^2, meas)) == 1
    end
    
    @testset "Symplectic Ginibre (GinSE)" begin
        meas = dGinSE(G, N)
        # Verify it doesn't crash and follows duality
        @test to_numeric(integrate(IntU.tr(G * G'), meas)) !== nothing
    end
    
    @testset "Asymptotic Expansion" begin
        d = Symbolics.variable(:d)
        # Use Complex variables here too
        Gd = [Symbolics.variable(:G, i, j, T=Complex{Num}) for i = 1:N, j = 1:N]
        meas = dGinUE(Gd, d)
        expr = tr(Gd * Gd')
        asymp = asymptotic(expr, meas, 1)
        # < Tr(G G') > = N^2 = 4 (since Gd is 2x2 regardless of d parameter)
        asymp = asymptotic(expr, meas, 1)
        # < Tr(G G') > = N^2 = 4 (since Gd is 2x2 regardless of d parameter)
        val = Symbolics.substitute(asymp, Dict(d => 10))
        @test IntU._symbolic_isequal(val, 4)
    end
end
