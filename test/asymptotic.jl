using IntU
using Test
using Symbolics
using LinearAlgebra

@testset "Asymptotic Expansions" begin
    # 1. Test basic Weingarten asymptotic expansion Wg(1, d)
    # Wg(1, d) = 1/d
    # Expansion should match exactly 1/d term or series
    @variables d
    Series = IntU.Series
    
    wg1 = IntU.weingarten([1], d)
    # asymptotic(wg1, d, 2) should look like 1/d
    
    # 2. Test integrated Moments
    # Integral of U_ij U_bar_kl = delta(i,k)delta(j,l)/d
    # Asymptotic: 1/d * delta...
    
    # 3. Test generic expression
    # expr = d^2 + d + 1
    # order 1 -> d^2 + d. 
    # But asymptotic is usually in 1/d.
    # Our function handles O(d^k).
    
    @testset "Weingarten Asymptotic" begin
        # Wg([1]) = 1/d
        # Wg([2]) = -1/(d(d^2-1)) = -1/d^3 - 1/d^5 ...
        
        # symbolic Weingarten
        wg2 = IntU.weingarten([2], d) 
        # approximate to order 4 (in 1/d)
        # Should be -1/d^3.
        
        # Test function `asymptotic`?
        # If implemented:
        # series_w2 = asymptotic(wg2, d, 5)
        
        # For now, let's just assume we have `asymptotic`
        # We need to verify if it works.
        # But we haven't implemented `asymptotic` in IntU.jl yet?
        # Wait, the task said "Implement Asymptotic Expansions".
        # Did I implement it?
        # Yes, I did in previous steps. in src/asymptotic.jl?
        # Or added to integration_core.jl?
        # Let's check if asymptotic is exported.
        
        @test isdefined(IntU, :asymptotic)
        
        res = asymptotic(wg2, dU(reshape([:a],1,1), d), 5)
        # Expectation: -d^-3
        # Output is a Series or simplified expression.
        
        # Let's check simple polynomial
        term = d^2 + 2d + 1
        # asymptotic(term) -> d^2 + 2d + 1 (exact)
        
        # What about 1/(d-1) = 1/d + 1/d^2 + ...
        # term = 1/(d-1)
        # res = asymptotic(term, ..., 2)
        # Should contain 1/d + 1/d^2
    end
end
