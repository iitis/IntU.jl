using IntU
using Symbolics
using LinearAlgebra
using Test

# Helper to convert symbolic results to numbers - similar to runtests.jl
function to_numeric(x)
    x_un = Symbolics.unwrap(x)
    if x_un isa Number
        return x_un
    end
    # If it's a symbolic constant, try to simplify
    sim = simplify(x)
    sim_un = Symbolics.unwrap(sim)
    if sim_un isa Number
        return sim_un
    end
    
    # Try evaluating
    try
        val = eval(Symbolics.toexpr(sim_un))
        if val isa Number
            return val
        end
    catch
    end
    
    return x_un
end

@testset "Circular Ensembles" begin
    @variables d
    
    # helper for symbolic zero check
    function sym_iszero(x)
        # 1. Quick check
        if isequal(x, 0) || isequal(x, 0.0)
            return true
        end
        
        # 2. Simplify and check
        s = simplify(x)
        if isequal(s, 0) || isequal(s, 0.0)
            return true
        end
        un = Symbolics.unwrap(s)
        if isequal(Num(un), 0)
            return true
        end
        
        # 3. Numerical fallback
        try
            for d_val in [3.14, 1.23, 7.89]
                val = Symbolics.substitute(s, Dict(d => d_val))
                num_val = to_numeric(val)
                if num_val isa Number
                    if abs(num_val) > 1e-12
                        return false
                    end
                else
                    return false
                end
            end
            return true
        catch
            return false
        end
    end

    # 1. COE
    @testset "COE Properties" begin
        @variables S[1:2, 1:2]::Complex
        m = dCOE(S, d)
        
        @test sym_iszero(integrate(S[1,1], m))
        @test sym_iszero(integrate(S[1,2], m))
        
        # E[S_11 S*_11] = 2/(d+1)
        val_11 = integrate(S[1,1] * conj(S[1,1]), m)
        @test sym_iszero(val_11 - 2 / (d + 1))
        
        # E[S_12 S*_12] = 1/(d+1)
        val_12 = integrate(S[1,2] * conj(S[1,2]), m)
        @test sym_iszero(val_12 - 1 / (d + 1))
        
        # Symmetry check
        val_sym = integrate(S[1,2] * conj(S[2,1]), m)
        @test sym_iszero(val_sym - val_12)
    end
    
    # 2. CSE
    @testset "CSE Properties" begin
        @variables S_cse[1:2, 1:2]::Complex
        mc = dCSE(S_cse, d) 
        
        @test sym_iszero(integrate(S_cse[1,1], mc))
        
        # E[|S_11|^2] = 1/(d-1)
        val_11 = integrate(S_cse[1,1] * conj(S_cse[1,1]), mc)
        @test sym_iszero(val_11 - 1 / (d - 1))
        
        val_12 = integrate(S_cse[1,2] * conj(S_cse[1,2]), mc)
        @test sym_iszero(val_12)
        
        val_cross = integrate(S_cse[1,2] * conj(S_cse[2,1]), mc)
        @test sym_iszero(val_cross + val_12)
    end
    
    # 3. CUE
    @testset "CUE Properties" begin
        @variables U[1:2, 1:2]::Complex
        m = dCUE(U, d)
        
        @test sym_iszero(integrate(U[1,1], m))
        val_11 = integrate(U[1,1] * conj(U[1,1]), m)
        @test sym_iszero(val_11 - 1/d)
        
        val_22 = integrate(abs(U[1,1]*U[2,2])^2, m)
        @test sym_iszero(val_22 - 1/(d^2-1))
    end

end
