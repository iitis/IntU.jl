using IntU
using Test
using Symbolics

# Helper to convert symbolic results to numbers
function to_numeric(x)
    x_un = Symbolics.unwrap(x)
    if x_un isa Number
        return x_un
    end
    sim = Symbolics.simplify(x)
    sim_un = Symbolics.unwrap(sim)
    if sim_un isa Number
        return sim_un
    end
    if IntU._symbolic_isequal(sim, 0)
        return 0.0
    end
    # Simple substitution
    sim2 = Symbolics.substitute(sim, Dict())
    sim2_un = Symbolics.unwrap(sim2)
    if sim2_un isa Number; return sim2_un; end
    
    # Last resort: try evaluating
    try
        val = eval(Meta.parse(string(sim)))
        if val isa Number
            return val
        end
    catch
    end
    
    return x_un
end



@testset verbose=true "Unitary t-Designs" begin
    d_val = 3
    @variables U[1:3, 1:3]::Complex
    
    # Define a 2-design
    design2 = dDesign(U, d_val, 2)
    
    @testset verbose=true "Degree <= t (t=2)" begin
        # Degree 1 integrals (should match Haar)
        expr1 = abs(U[1,1])^2
        res1 = integrate(expr1, design2)
        @test to_numeric(real(res1)) ≈ 1/3
        
        # Degree 2 integrals (should match Haar)
        expr2 = abs(U[1,1] * U[2,2])^2
        res2 = integrate(expr2, design2)
        @test to_numeric(real(res2)) ≈ 1/8
        
        # Unitarity Check (Degree 1 in U, 1 in U_bar -> Total degree 1 <= 2)
        sum_val = 0//1
        for k in 1:d_val
            sum_val += integrate(U[1,k] * conj(U[1,k]), design2)
        end
        @test to_numeric(real(sum_val)) ≈ 1.0
    end
    
    @testset verbose=true "Degree > t (t=2, degree=3)" begin
        # Degree 3 integral (should fail)
        expr3 = abs(U[1,1])^6 # |u|^6 -> u^3 * conj(u)^3 -> degree 3
        @test_throws ErrorException integrate(expr3, design2)
    end
    
    # Define a 1-design
    design1 = dDesign(U, d_val, 1)
    
    @testset verbose=true "1-Design Constraints" begin
        # Degree 1 works
        expr1 = abs(U[1,1])^2
        res1 = integrate(expr1, design1)
        @test to_numeric(real(res1)) ≈ 1/3
        
        # Degree 2 fails
        expr2 = abs(U[1,1])^4 # degree 2
        @test_throws ErrorException integrate(expr2, design1)
    end
end
