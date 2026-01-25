using IntU
using Test
using Symbolics

# Helper to convert symbolic results to numbers
function to_numeric(x)
    x_un = Symbolics.unwrap(x)
    if x_un isa Number
        return x_un
    end
    # If it's a symbolic constant (like Num(0)), try to simplify/substitute
    sim = Symbolics.simplify(x)
    sim_un = Symbolics.unwrap(sim)
    if sim_un isa Number
        return sim_un
    end
    
    # Try brute force substitution
    sim2 = Symbolics.substitute(sim, Dict())
    sim2_un = Symbolics.unwrap(sim2)
    if sim2_un isa Number
        return sim2_un
    end
    
    if isequal(sim, 0) || isequal(sim2, 0)
        return 0.0
    end
    
    # Last resort: try evaluating
    try
        val = eval(Meta.parse(string(sim)))
        if val isa Number
            return val
        end
    catch
    end
    
    if isequal(sim, 0)
        return 0.0
    end
    return x_un
end

@testset "IntU.jl Suite" begin
    @testset verbose=true "Weingarten Calculus" begin
        include("weingarten.jl")
    end

    @testset verbose=true "Haar Measure" begin
        include("haar_measure.jl")
    end

    @testset verbose=true "Pure States" begin
        include("pure_states.jl")
    end

    @testset verbose=true "QI Helpers" begin
        include("qi_helpers.jl")
    end
end
