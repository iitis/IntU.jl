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

    if IntU._symbolic_isequal(sim, 0) || IntU._symbolic_isequal(sim2, 0)
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

    if IntU._symbolic_isequal(sim, 0)
        return 0.0
    end
    return x_un
end

@testset verbose=true "IntU.jl Suite" begin
    @testset verbose=true "Aqua Tests" begin
        include("aqua.jl")
    end

    @testset verbose=true "Weingarten Calculus" begin
        include("weingarten.jl")
    end

    @testset verbose=true "Haar Measure" begin
        include("haar_measure.jl")
    end

    @testset verbose=true "Pure States" begin
        include("pure_states.jl")
    end

    @testset verbose=true "Orthogonal Group" begin
        include("orthogonal.jl")
    end

    @testset verbose=true "QI Helpers" begin
        include("qi_helpers.jl")
    end

    @testset verbose=true "Asymptotic Expansions" begin
        include("asymptotic.jl")
    end

    @testset verbose=true "Symbolic Trace" begin
        include("symbolic_trace.jl")
    end

    @testset verbose=true "GUE Integration" begin
        include("gue.jl")
    end

    @testset verbose=true "GOE Integration" begin
        include("goe.jl")
    end

    @testset verbose=true "Gaussian Miscellaneous" begin
        include("gaussian_misc.jl")
    end

    @testset verbose=true "Integral Library" begin
        include("library.jl")
    end

    @testset verbose=true "Circular Ensembles" begin
        include("test_circular.jl")
    end

    @testset verbose=true "ITensors Bridge" begin
        include("itensors_bridge.jl")
    end

    @testset verbose=true "Unitary Designs" begin
        include("unitary_designs.jl")
    end
    
    @testset verbose=true "Permutation Groups" begin
        include("permutation_groups.jl")
    end

    @testset verbose=true "SU(d) Integration" begin
        include("su_tests.jl")
    end
end
