using IntU
using Test
using Symbolics

function to_numeric(x)
    ux = Symbolics.unwrap(x)
    xf = SymbolicUtils.Postwalk(
        t -> begin
            ut = Symbolics.unwrap(t)
            if Symbolics.iscall(ut) && (
                Symbolics.operation(ut) == complex ||
                Symbolics.operation(ut) == Base.complex
            )
                aa = Symbolics.arguments(ut)
                return aa[1] + im*aa[2]
            end
            return t
        end,
    )(
        ux,
    )
    sim = Symbolics.simplify(Symbolics.wrap(xf))
    v = Symbolics.value(sim)
    if v isa Number
        return v
    end
    if v isa Real
        return Float64(v)
    elseif v isa Number
        return ComplexF64(v)
    end
    return v
end

function is_really_zero(x)
    x = Symbolics.simplify(x)
    IntU._symbolic_isequal(x, 0) && return true
    x = Symbolics.expand(x)
    x = Symbolics.simplify(x)
    IntU._symbolic_isequal(x, 0) && return true

    vars = Symbolics.get_variables(x)
    if isempty(vars)
        v = Symbolics.value(x)
        return v isa Number && abs(v) < 1e-10
    end

    for i = 1:3
        subs = Dict(v => 0.1 + 0.3 * i + 0.17 * j for (j, v) in enumerate(vars))
        val_sub = Symbolics.substitute(x, subs)
        v = to_numeric(val_sub)
        if v isa Number
            if abs(v) < 1e-9
                continue
            end
            return false
        end
        return false
    end
    return true
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
        include("bounds_checking.jl")
        include("kron_integration.jl")
        include("kron_enhancements.jl")
    end

    @testset verbose=true "GUE Integration" begin
        include("gue.jl")
    end

    @testset verbose=true "GOE Integration" begin
        include("goe.jl")
    end

    @testset verbose=true "Gaussian Miscellaneous" begin
        include("gaussian_misc.jl")
        include("gaussian_symbolic_trace.jl")
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

    @testset verbose=true "Diagonal Unitary Integration" begin
        include("test_diagonal_unitary.jl")
    end

    @testset verbose=true "HCIZ Integration" begin
        include("hciz.jl")
    end

    @testset verbose=true "Stiefel Manifold" begin
        include("stiefel_tests.jl")
    end

    @testset verbose=true "Ginibre Ensembles" begin
        include("ginibre.jl")
    end

    @testset verbose=true "Float Dimension Rejection" begin
        include("float_dimension_rejection.jl")
    end

    @testset verbose=true "UX Improvements" begin
        include("ux_improvements.jl")
        include("macro_regression.jl")
        include("lazy_abs.jl")
    end
end
