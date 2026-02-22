# Real and Symplectic measures

# Dummy types to represent the measures
struct OrthogonalMeasure{D} <: AbstractMeasure
    dim::D
end

struct SymplecticMeasure{D} <: AbstractMeasure
    dim::D
end

@doc raw"""
    dO(dim)

Defines the Haar measure for the real Orthogonal group $O(d)$ with dimension `dim`.
Integration engine identifies variables via metadata tag `:O`.
"""
dO(dim) = OrthogonalMeasure(dim)

@doc raw"""
    dSp(dim)

Defines the Haar measure for the Symplectic group $Sp(d)$. 
The dimension `dim` must be even.
Integration engine identifies variables via metadata tag `:Sp`.
"""
dSp(dim) = SymplecticMeasure(dim)


function IntU.measure_info(measure::OrthogonalMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = MetadataMatcher(:O)
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, :O)
end

function _j_pair_sign(idx, n)
    if idx <= n
        return idx + n, 1
    else
        return idx - n, -1
    end
end

function IntU.measure_info(measure::SymplecticMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = MetadataMatcher(:Sp)
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, :Sp)
end

"""
    asymptotic(expr, measure::Union{OrthogonalMeasure, SymplecticMeasure}, order=1)
"""
function asymptotic(expr, measure::Union{OrthogonalMeasure,SymplecticMeasure}, order = 1)
    d = measure.dim
    # If d is symbolic or not an integer, we can proceed directly
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    # If d is integer, checking asymptotic might require symbolic d
    d_asymp = Symbolics.variable(:d_asymp)
    # Reconstruct measure with symbolic dim
    m_sym = if measure isa OrthogonalMeasure
        dO(d_asymp)
    else
        dSp(d_asymp)
    end

    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end
