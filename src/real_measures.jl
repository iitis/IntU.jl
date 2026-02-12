# Real and Symplectic measures

# Dummy types to represent the measures
struct OrthogonalMeasure{D}
    dim::D
end

struct SymplecticMeasure{D}
    dim::D
end

@doc raw"""
    dO(dim)
    dO(O::SymbolicMatrix)

Defines the Haar measure for the real Orthogonal group $O(d)$ with dimension `dim`.

If called with `dim`, it integrates entries tagged with `:O` via `SymbolicMatrix(:O, :O)`.

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
dO(dim) = OrthogonalMeasure(dim)
dO(O::SymbolicMatrix) = OrthogonalMeasure(O.dim)

@doc raw"""
    dSp(dim)
    dSp(S::SymbolicMatrix)

Defines the Haar measure for the Symplectic group $Sp(d)$. 
The dimension `dim` must be even.

If called with `dim`, it integrates entries tagged with `:Sp` via `SymbolicMatrix(:S, :Sp)`.

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
dSp(dim) = SymplecticMeasure(dim)
dSp(S::SymbolicMatrix) = SymplecticMeasure(S.dim)


"""
    integrate(expr, measure::OrthogonalMeasure)
"""
function integrate(expr::AbstractArray, measure::OrthogonalMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::SymplecticMeasure)
    return map(e -> integrate(e, measure), expr)
end

# Resolve ambiguity with SymbolicMatrixProduct
function integrate(P::SymbolicMatrixProduct, measure::OrthogonalMeasure)
    # create a method instance for the generic SymbolicMatrixProduct integration
    return invoke(integrate, Tuple{SymbolicMatrixProduct, Any}, P, measure)
end

function integrate(P::SymbolicMatrixProduct, measure::SymplecticMeasure)
    return invoke(integrate, Tuple{SymbolicMatrixProduct, Any}, P, measure)
end

function IntU.measure_info(measure::OrthogonalMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = MetadataMatcher(:O)
    return (subs_dict, matcher, measure.dim, :O)
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
    return (subs_dict, matcher, measure.dim, :Sp)
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
