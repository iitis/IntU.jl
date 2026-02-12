# Stiefel Manifold integration

@doc raw"""
    dStiefel(dim, k)
    dStiefel(V::SymbolicMatrix, k)

Defines the measure for integration over the Stiefel manifold $V_k(\mathbb{C}^d)$.
This manifold represents the set of $d \times k$ matrices with orthonormal columns.

The integration is performed by mapping $V$ to the first $k$ columns of a Haar-random unitary matrix $U(d)$.
If called with `dim`, it integrates entries tagged with `:U` via `SymbolicMatrix(:V, :U)`.

Reference:
- Edelman, A., Arias, T. A., & Smith, S. T. (1998). The geometry of algorithms with orthogonality constraints.
"""
dStiefel(dim, k) = StiefelMeasure(dim, k)
dStiefel(V::SymbolicMatrix, k) = StiefelMeasure(V.dim, k)


"""
    StiefelMeasure(dim, k)

Internal type representing the measure on the Stiefel manifold. 
Users should use `dStiefel` constructor.
"""
struct StiefelMeasure{D,K}
    dim::D
    k::K
end

"""
    integrate(expr, measure::StiefelMeasure)
"""
function integrate(expr::AbstractArray, measure::StiefelMeasure)
    return map(e -> integrate(e, measure), expr)
end

function IntU.measure_info(measure::StiefelMeasure)
    subs_dict = Dict{Any,Any}()
    # Use metadata-based matching for U. 
    # Stiefel is just first k columns of U, handled during entry matching if needed, 
    # but for symbolic matching we just tag as U.
    matcher = MetadataMatcher(:U)
    return (subs_dict, matcher, measure.dim, :U)
end

"""
    asymptotic(expr, measure::StiefelMeasure, order=1)

Returns the series expansion of the integral in powers of `1/d`.
"""
function asymptotic(expr, measure::StiefelMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dStiefel(d_asymp, measure.k)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end
