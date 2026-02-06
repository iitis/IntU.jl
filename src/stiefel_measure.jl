# Stiefel Manifold integration

"""
    dStiefel(V, dim, k)

Defines the measure for integration over the Stiefel manifold V_k(C^d).
This manifold represents the set of d x k matrices with orthonormal columns.

The integration is performed by mapping V to the first k columns of a
Haar-random unitary matrix U(d).

Reference:
- Edelman, A., Arias, T. A., & Smith, S. T. (1998). The geometry of algorithms with orthogonality constraints.
"""
dStiefel(V, dim, k) = StiefelMeasure(V, dim, k)


"""
    StiefelMeasure(V, dim, k)

Internal type representing the measure on the Stiefel manifold. 
Users should use `dStiefel` constructor.
"""
struct StiefelMeasure{M,D,K}
    V::M
    dim::D
    k::K
end

"""
    integrate(expr, measure::StiefelMeasure)
"""
function integrate(expr::AbstractArray, measure::StiefelMeasure)
    return map(e -> integrate(e, measure), expr)
end

function fallback_integrate(expr, measure::StiefelMeasure)
    V_input = measure.V
    dim = measure.dim
    k_dim = measure.k

    subs_dict = Dict{Any,Any}()
    V_atomic_lookup = Dict{Any,Tuple{Int,Int}}()
    V_bar_lookup = Dict{Any,Tuple{Int,Int}}()

    # V is a d x k matrix
    if V_input isa AbstractArray
        rows = size(V_input, 1)
        cols = size(V_input, 2)
        
        # We allow V to be symbolic or explicit, but loops logic assumes we traverse it
        for i = 1:rows
            for j = 1:cols
                v_ij_num = _safe_Num(V_input[i, j])
                v_ij_un = Symbolics.unwrap(v_ij_num)
                
                # Create atomic variables for mapping
                # Map V_{i,j} to U_{i,j}
                # However, since j is in range 1:k, this maps to the first k columns of U
                
                # Check bounds if k is explicitly given as integer
                if k_dim isa Integer && j > k_dim
                     error("Matrix column index $j exceeds Stiefel dimension k=$k_dim")
                end

                v_atomic = Symbolics.variable(:V_atomic, i, j)
                v_bar_atomic = Symbolics.variable(:V_bar_atomic, i, j)

                # Store mapping to U indices (i, j)
                V_atomic_lookup[Symbolics.unwrap(v_atomic)] = (i, j)
                V_bar_lookup[Symbolics.unwrap(v_bar_atomic)] = (i, j)

                subs_dict[v_ij_un] = v_atomic

                c_ij_un = Symbolics.unwrap(conj(v_ij_num))
                subs_dict[c_ij_un] = v_bar_atomic
                
                bc_ij_un = Symbolics.unwrap(Base.conj(v_ij_num))
                subs_dict[bc_ij_un] = v_bar_atomic
            end
        end
    end

    matcher = LookupMatcher(V_atomic_lookup, V_bar_lookup)
    return _robust_real_num(_integrate_core(expr, dim, subs_dict, matcher))
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
    m_sym = dStiefel(measure.V, d_asymp, measure.k)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end
