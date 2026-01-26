# Gaussian Random Matrix measures (GUE, GOE, GSE)

struct GUEMeasure{T, N, D}
    H::AbstractArray{T, N}
    dim::D
end

"""
    dGUE(H, dim)

Define the measure for the Gaussian Unitary Ensemble (GUE).
`H` is the symbolic matrix representing the Hermitian Gaussian random matrix.
`dim` is the dimension (symbolic or integer).

Expectation values are defined by Wick's theorem with the contraction:
`< H_{ij} H_{kl} > = delta_{il} * delta_{jk}`

This normalization corresponds to `< Tr(H^2) > = dim^2`.
Note: H is Hermitian, so `conj(H_{ij})` is treated as `H_{ji}`.
"""
dGUE(H::AbstractArray{T,N}, dim) where {T,N} = GUEMeasure{T,N,typeof(dim)}(H, dim)

"""
    integrate(expr, measure::GUEMeasure)
"""
function integrate(expr::AbstractArray, measure::GUEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr, measure::GUEMeasure)
    H_sym = measure.H
    dim = measure.dim
    
    subs_dict = Dict{Any, Any}()
    H_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    
    if H_sym isa AbstractArray
        for i in 1:size(H_sym, 1)
            for j in 1:size(H_sym, 2)
                h_ij_num = _safe_Num(H_sym[i,j])
                h_ij_un = Symbolics.unwrap(h_ij_num)
                h_atomic = Symbolics.variable(:H_atomic, i, j)
                
                H_atomic_lookup[Symbolics.unwrap(h_atomic)] = (i, j)
                
                subs_dict[h_ij_un] = h_atomic
                
                # Handle conjugates
                # conj(H_{ij}) = H_{ji}
                # We map conj(H) to H_bar_atomic which we will identify later as H_{ji}
                # Actually, simpler: map conj(H_{ij}) to H_atomic(j, i)
                # But we can't easily create a variable "pointing" to another variable's indices
                # implicitly without separate lookup.
                
                # So we map conj(h_{ij}) to a new atomic variable h_bar_{ij}
                # and in the lookup we store (j, i) for it!
                
                hb_atomic = Symbolics.variable(:H_bar_atomic, i, j)
                # Note: For H_{ij}, the conjugate is H_{ji}.
                # So if we see conj(H_{ij}), we treat it as H at indices (j, i).
                H_atomic_lookup[Symbolics.unwrap(hb_atomic)] = (j, i)
                
                subs_dict[Symbolics.unwrap(conj(h_ij_un))] = hb_atomic
                subs_dict[Symbolics.unwrap(Base.conj(h_ij_un))] = hb_atomic
            end
        end
    end

    # For GUE, we don't separate U and U_bar. All are H.
    # We pass empty U_bar_lookup because we mapped everything to H_atomic_lookup (or compatible)
    return _integrate_core(expr, dim, subs_dict, H_atomic_lookup, Dict(), :GUE)
end

"""
    asymptotic(expr, measure::GUEMeasure, order=1)
"""
function asymptotic(expr, measure::GUEMeasure, order=1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end
    
    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dGUE(measure.H, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end
