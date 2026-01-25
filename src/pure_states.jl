# Pure state integration

# Type to represent integration over pure states |psi>
struct PureStateMeasure{T, N, D}
    psi::AbstractArray{T, N}
    dim::D
end
"""
    dPsi(psi, dim)

Define the Haar measure for random pure states (vectors) in dimension `dim`.
`psi` is the symbolic vector representing the state, and `dim` is the dimension (symbolic or integer).
"""
dPsi(psi::AbstractArray{T,N}, dim) where {T,N} = PureStateMeasure{T,N,typeof(dim)}(psi, dim)

"""
    integrate(expr, measure::PureStateMeasure)
"""
function integrate(expr::AbstractArray, measure::PureStateMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr, measure::PureStateMeasure)
    psi = measure.psi
    dim = measure.dim
    n = length(psi)
    
    
    subs_dict = Dict()
    psi_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    psi_bar_lookup = Dict{Any, Tuple{Int, Int}}()
    
    psi_vec = collect(psi)
    for i in 1:n
        u_elem = psi_vec[i]
        u_atomic = Symbolics.variable(:psi_atomic, i)
        u_bar_atomic = Symbolics.variable(:psi_bar_atomic, i)
        
        psi_atomic_lookup[Symbolics.unwrap(u_atomic)] = (i, 1)
        psi_bar_lookup[Symbolics.unwrap(u_bar_atomic)] = (i, 1)
        
        subs_dict[u_elem] = u_atomic
        subs_dict[conj(u_elem)] = u_bar_atomic
        subs_dict[real(u_elem)] = (1//2) * (u_atomic + u_bar_atomic)
        subs_dict[imag(u_elem)] = (1//(2im)) * (u_atomic - u_bar_atomic)
    end
    
    return _integrate_core(expr, dim, subs_dict, psi_atomic_lookup, psi_bar_lookup)
end

"""
    asymptotic(expr, measure::PureStateMeasure, order=1)

Returns the series expansion of the integral in powers of `1/d`.
"""
function asymptotic(expr, measure::PureStateMeasure, order=1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end
    
    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dPsi(measure.psi, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end
