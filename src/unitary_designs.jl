# Unitary t-designs

struct UnitaryDesign{M, D, T}
    U::M
    dim::D
    t::T
end

"""
    dDesign(U, dim, t)

Define a unitary t-design measure.
`U` is the symbolic matrix, `dim` is the dimension, and `t` is the design order.
Integrals of polynomials in elements of `U` and `conj(U)` will match Haar measure results
if the degree in `U` (and `conj(U)`) is at most `t`.
If the degree exceeds `t`, `integrate` will throw an error.
"""
dDesign(U, dim, t) = UnitaryDesign(U, dim, t)

"""
    integrate(expr, measure::UnitaryDesign)
"""
function integrate(expr::AbstractArray, measure::UnitaryDesign)
    return map(e -> integrate(e, measure), expr)
end

function fallback_integrate(expr, measure::UnitaryDesign)
    U_sym = measure.U
    dim = measure.dim
    t_val = measure.t
    
    # Substitute Re(U) and Im(U) - reusing logic similar to HaarMeasure
    subs_dict = Dict{Any, Any}()
    U_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    U_bar_lookup = Dict{Any, Tuple{Int, Int}}()
    
    if U_sym isa AbstractArray
        for i in 1:size(U_sym, 1)
            for j in 1:size(U_sym, 2)
                u_ij_num = _safe_Num(U_sym[i,j])
                u_ij_un = Symbolics.unwrap(u_ij_num)
                u_atomic = Symbolics.variable(:U_atomic, i, j)
                u_bar_atomic = Symbolics.variable(:U_bar_atomic, i, j)
                
                U_atomic_lookup[Symbolics.unwrap(u_atomic)] = (i, j)
                U_bar_lookup[Symbolics.unwrap(u_bar_atomic)] = (i, j)
                
                subs_dict[u_ij_un] = u_atomic
                
                c_ij_un = Symbolics.unwrap(conj(u_ij_num))
                subs_dict[c_ij_un] = u_bar_atomic
                
                bc_ij_un = Symbolics.unwrap(Base.conj(u_ij_num))
                subs_dict[bc_ij_un] = u_bar_atomic
            end
        end
    end

    # Pass measure_type as a tuple (:Design, t)
    matcher = LookupMatcher(U_atomic_lookup, U_bar_lookup)
    return _robust_real_num(_integrate_core(expr, dim, subs_dict, matcher, (:Design, t_val)))
end
