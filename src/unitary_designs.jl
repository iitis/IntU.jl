# Unitary t-designs

struct UnitaryDesign{M,D,T}
    U::M
    dim::D
    t::T
end

"""
    dDesign(U, dim, t)

Defines a measure representing a **unitary t-design**. 

A t-design is a set of unitaries that reproduces the first t moments 
of the Haar measure. Integration of any polynomial P(U, \\bar{U}) of 
degree (q, q) with q \\le t yields the same result as the Haar measure.

Reference:
- Gross, D., Audenaert, K., & Eisert, J. (2007). Evenly distributed unitaries: On the structure of unitary designs.
"""
dDesign(U, dim, t) = UnitaryDesign(U, dim, t)

"""
    integrate(expr, measure::UnitaryDesign)
"""
function integrate(expr::AbstractArray, measure::UnitaryDesign)
    return map(e -> integrate(e, measure), expr)
end

function IntU.measure_info(measure::UnitaryDesign)
    U_sym = measure.U
    dim = measure.dim
    t_val = measure.t

    if U_sym isa SymbolicMatrix
        subs_dict = Dict{Any,Any}()
        matcher = SymbolicMatcher(:U, Regex("^$(U_sym.name)_(\\d+)_(\\d+)\$"))
        return (subs_dict, matcher, dim, (:Design, t_val))
    end

    # Substitute Re(U) and Im(U) - reusing logic similar to HaarMeasure
    subs_dict = Dict{Any,Any}()
    U_atomic_lookup = Dict{Any,Tuple}()
    U_bar_lookup = Dict{Any,Tuple}()

    if U_sym isa AbstractArray
        for i = 1:size(U_sym, 1)
            for j = 1:size(U_sym, 2)
                u_ij_num = _safe_Num(U_sym[i, j])
                u_ij_un = Symbolics.unwrap(u_ij_num)
                u_atomic = Symbolics.variable(:U_atomic, i, j)
                u_bar_atomic = Symbolics.variable(:U_bar_atomic, i, j)

                U_atomic_lookup[Symbolics.unwrap(u_atomic)] = (i, j)
                U_bar_lookup[Symbolics.unwrap(u_bar_atomic)] = (i, j)

                subs_dict[u_ij_un] = u_atomic
                subs_dict[Symbolics.unwrap(conj(u_ij_num))] = u_bar_atomic
                subs_dict[Symbolics.unwrap(Base.conj(u_ij_num))] = u_bar_atomic
            end
        end
    end

    matcher = LookupMatcher(U_atomic_lookup, U_bar_lookup)
    return (subs_dict, matcher, dim, (:Design, t_val))
end
