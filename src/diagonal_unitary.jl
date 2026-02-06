# Diagonal Unitary Matrices (Torus group) Integration

@doc raw"""
    DiagonalUnitaryMeasure{M, D}

Represents integration over the group of diagonal unitary matrices (the torus \$T^d\$).
For a diagonal unitary matrix \$V\$, the only non-zero entries are \$V_{ii} = e^{i\theta_i}\$.
Integration over \$T^d\$ is equivalent to independent phase integrations for each diagonal entry.
"""
struct DiagonalUnitaryMeasure{M,D}
    V::M
    dim::D
end

"""
    dDiagUnitary(V, dim)
    dDiagUnitary(dim)

Defines the measure for the group of diagonal unitary matrices of dimension `dim`.
"""
dDiagUnitary(V, dim) = DiagonalUnitaryMeasure(V, dim)
dDiagUnitary(dim) = DiagonalUnitaryMeasure(nothing, dim)

"""
    integrate(expr, measure::DiagonalUnitaryMeasure)
"""
function integrate(expr::AbstractArray, measure::DiagonalUnitaryMeasure)
    return map(e -> integrate(e, measure), expr)
end

function fallback_integrate(expr, measure::DiagonalUnitaryMeasure)
    V_sym = measure.V
    dim = measure.dim

    subs_dict = Dict{Any,Any}()
    V_atomic_lookup = Dict{Any,Tuple{Int,Int}}()
    V_bar_lookup = Dict{Any,Tuple{Int,Int}}()

    if V_sym isa AbstractArray
        for i = 1:size(V_sym, 1)
            for j = 1:size(V_sym, 2)
                v_ij_num = _safe_Num(V_sym[i, j])
                v_ij_un = Symbolics.unwrap(v_ij_num)
                
                v_atomic = Symbolics.variable(:V_atomic, i, j)
                v_bar_atomic = Symbolics.variable(:V_bar_atomic, i, j)

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
    return _robust_real_num(_integrate_core(expr, dim, subs_dict, matcher, :DiagUnitary))
end
