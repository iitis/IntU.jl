# Haar measure integration

# Dummy type to represent the measure
struct HaarMeasure{T, N, D}
    U::AbstractArray{T, N}
    dim::D
end
dU(U::AbstractArray{T,N}, dim) where {T,N} = HaarMeasure{T,N,typeof(dim)}(U, dim)

"""
    integrate(expr, measure::HaarMeasure)
"""
function integrate(expr::AbstractArray, measure::HaarMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr, measure::HaarMeasure)
    U_sym = measure.U
    dim = measure.dim
    
    # Substitute Re(U) and Im(U)
    subs_dict = Dict{Any, Any}()
    U_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    U_bar_lookup = Dict{Any, Tuple{Int, Int}}()
    

    if U_sym isa AbstractArray
        for i in 1:size(U_sym, 1)
            for j in 1:size(U_sym, 2)
                u_ij = U_sym[i,j]
                u_atomic = Symbolics.variable(:U_atomic, i, j)
                u_bar_atomic = Symbolics.variable(:U_bar_atomic, i, j)
                
                U_atomic_lookup[Symbolics.unwrap(u_atomic)] = (i, j)
                U_bar_lookup[Symbolics.unwrap(u_bar_atomic)] = (i, j)
                
                subs_dict[u_ij] = u_atomic
                subs_dict[conj(u_ij)] = u_bar_atomic
                subs_dict[real(u_ij)] = (1//2) * (u_atomic + u_bar_atomic)
                subs_dict[imag(u_ij)] = (1//(2im)) * (u_atomic - u_bar_atomic)
            end
        end
    end

    return _integrate_core(expr, dim, subs_dict, U_atomic_lookup, U_bar_lookup)
end
