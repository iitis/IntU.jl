# Permutation Group measures

# Dummy types to represent the measures
struct PermutationMeasure{M,D}
    P::M
    dim::D
end

struct CenteredPermutationMeasure{M,D}
    Y::M
    dim::D
end

"""
    dPerm(P, dim)

Defines the Haar measure for the Symmetric group S_d (permutation matrices).

The integration of a monomial of entries is given by:
```math
\\int_{S_d} P_{i_1 j_1} \\dots P_{i_n j_n} dP = \\begin{cases} \\frac{(d-k)!}{d!} & \\text{if indices are consistent} \\\\ 0 & \\text{otherwise} \\end{cases}
```
where k is the number of distinct pairs (i, j) in the product.
"""
dPerm(P, dim) = PermutationMeasure(P, dim)
dPerm(dim) = PermutationMeasure(nothing, dim)
dPerm(P::LazySymbolicMatrix) = PermutationMeasure(P, P.dim)

"""
    dCPerm(Y, dim)

Defines the measure for Centered Permutation matrices Y = P - J/d.
"""
dCPerm(Y, dim) = CenteredPermutationMeasure(Y, dim)
dCPerm(dim) = CenteredPermutationMeasure(nothing, dim)
dCPerm(Y::LazySymbolicMatrix) = CenteredPermutationMeasure(Y, Y.dim)


function integrate(expr::AbstractArray, measure::PermutationMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::CenteredPermutationMeasure)
    return map(e -> integrate(e, measure), expr)
end

function IntU.measure_info(measure::PermutationMeasure)
    P_sym = measure.P
    dim = measure.dim

    if P_sym isa LazySymbolicMatrix
        error("Symbolic dimension support is not available for the permutation group.")
    end

    subs_dict = Dict{Any,Any}()
    P_atomic_lookup = Dict{Any,Tuple}()

    if P_sym isa AbstractArray
        for i = 1:size(P_sym, 1)
            for j = 1:size(P_sym, 2)
                p_ij_num = _safe_Num(P_sym[i, j])
                p_ij_un = Symbolics.unwrap(p_ij_num)
                p_atomic = Symbolics.variable(:P_atomic, i, j)

                P_atomic_lookup[Symbolics.unwrap(p_atomic)] = (i, j)

                subs_dict[p_ij_un] = p_atomic
                # P is real
                subs_dict[Symbolics.unwrap(conj(p_ij_un))] = p_atomic
                subs_dict[Symbolics.unwrap(Base.conj(p_ij_un))] = p_atomic
            end
        end
    end

    matcher = LookupMatcher(P_atomic_lookup, Dict{Any,Tuple}())
    return (subs_dict, matcher, dim, :Perm)
end

function IntU.measure_info(measure::CenteredPermutationMeasure)
    Y_sym = measure.Y
    dim = measure.dim
    
    subs_dict = Dict{Any,Any}()
    P_atomic_lookup = Dict{Any,Tuple}()

    if Y_sym isa AbstractArray
        for i = 1:size(Y_sym, 1)
            for j = 1:size(Y_sym, 2)
                y_ij_num = _safe_Num(Y_sym[i, j])
                y_ij_un = Symbolics.unwrap(y_ij_num)

                # Y_ij = P_ij - 1/dim
                p_atomic = Symbolics.variable(:P_atomic, i, j)
                P_atomic_lookup[Symbolics.unwrap(p_atomic)] = (i, j)

                subs_dict[y_ij_un] = p_atomic - 1/dim
                subs_dict[Symbolics.unwrap(conj(y_ij_un))] = p_atomic - 1/dim
                subs_dict[Symbolics.unwrap(Base.conj(y_ij_un))] = p_atomic - 1/dim
            end
        end
    end

    matcher = LookupMatcher(P_atomic_lookup, Dict{Any,Tuple}())
    return (subs_dict, matcher, dim, :Perm)
end

"""
    integrate_indices_permutation(indices, dim)

Integration over the symmetric group S_d.
"""
function integrate_indices_permutation(indices::Vector{Tuple{Int,Int}}, dim)
    if isempty(indices)
        return 1
    end

    unique_pairs = unique(indices)
    k = length(unique_pairs)

    # Check consistency
    rows = [p[1] for p in unique_pairs]
    cols = [p[2] for p in unique_pairs]

    if length(unique(rows)) != k || length(unique(cols)) != k
        return 0
    end

    # Result is 1 / (d * (d-1) * ... * (d-k+1))
    # which is (d-k)! / d!

    # Handle symbolic or large numeric dim
    res = (dim isa Integer) ? BigInt(1) // BigInt(1) : 1 // 1
    for m = 0:(k-1)
        res /= (dim - m)
    end

    return res
end
