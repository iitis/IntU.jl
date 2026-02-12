# Permutation Group measures

# Dummy types to represent the measures
struct PermutationMeasure{D}
    dim::D
end

struct CenteredPermutationMeasure{D}
    dim::D
end

@doc raw"""
    dPerm(dim)
    dPerm(P::SymbolicMatrix)

Defines the Haar measure for the Symmetric group $S_d$ (permutation matrices) of dimension `dim`.

If called with `dim`, it integrates entries tagged with `:Perm` via `SymbolicMatrix(:P, :Perm)`.

The integration of a monomial of entries is given by:
```math
\int_{S_d} P_{i_1 j_1} \dots P_{i_n j_n} dP = \begin{cases} \frac{(d-k)!}{d!} & \text{if indices are consistent} \\ 0 & \text{otherwise} \end{cases}
```
where $k$ is the number of distinct pairs $(i, j)$ in the product.
"""
dPerm(dim) = PermutationMeasure(dim)
dPerm(P::SymbolicMatrix) = PermutationMeasure(P.dim)

"""
    dCPerm(dim)
    dCPerm(Y::SymbolicMatrix)

Defines the measure for Centered Permutation matrices $Y = P - J/d$ where $P \in S_d$.

If called with `dim`, it integrates entries tagged with `:CPerm` (or `:Perm`) via `SymbolicMatrix(:Y, :CPerm)`.
"""
dCPerm(dim) = CenteredPermutationMeasure(dim)
dCPerm(Y::SymbolicMatrix) = CenteredPermutationMeasure(Y.dim)

function integrate(expr::AbstractArray, measure::PermutationMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::CenteredPermutationMeasure)
    return map(e -> integrate(e, measure), expr)
end

function IntU.measure_info(measure::PermutationMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = MetadataMatcher(:Perm)
    return (subs_dict, matcher, measure.dim, :Perm)
end

function IntU.measure_info(measure::CenteredPermutationMeasure)
    # Centered permutations are handled by substituting Y_ij = P_ij - 1/dim
    # We need a matcher for P, which is the Permutation group.
    # MetadataMatcher(:Perm) will match P_ij.
    # But wait, centered permutation entries Y_ij might not be tagged by the user.
    # Actually, if the user creates Y = SymbolicMatrix(:Y, :Perm), getindex will tag it.
    
    subs_dict = Dict{Any,Any}()
    # If expr contains Y_ij, we need to substitute it with P_ij - 1/dim.
    # This might require some more elaborate logic in process_term if we want it fully automatic.
    # For now, let's assume Y_ij is tagged as :Perm.
    
    matcher = MetadataMatcher(:Perm)
    return (subs_dict, matcher, measure.dim, :Perm)
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
