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

Defines the Haar measure for the Symmetric group $S_d$ (permutation matrices) of dimension `dim`.

Integration engine identifies variables via metadata tag `:Perm`.
"""
dPerm(dim) = PermutationMeasure(dim)

@doc raw"""
    dCPerm(dim)

Defines the measure for Centered Permutation matrices $Y = P - J/d$ where $P \in S_d$.
"""
dCPerm(dim) = CenteredPermutationMeasure(dim)

function integrate(expr::AbstractArray, measure::PermutationMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::CenteredPermutationMeasure)
    return map(e -> integrate(e, measure), expr)
end

# Resolve ambiguities with SymbolicMatrix/SymbolicMatrixProduct
function integrate(expr::SymbolicMatrix, measure::PermutationMeasure)
    return invoke(integrate, Tuple{SymbolicMatrix, Any}, expr, measure)
end
function integrate(expr::SymbolicMatrixProduct, measure::PermutationMeasure)
    return invoke(integrate, Tuple{SymbolicMatrixProduct, Any}, expr, measure)
end
function integrate(expr::SymbolicMatrix, measure::CenteredPermutationMeasure)
    return invoke(integrate, Tuple{SymbolicMatrix, Any}, expr, measure)
end
function integrate(expr::SymbolicMatrixProduct, measure::CenteredPermutationMeasure)
    return invoke(integrate, Tuple{SymbolicMatrixProduct, Any}, expr, measure)
end

function IntU.measure_info(measure::PermutationMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = MetadataMatcher(:Perm)
    return (subs_dict, matcher, measure.dim, :Perm)
end

function IntU.measure_info(measure::CenteredPermutationMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = MetadataMatcher(:Perm)
    return (subs_dict, matcher, measure.dim, :CPerm)
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
