# Permutation Group measures

# Dummy types to represent the measures
struct PermutationMeasure{D,M}
    dim::D
    matcher::M
end
PermutationMeasure(dim) = PermutationMeasure(dim, nothing)

struct CenteredPermutationMeasure{D,M}
    dim::D
    matcher::M
end
CenteredPermutationMeasure(dim) = CenteredPermutationMeasure(dim, nothing)

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
    matcher = measure.matcher === nothing ? MetadataMatcher(:Perm) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, :Perm)
end

function IntU.measure_info(measure::CenteredPermutationMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = measure.matcher === nothing ? MetadataMatcher(:Perm) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, :CPerm)
end

function fallback_integrate(t::LazyTrace, measure::PermutationMeasure)
    # E[tr(PA)] = sum(A) / d for one P.
    # For now, let's implement a very basic expansion for tr(PA)
    if length(t.cycles) == 1 && length(t.cycles[1]) == 2
        factors = t.cycles[1]
        matcher = measure.matcher === nothing ? MetadataMatcher(:Perm) : measure.matcher
        
        # Check if one of factors is P
        P_idx = nothing
        for (i, f) in enumerate(factors)
            if match_index(matcher, f) !== nothing
                P_idx = i; break
            end
        end
        
        if P_idx !== nothing
            A = factors[P_idx == 1 ? 2 : 1]
            # Result is t.prefactor * sum(A) / measure.dim
            # We can represent sum(A) as a new symbolic variable or expand if small.
            # For symbolic A, let's return a trace-like name or sum(A)
            return t.prefactor * (Symbolics.variable(Symbol("sum(" * string(A) * ")")) / measure.dim)
        end
    end
    error("Graphical integration for Permutations only supported for tr(PA) currently.")
end

function fallback_integrate(t::LazyTrace, measure::CenteredPermutationMeasure)
    # E[tr(YA)] = 0 because E[P - J/d] = 0.
    if length(t.cycles) == 1 && length(t.cycles[1]) == 2
         return 0
    end
    error("Graphical integration for Centered Permutations only supported for tr(YA) currently.")
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
