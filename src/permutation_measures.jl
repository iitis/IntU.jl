# Permutation Group measures

# Dummy types to represent the measures
struct PermutationMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end
PermutationMeasure(dim) = PermutationMeasure(dim, nothing)

struct CenteredPermutationMeasure{D,M} <: AbstractMeasure
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
    matcher = measure.matcher === nothing ? MetadataMatcher(:CPerm) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, :CPerm)
end

function fallback_integrate(t::LazyTrace, measure::PermutationMeasure)
    # E[tr(PA)] = sum(A) / d for one P.
    # For more complex terms, we expand to element-wise integration.
    matcher = measure.matcher === nothing ? MetadataMatcher(:Perm) : measure.matcher
    dim = measure.dim

    # If it's a simple tr(PA), keep the symbolic result if A is not concrete
    if length(t.cycles) == 1 && length(t.cycles[1]) == 2
        factors = t.cycles[1]
        P_idx = nothing
        for (i, f) in enumerate(factors)
            if match_index(matcher, f) !== nothing
                P_idx = i;
                break
            end
        end
        if P_idx !== nothing
            A = factors[P_idx == 1 ? 2 : 1]
            # If A is a SymbolicMatrix or a simple numeric matrix, we return a symbolic sum result
            if A isa SymbolicMatrix || !(A isa AbstractMatrix && !(eltype(A) <: Num))
                return t.prefactor *
                       (Symbolics.variable(Symbol("sum(" * string(A) * ")")) / measure.dim)
            end
        end
    end

    # General expansion
    expr = t.prefactor
    for cycle in t.cycles
        # Each cycle is tr(ABC...)
        # Expand tr(ABC...) as sum_{i,j,k...} A_ij B_jk C_ki
        n = length(cycle)
        dims = [size(f) for f in cycle]

        # We need a shared dimension for all factors in the cycle for trace to exist.
        # But for symbolic ones it is dim.
        d_val = dim isa Integer ? Int(dim) : 0
        if d_val == 0
            # Try to find a concrete dimension from any factor
            for f in cycle
                s = size(f, 1)
                if s isa Integer && s < 1000 # heuristic
                    d_val = s;
                    break
                end
            end
        end

        if d_val == 0
            error(
                "Cannot expand LazyTrace for Permutations: dimension is not concrete and term is not linear.",
            )
        end

        # Manual expansion of tr(C1 * C2 * ... * Cn)
        # sum_{i1, i2, ..., in} C1[i1, i2] * C2[i2, i3] * ... * Cn[in, i1]
        indices = Symbolics.variable.(Symbol.("i", 1:n), T = Int) # Not really needed as symbols if we just use loop

        term_sum = 0
        for idxs in Iterators.product(fill(1:d_val, n)...)
            prod_val = 1
            for k = 1:n
                next_k = (k % n) + 1
                prod_val *= cycle[k][idxs[k], idxs[next_k]]
            end
            term_sum += prod_val
        end
        expr *= term_sum
    end

    return integrate(expr, measure)
end

function fallback_integrate(t::LazyTrace, measure::CenteredPermutationMeasure)
    # E[tr(YA)] = 0 because E[P - J/d] = 0.
    if length(t.cycles) == 1 && length(t.cycles[1]) == 2
        return 0
    end
    error(
        "Graphical integration for Centered Permutations only supported for tr(YA) currently.",
    )
end

"""
    integrate_indices_permutation(indices, dim)

Integration over the symmetric group S_d.
"""
function integrate_indices_permutation(indices::AbstractVector{Tuple{Int,Int}}, dim)
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
