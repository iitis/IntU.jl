
"""
    GraphicalUnitary

Represents a random unitary tensor in a graphical network.
`out_indices` and `in_indices` are lists of index objects (e.g. ITensor Index).
"""
struct GraphicalUnitary
    out_indices::Vector{Any}
    in_indices::Vector{Any}
    is_adj::Bool
end
"""
    ITensorUnitary(tensor; out_indices, in_indices, is_adj=false)

Wrapper to mark an ITensor (or any other tensor object) as a random unitary for integration.
"""
struct ITensorUnitary{T}
    tensor::T
    out_indices::AbstractVector
    in_indices::AbstractVector
    is_adj::Bool
end

# Outer constructor for keyword arguments
function ITensorUnitary(tensor; out_indices, in_indices, is_adj = false)
    return ITensorUnitary(tensor, out_indices, in_indices, is_adj)
end

function ITensorUnitary(; out_indices, in_indices, is_adj = false)
    return ITensorUnitary(nothing, out_indices, in_indices, is_adj)
end

"""
    integrate_graphical(constants, unitaries, measure::AbstractMeasure)

Integrates a tensor network.
`constants` is a list of tensors (e.g. ITensors).
`unitaries` is a list of `GraphicalUnitary`.
`measure` is an `AbstractMeasure` (e.g. `HaarMeasure`, `OrthogonalMeasure`).

Returns a sum of terms, where each term is a product of `constants` and deltas.
"""
function integrate_graphical(constants, unitaries, measure::AbstractMeasure)
    # Default fallback or error
    throw(
        ArgumentError(
            "Graphical integration not implemented for measure $(typeof(measure))",
        ),
    )
end

integrate_graphical(constants, unitaries, measure::HaarMeasure) =
    _integrate_graphical_unitary(constants, unitaries, measure.dim)

integrate_graphical(constants, unitaries, measure::UnitaryDesign) =
    _integrate_graphical_unitary(constants, unitaries, measure.dim; design_t = measure.t)

integrate_graphical(constants, unitaries, measure::OrthogonalMeasure) =
    _integrate_graphical_real(constants, unitaries, measure)

integrate_graphical(constants, unitaries, measure::SymplecticMeasure) =
    _integrate_graphical_real(constants, unitaries, measure)

const _UNITARY_GRAPH_CACHE = Dict{Int,NamedTuple}()

function _get_unitary_graph_cache(n_u::Int)
    return get!(_UNITARY_GRAPH_CACHE, n_u) do
        perms = collect(permutations(1:n_u))
        n_perm = length(perms)
        inv_perms = [invperm(p) for p in perms]

        composed = Matrix{Vector{Int}}(undef, n_perm, n_perm)
        cycle_ids = Matrix{Int}(undef, n_perm, n_perm)
        cycle_types = Vector{Vector{Int}}()
        cycle_lookup = Dict{Tuple{Vararg{Int}},Int}()

        for i = 1:n_perm
            sigma = perms[i]
            for j = 1:n_perm
                tau_inv = inv_perms[j]
                comp = [sigma[tau_inv[k]] for k = 1:n_u]
                composed[i, j] = comp
                ct = get_cycle_type(comp)
                key = Tuple(ct)
                cid = get(cycle_lookup, key, 0)
                if cid == 0
                    push!(cycle_types, ct)
                    cid = length(cycle_types)
                    cycle_lookup[key] = cid
                end
                cycle_ids[i, j] = cid
            end
        end

        return (
            perms = perms,
            composed = composed,
            cycle_ids = cycle_ids,
            cycle_types = cycle_types,
        )
    end
end

function _create_unitary_delta_pairs(u_idxs, u_bar_idxs)
    n = length(u_idxs)
    pairs = Matrix{Any}(undef, n, n)
    for i = 1:n
        for j = 1:n
            pairs[i, j] = _create_deltas(u_idxs[i], u_bar_idxs[j])
        end
    end
    return pairs
end

_delta_elem_type(constants) = Any

function _new_deltas_buffer(constants, n_u::Int)
    T = _delta_elem_type(constants)
    deltas = T === Any ? Any[] : Vector{T}()
    sizehint!(deltas, max(2 * n_u, 4))
    return deltas
end

function _fill_unitary_deltas!(
    deltas::AbstractVector,
    out_delta_pairs,
    in_delta_pairs,
    sigma_p::AbstractVector{<:Integer},
    tau_p::AbstractVector{<:Integer},
)
    empty!(deltas)
    n_u = length(sigma_p)
    for k = 1:n_u
        append!(deltas, out_delta_pairs[k, sigma_p[k]])
        append!(deltas, in_delta_pairs[k, tau_p[k]])
    end
    return deltas
end

_supports_scalar_fastpath(constants, deltas) = false

function _contract_scalar_with_deltas(constants, deltas)
    throw(ArgumentError("Scalar contraction fast path is not implemented for constants of type $(typeof(constants))"))
end

_wrap_scalar_graphical_result(constants, scalar) = scalar
_scalar_coeff_constant_across_pairs(constants, u_out, u_in, u_dag_out, u_dag_in) = false

function _integrate_graphical_unitary(constants, unitaries, dim; design_t = nothing)
    # 1. Separate U and U_dag
    u_list = filter(u -> !u.is_adj, unitaries)
    u_dag_list = filter(u -> u.is_adj, unitaries)

    n_u = length(u_list)
    n_dag = length(u_dag_list)

    if n_u != n_dag
        return 0 # Or a symbolic zero
    end

    if design_t !== nothing && n_u > design_t
        throw(
            ArgumentError("Integrand degree ($n_u, $n_u) exceeds design order t=$design_t"),
        )
    end


    if n_u == 0
        # Just return the contraction of constants
        return _contract_all(constants)
    end

    # Pre-extract indices to avoid repeated access and allocations
    u_out = [u.out_indices for u in u_list]
    u_in = [u.in_indices for u in u_list]
    u_dag_out = [u.out_indices for u in u_dag_list]
    u_dag_in = [u.in_indices for u in u_dag_list]

    cache = _get_unitary_graph_cache(n_u)
    perms = cache.perms
    cycle_ids = cache.cycle_ids
    cycle_types = cache.cycle_types
    n_perm = length(perms)

    # Precompute Wg value for each unique cycle type once per call.
    wg_by_cycle = Vector{Any}(undef, length(cycle_types))
    wg_is_nonzero = Vector{Bool}(undef, length(cycle_types))
    for cid = 1:length(cycle_types)
        wg_val = weingarten(cycle_types[cid], dim)
        wg_by_cycle[cid] = wg_val
        wg_is_nonzero[cid] = !_iszero(wg_val)
    end

    active_pairs = Tuple{Int,Int,Int}[]
    for i = 1:n_perm
        for j = 1:n_perm
            cid = cycle_ids[i, j]
            wg_is_nonzero[cid] || continue
            push!(active_pairs, (i, j, cid))
        end
    end

    isempty(active_pairs) && return 0

    # Precompute delta tensors for every (U_k, U†_l) matching and reuse.
    out_delta_pairs = _create_unitary_delta_pairs(u_out, u_dag_out)
    in_delta_pairs = _create_unitary_delta_pairs(u_in, u_dag_in)

    deltas = _new_deltas_buffer(constants, n_u)
    first_i, first_j, _ = active_pairs[1]
    _fill_unitary_deltas!(deltas, out_delta_pairs, in_delta_pairs, perms[first_i], perms[first_j])

    # Closed-scalar networks can be accumulated on scalars and wrapped at the end.
    if _supports_scalar_fastpath(constants, deltas)
        coeff_by_cycle = Dict{Int,Any}()
        if _scalar_coeff_constant_across_pairs(constants, u_out, u_in, u_dag_out, u_dag_in)
            coeff0 = _contract_scalar_with_deltas(constants, deltas)
            cycle_counts = Dict{Int,Int}()
            for (_, _, cid) in active_pairs
                cycle_counts[cid] = get(cycle_counts, cid, 0) + 1
            end
            for (cid, cnt) in cycle_counts
                coeff_by_cycle[cid] = cnt * coeff0
            end
        else
            for (i, j, cid) in active_pairs
                _fill_unitary_deltas!(deltas, out_delta_pairs, in_delta_pairs, perms[i], perms[j])
                coeff = _contract_scalar_with_deltas(constants, deltas)
                if haskey(coeff_by_cycle, cid)
                    coeff_by_cycle[cid] += coeff
                else
                    coeff_by_cycle[cid] = coeff
                end
            end
        end

        total_scalar = nothing
        for (cid, coeff) in coeff_by_cycle
            term = coeff * wg_by_cycle[cid]
            if total_scalar === nothing
                total_scalar = term
            else
                total_scalar += term
            end
        end
        return _wrap_scalar_graphical_result(constants, total_scalar === nothing ? 0 : total_scalar)
    end

    total_result = nothing
    for (i, j, cid) in active_pairs
        _fill_unitary_deltas!(deltas, out_delta_pairs, in_delta_pairs, perms[i], perms[j])
        term = _contract_with_deltas(constants, deltas, wg_by_cycle[cid])
        if total_result === nothing
            total_result = term
        else
            total_result = total_result + term
        end
    end

    return total_result === nothing ? 0 : total_result
end

function _integrate_graphical_real(constants, unitaries, measure::AbstractMeasure)
    # Orthogonal/Symplectic: all unitaries are treated same (O = O_bar)
    # n_total must be even
    n_total = length(unitaries)
    if n_total == 0
        return _contract_all(constants)
    end
    if n_total % 2 != 0
        return 0
    end

    dim = measure.dim
    k = div(n_total, 2)
    # Get all pair partitions of 1..2k
    partitions = get_pair_partitions(n_total)

    # Pre-extract indices
    all_out = [u.out_indices for u in unitaries]
    all_in = [u.in_indices for u in unitaries]

    # Pre-canonicalize partitions for faster Wg lookup
    c_partitions = [canonicalize_pair_partition(p) for p in partitions]
    total_result = nothing

    for i = 1:length(partitions)
        pi = partitions[i]
        c_pi = c_partitions[i]
        for j = 1:length(partitions)
            sigma = partitions[j]
            c_sigma = c_partitions[j]

            wg_val = _weingarten_real(measure, c_pi, c_sigma, dim)

            if _iszero(wg_val)
                continue
            end

            deltas = []
            # pi specifies which out legs to match
            for (a, b) in pi
                append!(deltas, _create_deltas_general(measure, all_out[a], all_out[b]))
            end
            # sigma specifies which in legs to match
            for (a, b) in sigma
                append!(deltas, _create_deltas_general(measure, all_in[a], all_in[b]))
            end

            term = _contract_with_deltas(constants, deltas, wg_val)
            if total_result === nothing
                total_result = term
            else
                total_result = total_result + term
            end
        end
    end
    return total_result === nothing ? 0 : total_result
end

# Helpers for real integration dispatch
_weingarten_real(::OrthogonalMeasure, c_pi, c_sigma, dim) =
    weingarten_orthogonal_val_canonical(c_pi, c_sigma, dim)
_weingarten_real(::SymplecticMeasure, c_pi, c_sigma, dim) =
    weingarten_symplectic_val(c_pi, c_sigma, dim)

_create_deltas_general(::AbstractMeasure, idxs1, idxs2) = _create_deltas(idxs1, idxs2)
_create_deltas_general(m::SymplecticMeasure, idxs1, idxs2) =
    _create_deltas_symplectic(idxs1, idxs2, m.dim)

# Additional symplectic helper
function _create_deltas_symplectic end

# These helpers will be specialized or provided by the extension if they depend on ITensors
function _contract_all end
function _create_deltas end
function _contract_with_deltas end
