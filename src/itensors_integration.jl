# src/itensors_integration.jl

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
    integrate_graphical(constants, unitaries, dim, measure_type=:U)

Integrates a tensor network.
`constants` is a list of tensors (e.g. ITensors).
`unitaries` is a list of `GraphicalUnitary`.
`dim` is the dimension of the space.
`measure_type` can be `:U`, `:O`, or `:Sp`.

Returns a sum of terms, where each term is a product of `constants` and deltas.
"""
function integrate_graphical(constants, unitaries, dim, measure_type = :U)
    if measure_type == :U
        return _integrate_graphical_unitary(constants, unitaries, dim)
    elseif measure_type == :O
        return _integrate_graphical_orthogonal(constants, unitaries, dim)
    elseif measure_type == :Sp
        return _integrate_graphical_symplectic(constants, unitaries, dim)
    elseif measure_type isa Tuple && measure_type[1] == :Design
        return _integrate_graphical_unitary(
            constants,
            unitaries,
            dim;
            design_t = measure_type[2],
        )
    else
        error("Unsupported measure type: $measure_type")
    end
end

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
        error("Integrand degree ($n_u, $n_u) exceeds design order t=$design_t")
    end


    if n_u == 0
        # Just return the contraction of constants
        return _contract_all(constants)
    end

    # 2. Iterate over permutations sigma, tau in S_n
    perms = collect(permutations(1:n_u))

    total_result = nothing

    for sigma_p in perms
        for tau_p in perms
            # sigma and tau are permutations of 1..n
            # Cycle type of sigma * tau^-1
            # We can use IntU's existing get_cycle_type if we have the combined permutation
            # sigma_p * inv(tau_p)
            P = [sigma_p[invperm(tau_p)[i]] for i = 1:n_u]
            cycle_type = get_cycle_type(P)
            wg_val = weingarten(cycle_type, dim)

            if _iszero(wg_val)
                continue
            end

            # Create deltas
            # For each k, match u_list[k].out with u_dag_list[sigma_p[k]].out
            # and u_list[k].in with u_dag_list[tau_p[k]].in

            deltas = []
            for k = 1:n_u
                # Out matchings
                append!(
                    deltas,
                    _create_deltas(
                        u_list[k].out_indices,
                        u_dag_list[sigma_p[k]].out_indices,
                    ),
                )
                # In matchings
                append!(
                    deltas,
                    _create_deltas(u_list[k].in_indices, u_dag_list[tau_p[k]].in_indices),
                )
            end

            # The result for this permutation pair is wg_val * constants * deltas
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

function _integrate_graphical_orthogonal(constants, unitaries, dim)
    # Orthogonal: all unitaries are treated same (O = O_bar)
    # n_total must be even
    n_total = length(unitaries)
    if n_total == 0
        return _contract_all(constants)
    end
    if n_total % 2 != 0
        return 0
    end

    k = div(n_total, 2)
    # Get all pair partitions of 1..2k
    partitions = get_pair_partitions(n_total)

    total_result = nothing

    for pi in partitions
        for sigma in partitions
            wg_val = weingarten_orthogonal_val(pi, sigma, dim)
            if _iszero(wg_val)
                continue
            end

            # Create deltas for pi (out) and sigma (in)
            # Actually for Orthogonal, we match legs within members of the partition
            deltas = []
            # Combine all unitaries' legs into a single list
            all_out = []
            all_in = []
            for u in unitaries
                push!(all_out, u.out_indices)
                push!(all_in, u.in_indices)
            end

            # pi specifies which out legs to match
            for (a, b) in pi
                append!(deltas, _create_deltas(all_out[a], all_out[b]))
            end
            # sigma specifies which in legs to match
            for (a, b) in sigma
                append!(deltas, _create_deltas(all_in[a], all_in[b]))
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

function _integrate_graphical_symplectic(constants, unitaries, dim)
    # Symplectic: similar to orthogonal but with symplectic form J
    n_total = length(unitaries)
    if n_total == 0
        return _contract_all(constants)
    end
    if n_total % 2 != 0
        return 0
    end

    partitions = get_pair_partitions(n_total)
    total_result = nothing

    for pi in partitions
        for sigma in partitions
            wg_val = weingarten_symplectic_val(pi, sigma, dim)
            if _iszero(wg_val)
                continue
            end

            # Matchings involve the symplectic form J
            # In graphical terms, this means adding a J tensor between legs
            # For simplicity, we can treat J as part of _create_deltas_symplectic

            deltas = []
            all_out = [u.out_indices for u in unitaries]
            all_in = [u.in_indices for u in unitaries]

            for (a, b) in pi
                append!(deltas, _create_deltas_symplectic(all_out[a], all_out[b], dim))
            end
            for (a, b) in sigma
                append!(deltas, _create_deltas_symplectic(all_in[a], all_in[b], dim))
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

# Additional symplectic helper
function _create_deltas_symplectic end

# These helpers will be specialized or provided by the extension if they depend on ITensors
function _contract_all end
function _create_deltas end
function _contract_with_deltas end
