
module IntUITensorsExt

using IntU
using ITensors
using Symbolics
import IntU:
    integrate,
    _contract_all,
    _create_deltas,
    _contract_with_deltas,
    _contract_scalar_with_deltas,
    _delta_elem_type,
    _scalar_coeff_constant_across_pairs,
    _supports_scalar_fastpath,
    _wrap_scalar_graphical_result,
    GraphicalUnitary,
    ITensorUnitary

# Integration for a single ITensorUnitary
function integrate(u::ITensorUnitary, measure::IntU.HaarMeasure)
    return integrate([u], measure)
end

# Typed graphical entry points for ITensor constants (avoid generic Any dispatch).
IntU.integrate_graphical(
    constants::AbstractVector{ITensor},
    unitaries,
    measure::IntU.HaarMeasure,
) = IntU._integrate_graphical_unitary(constants, unitaries, measure.dim)

IntU.integrate_graphical(
    constants::AbstractVector{ITensor},
    unitaries,
    measure::IntU.UnitaryDesign,
) = IntU._integrate_graphical_unitary(constants, unitaries, measure.dim; design_t = measure.t)

# Main entry points for a network of tensors
function integrate(tensors::AbstractVector, measure::IntU.AbstractMeasure)
    # Guard: only treat as a tensor network if it contains ITensors or ITensorUnitary
    is_tensor_network = any(t -> t isa ITensor || t isa ITensorUnitary, tensors)

    if is_tensor_network
        return _integrate_tensor_network(tensors, measure)
    else
        # Fallback to standard element-wise integration from IntU
        # we use invoke to call the more general AbstractArray method from IntU
        return invoke(
            integrate,
            Tuple{AbstractArray,IntU.AbstractMeasure},
            tensors,
            measure,
        )
    end
end

function _integrate_tensor_network(tensors::AbstractVector, measure)
    # Identify random unitaries and constants
    unitaries = GraphicalUnitary[]
    constants = ITensor[]

    for t in tensors
        if t isa ITensorUnitary
            push!(
                unitaries,
                GraphicalUnitary(
                    collect(Any, t.out_indices),
                    collect(Any, t.in_indices),
                    t.is_adj,
                ),
            )
        elseif t isa ITensor
            push!(constants, t)
        elseif t isa Union{Number,IntU.Num}
            push!(constants, ITensor(t))
        else
            # Try to handle other types if possible, or error
            error("Unknown tensor type: $(typeof(t))")
        end
    end

    return IntU.integrate_graphical(constants, unitaries, measure)
end

# Specific overloads for ITensor vectors to ensure they hit the extension
function integrate(tensors::AbstractVector{<:ITensor}, measure::IntU.AbstractMeasure)
    return integrate(collect(Any, tensors), measure)
end

# Implementation of internal hooks for ITensors

function _contract_all(cs::AbstractVector{ITensor})
    if isempty(cs)
        return 1.0
    end
    res = cs[1]
    for i = 2:length(cs)
        res = res * cs[i]
    end
    return res
end

function _create_deltas(idxs1, idxs2)
    # idxs1 and idxs2 are Vector{Index}
    # Create a delta for each pair
    if length(idxs1) != length(idxs2)
        error("Index mismatch in delta creation")
    end
    return [delta(idxs1[i], idxs2[i]) for i = 1:length(idxs1)]
end

_delta_elem_type(::AbstractVector{ITensor}) = ITensor

function _merge_tensor_lists(
    cs::AbstractVector{ITensor},
    ds::AbstractVector{<:ITensor},
)
    n_cs = length(cs)
    n_ds = length(ds)
    n_all = n_cs + n_ds
    if n_all == 0
        return ITensor[]
    end
    all_tensors = Vector{ITensor}(undef, n_all)
    copyto!(all_tensors, 1, cs, 1, n_cs)
    copyto!(all_tensors, n_cs + 1, ds, 1, n_ds)
    return all_tensors
end

function _contract_with_deltas(
    cs::AbstractVector{ITensor},
    ds::AbstractVector{<:ITensor},
    wg,
)
    all_tensors = _merge_tensor_lists(cs, ds)
    if isempty(all_tensors)
        return wg
    end

    res = ITensors.contract(all_tensors)

    return wg * res
end

function _contract_with_deltas(cs::AbstractVector{ITensor}, ds::AbstractVector, wg)
    typed_ds = ITensor[delta_t for delta_t in ds]
    return _contract_with_deltas(cs, typed_ds, wg)
end

function _is_scalar_itensor(t::ITensor)
    return isempty(inds(t))
end

function _canonicalize_scalar_coeff(x)
    x_unwrapped = x isa IntU.Num ? Symbolics.unwrap(x) : x

    if x_unwrapped isa Integer || x_unwrapped isa Rational
        return x_unwrapped
    end

    if x_unwrapped isa AbstractFloat
        if isfinite(x_unwrapped) && isinteger(x_unwrapped)
            return round(BigInt, x_unwrapped)
        end
        return x
    end

    if x_unwrapped isa Real
        if isinteger(x_unwrapped)
            return round(BigInt, BigFloat(x_unwrapped))
        end
        return x
    end

    return x
end

function _supports_scalar_fastpath(
    cs::AbstractVector{ITensor},
    ds::AbstractVector{<:ITensor},
)
    res = ITensors.contract(_merge_tensor_lists(cs, ds))
    return _is_scalar_itensor(res)
end

function _supports_scalar_fastpath(cs::AbstractVector{ITensor}, ds::AbstractVector)
    typed_ds = ITensor[delta_t for delta_t in ds]
    return _supports_scalar_fastpath(cs, typed_ds)
end

function _contract_scalar_with_deltas(
    cs::AbstractVector{ITensor},
    ds::AbstractVector{<:ITensor},
)
    res = ITensors.contract(_merge_tensor_lists(cs, ds))
    _is_scalar_itensor(res) ||
        throw(ArgumentError("Scalar fast path requested for a non-scalar ITensor contraction"))
    return _canonicalize_scalar_coeff(ITensors.scalar(res))
end

function _contract_scalar_with_deltas(cs::AbstractVector{ITensor}, ds::AbstractVector)
    typed_ds = ITensor[delta_t for delta_t in ds]
    return _contract_scalar_with_deltas(cs, typed_ds)
end

_wrap_scalar_graphical_result(::AbstractVector{ITensor}, scalar) = ITensor(scalar)

function _all_indices_dim_one(constants::AbstractVector{ITensor})
    for t in constants
        for idx in inds(t)
            dim(idx) == 1 || return false
        end
    end
    return true
end

function _all_legs_dim_one(legs)
    for leg_group in legs
        for idx in leg_group
            dim(idx) == 1 || return false
        end
    end
    return true
end

function _scalar_coeff_constant_across_pairs(
    constants::AbstractVector{ITensor},
    u_out,
    u_in,
    u_dag_out,
    u_dag_in,
)
    return _all_indices_dim_one(constants) &&
           _all_legs_dim_one(u_out) &&
           _all_legs_dim_one(u_in) &&
           _all_legs_dim_one(u_dag_out) &&
           _all_legs_dim_one(u_dag_in)
end

function IntU._create_deltas_symplectic(idxs1, idxs2, dim)
    # Symplectic contraction involves J
    # J_{ab} = δ_{a, b+n} - δ_{a+n, b} where n = dim/2
    # In ITensors, we can't easily express J as a single delta if it's not diagonal.
    # But J is a constant tensor. We can create it.

    if length(idxs1) != length(idxs2)
        error("Index mismatch in symplectic delta creation")
    end

    deltas = ITensor[]
    for i = 1:length(idxs1)
        # Create a J tensor for this pair of indices
        # J = delta(idxs1[i], idxs2[i] + n) - delta(idxs1[i] + n, idxs2[i])
        # This requires the indices to be of the form that allows adding n (Integer-like).
        # In ITensors, we usually use `replaceind` or just know the name.

        # A more robust way in ITensors:
        # Create a J tensor with indices idxs1[i] and idxs2[i]
        # We need to access index values.

        n_half = Int(dim ÷ 2)
        j_tensor = ITensor(idxs1[i], idxs2[i])
        for val = 1:Int(dim)
            if val <= n_half
                j_tensor[idxs1[i]=>val, idxs2[i]=>val+n_half] = 1.0
            else
                j_tensor[idxs1[i]=>val, idxs2[i]=>val-n_half] = -1.0
            end
        end
        push!(deltas, j_tensor)
    end
    return deltas
end

# Overload for a product of tensors (if passed as T1 * T2 * ...)
# ITensors doesn't have a specific type for "uncontracted product of ITensors" 
# other than Vector{ITensor} or just the result of contraction.
# But we usually want to integrate BEFORE contraction.

end # module
