
module IntUITensorsExt

using IntU
using ITensors
import IntU:
    integrate,
    _contract_all,
    _create_deltas,
    _contract_with_deltas,
    GraphicalUnitary,
    ITensorUnitary

# Integration for a single ITensorUnitary
function integrate(u::ITensorUnitary, measure::IntU.HaarMeasure)
    return integrate([u], measure)
end

# Main entry points for a network of tensors
function integrate(tensors::AbstractVector, measure::IntU.AbstractMeasure)
    _integrate_tensor_network(tensors, measure)
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

function _contract_with_deltas(cs::AbstractVector{ITensor}, ds::AbstractVector, wg)
    # ds is a list of deltas (which are ITensors)
    # Contract everything
    all_tensors = vcat(cs, ds)
    if isempty(all_tensors)
        return wg
    end

    res = ITensors.contract(all_tensors)

    return wg * res
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
