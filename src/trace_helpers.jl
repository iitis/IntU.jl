
"""
    _extract_trace_data(t::LazyTrace)

Extract cycle ranges and all factors from a `LazyTrace`.
Returns `(total_factors, cycle_ranges, all_factors)`.
"""
function _extract_trace_data(t::LazyTrace)
    total_factors = 0
    cycle_ranges = UnitRange{Int}[]
    all_factors = Any[]
    for cycle in t.cycles
        start_idx = total_factors + 1
        append!(all_factors, cycle)
        total_factors += length(cycle)
        push!(cycle_ranges, start_idx:total_factors)
    end
    return (total_factors, cycle_ranges, all_factors)
end

"""
    _build_wires(U_indices, U_bar_indices, cycle_ranges, all_factors)

Build forward and backward wire connectivity maps for graphical integration.
Returns `(wires, reverse_wires)` where each maps a slot index to
`(next_slot, constant_matrices_between)`.
"""
function _build_wires(U_indices, U_bar_indices, cycle_ranges, all_factors)
    wires = Dict{Int,Any}()
    reverse_wires = Dict{Int,Any}()
    all_slots = sort([U_indices; U_bar_indices])
    all_slots_set = Set{Int}(all_slots)

    for slot in all_slots
        cid = findfirst(rng -> slot in rng, cycle_ranges)
        rng = cycle_ranges[cid]

        curr = slot
        consts = Any[]
        while true
            curr = curr == last(rng) ? first(rng) : curr + 1
            if curr in all_slots_set
                wires[slot] = (curr, isempty(consts) ? nothing : consts)
                break
            end
            push!(consts, all_factors[curr])
        end

        curr = slot
        consts_rev = Any[]
        while true
            curr = curr == first(rng) ? last(rng) : curr - 1
            if curr in all_slots_set
                reverse_wires[slot] = (curr, isempty(consts_rev) ? nothing : consts_rev)
                break
            end
            push!(consts_rev, transpose(all_factors[curr]))
        end
    end
    return wires, reverse_wires
end

"""
    _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

Evaluate trace cycles that contain no integration variables (constant cycles).
Returns the product of the prefactor and all constant cycle traces.
"""
function _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    res = t.prefactor
    slot_set = Set{Int}(all_slots)
    for (cid, rng) in enumerate(cycle_ranges)
        has_U = any(idx -> idx in slot_set, rng)

        if !has_U
            cycle = t.cycles[cid]
            if isempty(cycle)
                res *= dim
            else
                res *= tr_val(cycle)
            end
        end
    end
    return res
end
