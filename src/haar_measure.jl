# Haar measure integration

# Matcher for SymbolicMatrix
struct SymbolicMatcher <: AbstractIndexMatcher
    tag::Symbol
    regex::Regex
end

function match_index(m::SymbolicMatcher, t)
    # Unwrap to Sym
    s = Symbolics.unwrap(t)
    # We check string representation of the symbol. Handle conj(U_i_j)
    is_conj = false
    if Symbolics.iscall(s) && (Symbolics.operation(s) == conj || Symbolics.operation(s) == Base.conj)
        is_conj = true
        s = Symbolics.arguments(s)[1]
    end
    
    s_str = string(s)
    mat = match(m.regex, s_str)
    if mat !== nothing
        try
            i = parse(Int, mat[1])
            j = parse(Int, mat[2])
            final_tag = is_conj ? :U_bar : :U
            return (final_tag, i, j)
        catch
        end
    end
    return nothing
end


# Dummy type to represent the measure
struct HaarMeasure{D}
    dim::D
end
@doc raw"""
    dU(dim)

Defines the Haar measure for the Unitary group $U(d)$.

Integration engine identifies variables via metadata tag `:U`.
"""
dU(dim) = HaarMeasure(dim)


"""
    measure_info(measure)

Returns a tuple `(matcher, dim, type)` for the given measure. 
Internal function used for dispatching integration logic.
"""
function IntU.measure_info(measure::HaarMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = MetadataMatcher(:U)
    return (subs_dict, matcher, measure.dim, :U)
end

function _manual_fallback(expr, measure::HaarMeasure)
    # LazyTrace expressions are handled by fallback_integrate dispatch, not here.
    error("HaarMeasure integration failed for: $(typeof(expr))")
end

"""
    asymptotic(expr, measure::HaarMeasure, order=1)

Returns the series expansion of the integral in powers of `1/d`.
"""
function asymptotic(expr, measure::HaarMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dU(d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end

# Integrate a product of traces of matrices over the Haar measure.
# Uses the graphical Weingarten calculus.
function fallback_integrate(t::LazyTrace, measure::HaarMeasure)
    dim = measure.dim
    constant_part = t.prefactor

    # Identify ALL cycles and factors
    total_factors = 0
    cycle_ranges = UnitRange{Int}[]
    all_factors = Any[]
    for cycle in t.cycles
        start_idx = total_factors + 1
        append!(all_factors, cycle)
        total_factors += length(cycle)
        push!(cycle_ranges, start_idx:total_factors)
    end

    U_indices = Int[]
    U_bar_indices = Int[]
    for (i, f) in enumerate(all_factors)
        if f isa SymbolicMatrix && f.special_type == :U
            if f.is_adj
                push!(U_bar_indices, i)
            else
                push!(U_indices, i)
            end
        end
    end
    n_U = length(U_indices)
    n_U_bar = length(U_bar_indices)
    if n_U != n_U_bar
        return Num(0)
    end

    all_slots = sort([U_indices; U_bar_indices])
    if isempty(all_slots)
        return _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    end

    # Build Wires and evaluate constants in U-bearing cycles
    wires, reverse_wires = _build_wires(U_indices, U_bar_indices, cycle_ranges, all_factors)
    # Re-evaluate constant part to include cycles that have no U
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    u_map = Dict(idx => m for (m, idx) in enumerate(U_indices))
    ub_map = Dict(idx => m for (m, idx) in enumerate(U_bar_indices))
    permutations = collect(Combinatorics.permutations(1:n_U))
    total_val = Num(0)

    for sigma in permutations
        for tau in permutations
            inv_sigma = invperm(sigma)
            inv_tau = invperm(tau)
            visited = falses(total_factors, 2)
            current_term_traces = Num[]
            
            for slot in all_slots
                for port = 1:2
                    if !visited[slot, port]
                        val = _traverse_trace_cycle_final(
                            slot,
                            port,
                            sigma,
                            tau,
                            inv_sigma,
                            inv_tau,
                            wires,
                            reverse_wires,
                            u_map,
                            ub_map,
                            visited,
                            U_indices,
                            U_bar_indices,
                            all_factors,
                            dim,
                        )
                        push!(current_term_traces, val)
                    end
                end
            end

            P = [sigma[inv_tau[i]] for i = 1:n_U]
            wg_val = weingarten(get_cycle_type(P), dim)
            term_prod = isempty(current_term_traces) ? Num(1) : prod(current_term_traces)
            total_val += term_prod * wg_val
        end
    end

    return constant_part * total_val
end

function _build_wires(U_indices, U_bar_indices, cycle_ranges, all_factors)
    wires = Dict{Int,Any}()
    reverse_wires = Dict{Int,Any}()
    all_slots = sort([U_indices; U_bar_indices])

    for slot in all_slots
        # Determine which cycle this slot belongs to
        cid = findfirst(rng -> slot in rng, cycle_ranges)
        rng = cycle_ranges[cid]

        # Forward wire from index 2 to next index 1
        curr = slot
        consts = Any[]
        while true
            curr = curr == last(rng) ? first(rng) : curr + 1
            if curr in all_slots
                wires[slot] = (curr, isempty(consts) ? nothing : consts)
                break
            end
            push!(consts, all_factors[curr])
        end

        # Backward wire from index 1 to previous index 2
        curr = slot
        consts_rev = Any[]
        while true
            curr = curr == first(rng) ? last(rng) : curr - 1
            if curr in all_slots
                reverse_wires[slot] = (curr, isempty(consts_rev) ? nothing : consts_rev)
                break
            end
            push!(consts_rev, adjoint(all_factors[curr]))
        end
    end
    return wires, reverse_wires
end

function _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    res = t.prefactor
    for (cid, rng) in enumerate(cycle_ranges)
        # Check if any slot in rng is in all_slots
        has_U = false
        for idx in rng
            if idx in all_slots
                has_U = true
                break
            end
        end

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

function _traverse_trace_cycle_final(
    s,
    p,
    sigma,
    tau,
    inv_sigma,
    inv_tau,
    wires,
    reverse_wires,
    u_map,
    ub_map,
    visited,
    U_indices,
    U_bar_indices,
    all_factors,
    dim,
)
    curr_factors = Any[]
    
    while !visited[s, p]
        visited[s, p] = true
        
        # 1. Weingarten Matching
        if haskey(u_map, s)
            u_m = u_map[s]
            if p == 1 # Row-port of U matches Row-port of U_bar
                ub_k = sigma[u_m]
                s = U_bar_indices[ub_k]
                p = all_factors[s].is_adj ? 2 : 1
            else # Col-port of U matches Col-port of U_bar
                ub_k = tau[u_m]
                s = U_bar_indices[ub_k]
                p = all_factors[s].is_adj ? 1 : 2
            end
        else # s is in U_bar
            ub_m = ub_map[s]
            is_adj = all_factors[s].is_adj
            is_row_port = (is_adj && p == 2) || (!is_adj && p == 1)
            if is_row_port
                u_k = inv_sigma[ub_m]
                s = U_indices[u_k]
                p = 1
            else
                u_k = inv_tau[ub_m]
                s = U_indices[u_k]
                p = 2
            end
        end
        visited[s, p] = true
        
        # 2. Wire Traversal
        if p == 2
            s, mat_segment = wires[s]
            if mat_segment !== nothing
                append!(curr_factors, mat_segment)
            end
            p = 1
        else
            s, mat_segment = reverse_wires[s]
            if mat_segment !== nothing
                append!(curr_factors, mat_segment)
            end
            p = 2
        end
    end
    
    if isempty(curr_factors)
        return dim
    else
        return tr_val(curr_factors)
    end
end
