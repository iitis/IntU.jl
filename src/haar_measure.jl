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
    dU(U::SymbolicMatrix)

Defines the Haar measure for the Unitary group $U(d)$.

Can be called with just the dimension $d$, in which case the integration engine will look for 
symbolic entries tagged with `:U` via `SymbolicMatrix(:U, :U)` or metadata mapping.

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
dU(dim) = HaarMeasure(dim)
dU(U::SymbolicMatrix) = HaarMeasure(U.dim)


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

"""
    integrate(t::LazyTrace, measure::HaarMeasure)

Integrate a product of traces of matrices over the Haar measure.
Uses the graphical Weingarten calculus.
"""
function fallback_integrate(t::LazyTrace, measure::HaarMeasure)
    dim = measure.dim

    if isempty(t.cycles)
        return t.prefactor
    end

    # Identify U and U_dag instances across ALL cycles
    U_indices = Int[]
    U_bar_indices = Int[]

    total_factors = 0
    cycle_ranges = UnitRange{Int}[]
    all_factors = SymbolicMatrix[]

    for cycle in t.cycles
        start_idx = total_factors + 1
        append!(all_factors, cycle)
        total_factors += length(cycle)
        push!(cycle_ranges, start_idx:total_factors)
    end

    for (i, f) in enumerate(all_factors)
        if f.special_type == :U
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
        return 0
    end

    # Handle case with no U/U_bar (all constant traces)
    all_slots = sort([U_indices; U_bar_indices])
    if isempty(all_slots)
        return _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    end

    # Build Wires
    wires = _build_wires(U_indices, U_bar_indices, cycle_ranges, all_factors)

    # Calculate constant part (cycles with no U/U_bar)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    # Weingarten Sum
    u_map = Dict{Int,Int}(idx => m for (m, idx) in enumerate(U_indices))
    ub_map = Dict{Int,Int}(idx => m for (m, idx) in enumerate(U_bar_indices))

    permutations = collect(Combinatorics.permutations(1:n_U))
    total_val = Num(0)

    for sigma in permutations
        for tau in permutations

            inv_tau = invperm(tau)
            P = [sigma[inv_tau[i]] for i = 1:n_U]
            cycle_type = get_cycle_type(P)
            wg_val = weingarten(cycle_type, dim)

            if _symbolic_isequal(wg_val, 0)
                continue
            end

            visited_U = falses(n_U)
            visited_Ub = falses(n_U)
            current_term_traces = Num[]

            # Check cycles (start from U nodes)
            # Since graph is bipartite (U <-> Ub), any cycle contains at least one U.
            for start_m = 1:n_U
                if !visited_U[start_m]
                    val = _traverse_trace_cycle(
                        start_m,
                        1,
                        sigma,
                        inv_tau,
                        wires,
                        u_map,
                        ub_map,
                        visited_U,
                        visited_Ub,
                        U_indices,
                        U_bar_indices,
                        dim,
                    )
                    push!(current_term_traces, val)
                end
            end

            term_prod = isempty(current_term_traces) ? Num(1) : prod(current_term_traces)
            total_val += term_prod * wg_val
        end
    end

    return constant_part * total_val
end

function _build_wires(U_indices, U_bar_indices, cycle_ranges, all_factors)
    wires = Dict{Int,Any}()
    all_slots = sort([U_indices; U_bar_indices])
    n_slots = length(all_slots)

    for k = 1:n_slots
        start_idx = all_slots[k]

        # Determine which cycle this slot belongs to
        cycle_id = 0
        for (cid, rng) in enumerate(cycle_ranges)
            if start_idx in rng
                cycle_id = cid
                break
            end
        end
        cycle_range = cycle_ranges[cycle_id]

        # Traverse forward from start_idx until we hit another U-slot OR wrap around
        consts = SymbolicMatrix[]
        curr = start_idx

        while true
            # Move to next in cycle
            if curr == last(cycle_range)
                curr = first(cycle_range)
            else
                curr += 1
            end

            # Check if we hit a U/U_bar
            if curr in all_slots
                # Destination found
                end_idx = curr
                wires[start_idx] = (end_idx, isempty(consts) ? nothing : consts)
                break
            end

            # Otherwise it's a constant
            push!(consts, all_factors[curr])

            # Safety break if cycle is full of constants (shouldn't happen as we started from a U slot)
            if curr == start_idx
                error("Cycle should contain at least one U/U_bar")
            end
        end
    end
    return wires
end

function _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    constant_part = t.prefactor
    for (cid, rng) in enumerate(cycle_ranges)
        # Check if any slot in rng is in all_slots
        has_U = false
        for idx in rng
            if idx in all_slots
                has_U = true;
                break
            end
        end

        if !has_U
            cycle = t.cycles[cid]
            if isempty(cycle)
                constant_part *= dim
            else
                constant_part *= tr_val(cycle)
            end
        end
    end
    return constant_part
end

function _traverse_trace_cycle(
    start_m,
    start_type,
    sigma,
    inv_tau,
    wires,
    u_map,
    ub_map,
    visited_U,
    visited_Ub,
    U_indices,
    U_bar_indices,
    dim,
)
    curr_trace_factors = SymbolicMatrix[]
    curr_type = start_type
    curr_idx = start_m

    while true
        if curr_type == 1
            if visited_U[curr_idx]
                break
            end
            visited_U[curr_idx] = true
            next_ub_m = sigma[curr_idx]
            # visited_Ub[next_ub_m] = true # Don't mark next node yet
            start_factor_idx = U_bar_indices[next_ub_m]
            dest_factor_idx, mat_segment = wires[start_factor_idx]
            if mat_segment !== nothing
                append!(curr_trace_factors, mat_segment)
            end
            if haskey(u_map, dest_factor_idx)
                curr_type = 1;
                curr_idx = u_map[dest_factor_idx]
            else
                curr_type = 2;
                curr_idx = ub_map[dest_factor_idx]
            end
        else
            if visited_Ub[curr_idx]
                break
            end
            visited_Ub[curr_idx] = true
            next_u_m = inv_tau[curr_idx]
            # visited_U[next_u_m] = true # Don't mark next node yet
            start_factor_idx = U_indices[next_u_m]
            dest_factor_idx, mat_segment = wires[start_factor_idx]
            if mat_segment !== nothing
                append!(curr_trace_factors, mat_segment)
            end
            if haskey(u_map, dest_factor_idx)
                curr_type = 1;
                curr_idx = u_map[dest_factor_idx]
            else
                curr_type = 2;
                curr_idx = ub_map[dest_factor_idx]
            end
        end
    end

    if isempty(curr_trace_factors)
        return dim
    else
        return tr_val(curr_trace_factors)
    end
end
