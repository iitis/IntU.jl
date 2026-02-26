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
    if Symbolics.iscall(s) &&
       (Symbolics.operation(s) == conj || Symbolics.operation(s) == Base.conj)
        is_conj = true
        s = Symbolics.arguments(s)[1]
    end

    s_str = string(s)
    mat = match(m.regex, s_str)
    if mat !== nothing
        i = tryparse(Int, mat[1])
        j = tryparse(Int, mat[2])
        if i !== nothing && j !== nothing
            final_tag = is_conj ? :U_bar : :U
            return (final_tag, i, j)
        end
    end
    return nothing
end


# Dummy type to represent the measure
struct HaarMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end

# Constructor for backward compatibility
HaarMeasure(dim) = HaarMeasure(dim, nothing)
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
    matcher = measure.matcher === nothing ? MetadataMatcher(:U) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, :U)
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
    if d isa SymbolicMatrix
        d = d.dim
    end

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

    matcher = measure.matcher === nothing ? MetadataMatcher(:U) : measure.matcher
    U_type = (matcher isa MetadataMatcher) ? matcher.type_tag : :U

    U_indices = Int[]
    U_bar_indices = Int[]
    for (i, f) in enumerate(all_factors)
        if f isa SymbolicMatrix && f.special_type == U_type
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
    
    # NEW: Convert maps to speed up inner loop. Max slot index is total_factors.
    # total_factors = n_U + n_U_bar (usually 2*n_U)
    u_map_vec = fill(0, total_factors)
    for (idx, m) in u_map; u_map_vec[idx] = m; end
    ub_map_vec = fill(0, total_factors)
    for (idx, m) in ub_map; ub_map_vec[idx] = m; end
    
    # Convert wires to vectors
    wires_s = fill(0, total_factors)
    wires_m = Vector{Any}(nothing, total_factors)
    rev_wires_s = fill(0, total_factors)
    rev_wires_m = Vector{Any}(nothing, total_factors)
    for s in all_slots
        ws, wm = wires[s]
        wires_s[s] = ws; wires_m[s] = wm
        rs, rm = reverse_wires[s]
        rev_wires_s[s] = rs; rev_wires_m[s] = rm
    end

    is_trans_vec = [f.is_trans for f in all_factors]

    perms = collect(Combinatorics.permutations(1:n_U))
    inv_perms = [invperm(p) for p in perms]
    n_perms = length(perms)

    d_un = Symbolics.unwrap(dim)
    is_numeric_dim = d_un isa Integer

    # Detect "pure trace" case: no constant matrices
    is_pure_trace = all(m -> m === nothing, wires_m) && all(m -> m === nothing, rev_wires_m)
    
    if is_pure_trace
        total_val_pi = is_numeric_dim ? zero(Rational{BigInt}) : Num(0)
        p_fast = Progress(n_perms; dt=0.5, desc="Integrating Haar (n=$n_U, fast path)... ")
        for Pi in perms
            ct = get_cycle_type(Pi)
            n_cycles = length(ct)
            wg_val = weingarten(ct, dim)
            term_prod = dim^n_cycles
            total_val_pi += term_prod * wg_val
            next!(p_fast)
        end
        return constant_part * n_perms * total_val_pi
    end

    total_iters = n_perms * n_perms
    desc = "Integrating Haar (n=$n_U)... "
    p = Progress(total_iters; dt=0.5, desc=desc)

    total_val = is_numeric_dim ? zero(Rational{BigInt}) : Num(0)
    P = Vector{Int}(undef, n_U)
    U_idx_vec = collect(U_indices)
    Ub_idx_vec = collect(U_bar_indices)

    for (i, sigma) in enumerate(perms)
        inv_sigma = inv_perms[i]
        for (j, tau) in enumerate(perms)
            inv_tau = inv_perms[j]
            # Use bitmask for visited if total_factors <= 32
            # 2 ports * total_factors
            visited = UInt64(0)
            term_prod = is_numeric_dim ? one(Rational{BigInt}) : Num(1)

            for start_slot in all_slots
                for start_port = 1:2
                    # bit index: (port-1)*total_factors + slot
                    bit_idx = (start_port - 1) * total_factors + start_slot
                    if (visited & (UInt64(1) << (bit_idx - 1))) == 0
                        # Traverse cycle
                        val, visited = _traverse_trace_cycle_fast(
                            start_slot,
                            start_port,
                            sigma,
                            tau,
                            inv_sigma,
                            inv_tau,
                            wires_s,
                            wires_m,
                            rev_wires_s,
                            rev_wires_m,
                            u_map_vec,
                            ub_map_vec,
                            is_trans_vec,
                            U_idx_vec,
                            Ub_idx_vec,
                            total_factors,
                            dim,
                            visited
                        )
                        term_prod *= val
                    end
                end
            end

            for k = 1:n_U
                P[k] = sigma[inv_tau[k]]
            end
            wg_val = weingarten(get_cycle_type(P), dim)
            total_val += term_prod * wg_val
            next!(p)
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
            push!(consts_rev, transpose(all_factors[curr]))
        end
    end
    return wires, reverse_wires
end

function _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    res = t.prefactor
    slot_set = Set{Int}(all_slots)
    for (cid, rng) in enumerate(cycle_ranges)
        # Check if any slot in rng is in all_slots (O(1) set lookup)
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

function _traverse_trace_cycle_fast(
    s,
    p,
    sigma,
    tau,
    inv_sigma,
    inv_tau,
    wires_s,
    wires_m,
    rev_wires_s,
    rev_wires_m,
    u_map_vec,
    ub_map_vec,
    is_trans_vec,
    U_idx_vec,
    Ub_idx_vec,
    total_factors,
    dim,
    visited
)
    curr_factors = Any[]

    while true
        bit_idx = (p - 1) * total_factors + s
        if (visited & (UInt64(1) << (bit_idx - 1))) != 0
            break
        end
        visited |= (UInt64(1) << (bit_idx - 1))

        # 1. Weingarten Matching
        u_m = u_map_vec[s]
        if u_m != 0
            if p == 1 # Row-port of U matches Row-port of U_bar
                ub_k = sigma[u_m]
                s = Ub_idx_vec[ub_k]
                p = is_trans_vec[s] ? 2 : 1
            else # Col-port of U matches Col-port of U_bar
                ub_k = tau[u_m]
                s = Ub_idx_vec[ub_k]
                p = is_trans_vec[s] ? 1 : 2
            end
        else # s is in U_bar
            ub_m = ub_map_vec[s]
            is_row_port = (is_trans_vec[s] && p == 2) || (!is_trans_vec[s] && p == 1)
            if is_row_port
                u_k = inv_sigma[ub_m]
                s = U_idx_vec[u_k]
                p = 1
            else
                u_k = inv_tau[ub_m]
                s = U_idx_vec[u_k]
                p = 2
            end
        end
        
        # Mark the other port of the newly reached factor as visited too
        bit_idx_other = (p - 1) * total_factors + s
        visited |= (UInt64(1) << (bit_idx_other - 1))

        # 2. Wire Traversal
        if p == 2
            mat_segment = wires_m[s]
            if mat_segment !== nothing
                append!(curr_factors, mat_segment)
            end
            s = wires_s[s]
            p = 1
        else
            mat_segment = rev_wires_m[s]
            if mat_segment !== nothing
                append!(curr_factors, mat_segment)
            end
            s = rev_wires_s[s]
            p = 2
        end
    end

    if isempty(curr_factors)
        return dim, visited
    else
        return tr_val(curr_factors), visited
    end
end
