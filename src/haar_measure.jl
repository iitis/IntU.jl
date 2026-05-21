
struct SymbolicMatcher <: AbstractIndexMatcher
    tag::Symbol
    regex::Regex
end

function match_index(m::SymbolicMatcher, t)
    s = Symbolics.unwrap(t)
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


struct HaarMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end

HaarMeasure(dim) = HaarMeasure(dim, nothing)
@doc raw"""
    dU(dim)

Defines the Haar measure for the Unitary group $U(d)$.

Integration engine identifies variables via metadata tag `:U`.
"""
function dU(dim)
    _assert_no_float_param(dim, "dim", "dU")
    return HaarMeasure(dim)
end


IntU._measure_tag(::HaarMeasure) = :U

function _manual_fallback(expr, measure::HaarMeasure)
    throw(ArgumentError("HaarMeasure integration failed for: $(typeof(expr))"))
end

IntU._reconstruct_symbolic(::HaarMeasure, d_asymp) = dU(d_asymp)


function fallback_integrate(t::LazyTrace, measure::HaarMeasure)
    dim = measure.dim
    constant_part = t.prefactor

    total_factors, cycle_ranges, all_factors = _extract_trace_data(t)

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

    wires, reverse_wires = _build_wires(U_indices, U_bar_indices, cycle_ranges, all_factors)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    u_map = Dict(idx => m for (m, idx) in enumerate(U_indices))
    ub_map = Dict(idx => m for (m, idx) in enumerate(U_bar_indices))

    u_map_vec = fill(0, total_factors)
    for (idx, m) in u_map
        ;
        u_map_vec[idx] = m;
    end
    ub_map_vec = fill(0, total_factors)
    for (idx, m) in ub_map
        ;
        ub_map_vec[idx] = m;
    end

    wires_s = fill(0, total_factors)
    wires_m = Vector{Any}(nothing, total_factors)
    rev_wires_s = fill(0, total_factors)
    rev_wires_m = Vector{Any}(nothing, total_factors)
    for s in all_slots
        ws, wm = wires[s]
        wires_s[s] = ws;
        wires_m[s] = wm
        rs, rm = reverse_wires[s]
        rev_wires_s[s] = rs;
        rev_wires_m[s] = rm
    end

    is_trans_vec = [f.is_trans for f in all_factors]

    perms = collect(Combinatorics.permutations(1:n_U))
    inv_perms = [invperm(p) for p in perms]
    n_perms = length(perms)

    d_un = Symbolics.unwrap(dim)
    is_numeric_dim = d_un isa Integer

    is_pure_trace = all(m -> m === nothing, wires_m) && all(m -> m === nothing, rev_wires_m)

    if is_pure_trace
        total_val_pi = is_numeric_dim ? zero(Rational{BigInt}) : Num(0)
        p_fast =
            Progress(n_perms; dt = 0.5, desc = "Integrating Haar (n=$n_U, fast path)... ")
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
    p = Progress(total_iters; dt = 0.5, desc = desc)

    total_val = is_numeric_dim ? zero(Rational{BigInt}) : Num(0)
    P = Vector{Int}(undef, n_U)
    U_idx_vec = collect(U_indices)
    Ub_idx_vec = collect(U_bar_indices)
    curr_factors_buf = Any[]

    for (i, sigma) in enumerate(perms)
        inv_sigma = inv_perms[i]
        for (j, tau) in enumerate(perms)
            inv_tau = inv_perms[j]
            visited = falses(2 * total_factors)
            term_prod = is_numeric_dim ? one(Rational{BigInt}) : Num(1)

            for start_slot in all_slots
                for start_port = 1:2
                    bit_idx = (start_port - 1) * total_factors + start_slot
                    if !visited[bit_idx]
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
                            visited,
                            curr_factors_buf,
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
    visited,
    curr_factors = Any[],
)
    empty!(curr_factors)

    while true
        bit_idx = (p - 1) * total_factors + s
        if visited[bit_idx]
            break
        end
        visited[bit_idx] = true

        u_m = u_map_vec[s]
        if u_m != 0
            if p == 1
                ub_k = sigma[u_m]
                s = Ub_idx_vec[ub_k]
                p = is_trans_vec[s] ? 2 : 1
            else
                ub_k = tau[u_m]
                s = Ub_idx_vec[ub_k]
                p = is_trans_vec[s] ? 1 : 2
            end
        else
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

        bit_idx_other = (p - 1) * total_factors + s
        visited[bit_idx_other] = true

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
