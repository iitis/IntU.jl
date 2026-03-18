
for (T, tag, ctor) in [
    (:GUEMeasure, :GUE, :dGUE),
    (:GOEMeasure, :GOE, :dGOE),
    (:GSEMeasure, :GSE, :dGSE),
    (:GinUEMeasure, :GinUE, :dGinUE),
    (:GinOEMeasure, :GinOE, :dGinOE),
    (:GinSEMeasure, :GinSE, :dGinSE),
]
    @eval begin
        struct $T{D,M} <: AbstractMeasure
            dim::D
            matcher::M
        end
        $T(dim) = $T(dim, nothing)

        function IntU.measure_info(measure::$T)
            subs_dict = Dict{Any,Any}()
            matcher =
                measure.matcher === nothing ? MetadataMatcher($(QuoteNode(tag))) :
                measure.matcher
            dim = measure.dim isa SymbolicMatrix ? measure.dim.dim : measure.dim
            return (subs_dict, matcher, dim, $(QuoteNode(tag)))
        end

        IntU._reconstruct_symbolic(::$T, d_asymp) = $ctor(d_asymp)
    end
end


"""
    dGUE(dim)

Gaussian Unitary Ensemble (GUE) measure.
"""
dGUE(dim) = GUEMeasure(dim)

"""
    dGOE(dim)

Gaussian Orthogonal Ensemble (GOE) measure.
"""
dGOE(dim) = GOEMeasure(dim)

"""
    dGSE(dim)

Gaussian Symplectic Ensemble (GSE) measure.
`dim` must be even.
"""
function dGSE(dim)
    if dim isa Integer && isodd(dim)
        throw(ArgumentError("GSE dimension must be even, got $dim"))
    end
    return GSEMeasure(dim)
end

"""
    dGinUE(dim)

Complex Ginibre Ensemble (GinUE) measure.
"""
dGinUE(dim) = GinUEMeasure(dim)

"""
    dGinOE(dim)

Real Ginibre Ensemble (GinOE) measure.
"""
dGinOE(dim) = GinOEMeasure(dim)

"""
    dGinSE(dim)

Quaternionic/Symplectic Ginibre Ensemble (GinSE) measure.
`dim` must be even.
"""
function dGinSE(dim)
    if dim isa Integer && isodd(dim)
        throw(ArgumentError("GinSE dimension must be even, got $dim"))
    end
    return GinSEMeasure(dim)
end


"""
    _find_tagged_indices(all_factors, tag; separate_adj=false)

Find indices of `SymbolicMatrix` entries with the given `special_type` tag.

If `separate_adj=false`, returns a single `Vector{Int}` of all matching indices.
If `separate_adj=true`, returns `(indices, bar_indices)` separating non-adjoint
and adjoint entries.
"""
function _find_tagged_indices(all_factors, tag; separate_adj = false)
    if separate_adj
        indices = Int[]
        bar_indices = Int[]
        for (i, f) in enumerate(all_factors)
            if f isa SymbolicMatrix && f.special_type == tag
                push!(f.is_adj ? bar_indices : indices, i)
            end
        end
        return indices, bar_indices
    else
        indices = Int[]
        for (i, f) in enumerate(all_factors)
            if f isa SymbolicMatrix && f.special_type == tag
                push!(indices, i)
            end
        end
        return indices
    end
end

"""
    _wick_hermitian_integrate(t::LazyTrace, dim, tag, matcher; symmetric=false)

Unified Wick contraction for Hermitian-type (GUE, GOE) ensembles.

- `symmetric=false` (GUE): E[H_ij H_kl] = δ_il δ_jk — single contraction per pair
- `symmetric=true`  (GOE): E[H_ij H_kl] = δ_il δ_jk + δ_ik δ_jl — two contractions per pair
"""
function _wick_hermitian_integrate(t::LazyTrace, dim, tag, matcher; symmetric = false)
    if isempty(t.cycles)
        return t.prefactor
    end

    H_type = (matcher isa MetadataMatcher) ? matcher.type_tag : tag

    total_factors, cycle_ranges, all_factors = _extract_trace_data(t)
    H_indices = _find_tagged_indices(all_factors, H_type)

    n_H = length(H_indices)
    if isodd(n_H)
        return 0
    end

    all_slots = sort(H_indices)
    if isempty(all_slots)
        return _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    end

    wires, reverse_wires = _build_wires(H_indices, Int[], cycle_ranges, all_factors)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    partitions = get_pair_partitions(n_H)
    total_val = 0
    h_map = Dict(idx => m for (m, idx) in enumerate(H_indices))

    for pi in partitions
        partner_map = Dict{Int,Int}()
        for (u, v) in pi
            partner_map[u] = v
            partner_map[v] = u
        end

        if symmetric
            choice_combos = Iterators.product(fill([1, 2], n_H ÷ 2)...)
        else
            choice_combos = (ntuple(_ -> 1, n_H ÷ 2),)
        end

        for choices in choice_combos
            pair_choices = Dict{Int,Int}()
            for (p_idx, (u, v)) in enumerate(pi)
                pair_choices[u] = choices[p_idx]
                pair_choices[v] = choices[p_idx]
            end

            visited = falses(total_factors, 2)
            current_partition_traces = []

            for slot in all_slots
                for port = 1:2
                    if !visited[slot, port]
                        curr_factors = Any[]
                        s, p = slot, port

                        while !visited[s, p]
                            visited[s, p] = true

                            m = h_map[s]
                            partner_m = partner_map[m]
                            choice = pair_choices[m]

                            s = H_indices[partner_m]
                            if choice == 1
                                p = (p == 1 ? 2 : 1)
                            end
                            visited[s, p] = true

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
                            push!(current_partition_traces, dim)
                        else
                            push!(current_partition_traces, tr_val(curr_factors))
                        end
                    end
                end
            end

            total_val +=
                isempty(current_partition_traces) ? 1 : prod(current_partition_traces)
        end
    end

    return constant_part * total_val
end


function fallback_integrate(t::LazyTrace, measure::GUEMeasure)
    matcher = measure.matcher === nothing ? MetadataMatcher(:GUE) : measure.matcher
    return _wick_hermitian_integrate(t, measure.dim, :GUE, matcher; symmetric = false)
end

function fallback_integrate(t::LazyTrace, measure::GOEMeasure)
    matcher = measure.matcher === nothing ? MetadataMatcher(:GOE) : measure.matcher
    return _wick_hermitian_integrate(t, measure.dim, :GOE, matcher; symmetric = true)
end

"""
    _duality_integrate(t, measure, partner_ctor, tag)

Shared duality pattern for GSE↔GOE and GinSE↔GinOE:
  1. Count tagged matrices
  2. Integrate via partner measure with MetadataMatcher
  3. Substitute dim → -dim and apply sign (-1)^(n/2+1)
"""
function _duality_integrate(t::LazyTrace, measure, partner_ctor, tag::Symbol)
    all_factors = vcat(t.cycles...)
    n = count(f -> f isa SymbolicMatrix && f.special_type == tag, all_factors)

    if isodd(n)
        return 0
    end

    partner_m = partner_ctor(measure.dim, MetadataMatcher(tag))
    partner_res = integrate(t, partner_m)

    dim = measure.dim
    res_subbed = Symbolics.substitute(partner_res, Dict(dim => -dim))
    final_sign = ((-1)^(n ÷ 2 + 1))
    return final_sign * res_subbed
end

function fallback_integrate(t::LazyTrace, measure::GSEMeasure)
    return _duality_integrate(t, measure, GOEMeasure, :GSE)
end

function fallback_integrate(t::LazyTrace, measure::GinUEMeasure)
    dim = measure.dim
    if isempty(t.cycles)
        return t.prefactor
    end

    matcher = measure.matcher === nothing ? MetadataMatcher(:GinUE) : measure.matcher
    G_type = (matcher isa MetadataMatcher) ? matcher.type_tag : :GinUE

    total_factors, cycle_ranges, all_factors = _extract_trace_data(t)
    G_indices, G_bar_indices =
        _find_tagged_indices(all_factors, G_type; separate_adj = true)

    n_G = length(G_indices)
    n_G_bar = length(G_bar_indices)

    if n_G != n_G_bar
        return 0
    end

    all_slots = sort([G_indices; G_bar_indices])
    if isempty(all_slots)
        return _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    end

    wires, reverse_wires = _build_wires(G_indices, G_bar_indices, cycle_ranges, all_factors)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    perms = permutations(1:n_G_bar)
    total_val = 0
    pos_map = Dict(idx => i for (i, idx) in enumerate(all_slots))

    for p in perms
        visited = falses(length(all_slots))
        current_partition_traces = []

        for start_pos = 1:length(all_slots)
            if !visited[start_pos]
                curr_trace_factors = Any[]
                curr_pos = start_pos
                while !visited[curr_pos]
                    visited[curr_pos] = true
                    curr_factor_idx = all_slots[curr_pos]

                    dest_factor_idx = 0
                    mat_segment = nothing

                    if haskey(wires, curr_factor_idx)
                        dest_factor_idx, mat_segment = wires[curr_factor_idx]
                    elseif haskey(reverse_wires, curr_factor_idx)
                        dest_factor_idx, mat_segment = reverse_wires[curr_factor_idx]
                    else
                        throw(ErrorException("Connectivity error in GinUE integration"))
                    end

                    if mat_segment !== nothing
                        append!(curr_trace_factors, mat_segment)
                    end

                    if dest_factor_idx in G_indices
                        m = findfirst(==(dest_factor_idx), G_indices)
                        partner_idx = G_bar_indices[p[m]]
                    else
                        m_bar = findfirst(==(dest_factor_idx), G_bar_indices)
                        m = findfirst(==(m_bar), p)
                        partner_idx = G_indices[m]
                    end

                    curr_pos = pos_map[partner_idx]
                end

                if isempty(curr_trace_factors)
                    push!(current_partition_traces, dim)
                else
                    push!(current_partition_traces, tr_val(curr_trace_factors))
                end
            end
        end
        total_val += isempty(current_partition_traces) ? 1 : prod(current_partition_traces)
    end

    return constant_part * total_val
end

function fallback_integrate(t::LazyTrace, measure::GinOEMeasure)
    dim = measure.dim
    if isempty(t.cycles)
        return t.prefactor
    end

    matcher = measure.matcher === nothing ? MetadataMatcher(:GinOE) : measure.matcher
    G_type = (matcher isa MetadataMatcher) ? matcher.type_tag : :GinOE

    total_factors, cycle_ranges, all_factors = _extract_trace_data(t)
    G_indices = _find_tagged_indices(all_factors, G_type)

    n_G = length(G_indices)
    if isodd(n_G)
        return 0
    end

    all_slots = sort(G_indices)
    if isempty(all_slots)
        return _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    end

    wires, reverse_wires = _build_wires(G_indices, Int[], cycle_ranges, all_factors)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    partitions = get_pair_partitions(n_G)
    total_val = 0
    pos_map = Dict(idx => i for (i, idx) in enumerate(G_indices))

    for pi in partitions
        visited_ports = falses(n_G, 2)
        current_partition_traces = []

        for start_m = 1:n_G
            for start_port in [1, 2]
                if !visited_ports[start_m, start_port]
                    curr_trace_factors = Any[]
                    curr_m = start_m
                    curr_port = start_port

                    while !visited_ports[curr_m, curr_port]
                        visited_ports[curr_m, curr_port] = true

                        curr_factor_idx = G_indices[curr_m]

                        dest_factor_idx = 0
                        mat_segment = nothing
                        if curr_port == 2
                            dest_factor_idx, mat_segment = wires[curr_factor_idx]
                            landed_port = 1
                        else
                            dest_factor_idx, mat_segment = reverse_wires[curr_factor_idx]
                            landed_port = 2
                        end

                        if mat_segment !== nothing
                            append!(curr_trace_factors, mat_segment)
                        end

                        landed_m = pos_map[dest_factor_idx]
                        f_landed = all_factors[dest_factor_idx]

                        visited_ports[landed_m, landed_port] = true

                        partner_m = 0
                        for (u, v) in pi
                            if u == landed_m
                                partner_m = v
                                break
                            elseif v == landed_m
                                partner_m = u
                                break
                            end
                        end

                        f_partner = all_factors[G_indices[partner_m]]
                        if f_landed.is_trans == f_partner.is_trans
                            partner_port = landed_port
                        else
                            partner_port = (landed_port == 1 ? 2 : 1)
                        end

                        curr_m = partner_m
                        curr_port = partner_port
                    end

                    if isempty(curr_trace_factors)
                        push!(current_partition_traces, dim)
                    else
                        push!(current_partition_traces, tr_val(curr_trace_factors))
                    end
                end
            end
        end
        total_val += isempty(current_partition_traces) ? 1 : prod(current_partition_traces)
    end
    return constant_part * total_val
end

function fallback_integrate(t::LazyTrace, measure::GinSEMeasure)
    return _duality_integrate(t, measure, GinOEMeasure, :GinSE)
end
