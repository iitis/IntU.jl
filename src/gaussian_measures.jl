# Gaussian Random Matrix measures (GUE, GOE, GSE)

struct GUEMeasure{D,M}
    dim::D
    matcher::M
end
GUEMeasure(dim) = GUEMeasure(dim, nothing)

struct GOEMeasure{D,M}
    dim::D
    matcher::M
end
GOEMeasure(dim) = GOEMeasure(dim, nothing)

struct GSEMeasure{D,M}
    dim::D
    matcher::M
end
GSEMeasure(dim) = GSEMeasure(dim, nothing)

struct GinUEMeasure{D,M}
    dim::D
    matcher::M
end
GinUEMeasure(dim) = GinUEMeasure(dim, nothing)

struct GinOEMeasure{D,M}
    dim::D
    matcher::M
end
GinOEMeasure(dim) = GinOEMeasure(dim, nothing)

struct GinSEMeasure{D,M}
    dim::D
    matcher::M
end
GinSEMeasure(dim) = GinSEMeasure(dim, nothing)

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
    integrate(expr, measure::GUEMeasure)
"""
# Generate element-wise integration for all Gaussian measure types
for T_measure in
    [GUEMeasure, GOEMeasure, GSEMeasure, GinUEMeasure, GinOEMeasure, GinSEMeasure]
    @eval function integrate(expr::AbstractArray, measure::$T_measure)
        return map(e -> integrate(e, measure), expr)
    end
    @eval function integrate(expr::SymbolicMatrix, measure::$T_measure)
        return invoke(integrate, Tuple{SymbolicMatrix,Any}, expr, measure)
    end
    @eval function integrate(expr::SymbolicMatrixProduct, measure::$T_measure)
        return invoke(integrate, Tuple{SymbolicMatrixProduct,Any}, expr, measure)
    end
end

# Generate measure_info methods for all Gaussian measure types
for (T_measure, tag) in [
    (GUEMeasure, :GUE), (GOEMeasure, :GOE), (GSEMeasure, :GSE),
    (GinUEMeasure, :GinUE), (GinOEMeasure, :GinOE), (GinSEMeasure, :GinSE),
]
    @eval function IntU.measure_info(measure::$T_measure)
        subs_dict = Dict{Any,Any}()
        matcher = measure.matcher === nothing ? MetadataMatcher($(QuoteNode(tag))) : measure.matcher
        dim = measure.dim isa SymbolicMatrix ? measure.dim.dim : measure.dim
        return (subs_dict, matcher, dim, $(QuoteNode(tag)))
    end
end

function fallback_integrate(t::LazyTrace, measure::GUEMeasure)
    factors = t.factors
    # We look for matrices tagged correctly
    matcher = measure.matcher === nothing ? MetadataMatcher(:GUE) : measure.matcher
    H_type = (matcher isa MetadataMatcher) ? matcher.type_tag : :GUE


    H_indices = Int[]
    for (i, f) in enumerate(factors)
        if f isa SymbolicMatrix && f.special_type == H_type
            push!(H_indices, i)
        end
    end

    n_H = length(H_indices)
    if isodd(n_H)
        return 0
    end
    if n_H == 0
        return tr_val(factors)
    end

    dim = measure.dim
    n_factors = length(factors)


    wires = Dict{Int,Any}()
    for k = 1:n_H
        idx = H_indices[k]
        next_h_idx = H_indices[mod1(k+1, n_H)]

        consts = SymbolicMatrix[]
        curr = mod1(idx + 1, n_factors)
        while curr != next_h_idx
            push!(consts, factors[curr])
            curr = mod1(curr + 1, n_factors)
        end
        wires[idx] = (next_h_idx, isempty(consts) ? nothing : consts)
    end

    partitions = get_pair_partitions(n_H)
    total_val = 0

    for pi in partitions
        # Each pairing connects trace cycles through Wick contraction

        # Map pairing to permutation
        perm_map = Dict{Int,Int}()
        for (u, v) in pi
            perm_map[u] = v
            perm_map[v] = u
        end


        visited = falses(n_H)
        current_partition_traces = []

        for start_m = 1:n_H
            if !visited[start_m]
                curr_trace_factors = SymbolicMatrix[]
                curr_m = start_m
                while !visited[curr_m]
                    visited[curr_m] = true
                    paired_m = perm_map[curr_m]
                    dest_factor_idx, mat_segment = wires[H_indices[curr_m]]
                    if mat_segment !== nothing
                        append!(curr_trace_factors, mat_segment)
                    end

                    # Find landing slot's m-index
                    next_m = 1
                    while H_indices[next_m] != dest_factor_idx
                        next_m += 1
                    end

                    # Wick contraction jump to partner
                    curr_m = perm_map[next_m]
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

    return total_val
end



function fallback_integrate(t::LazyTrace, measure::GOEMeasure)
    factors = t.factors
    matcher = measure.matcher === nothing ? MetadataMatcher(:GOE) : measure.matcher
    H_type = (matcher isa MetadataMatcher) ? matcher.type_tag : :GOE

    H_indices = Int[]
    for (i, f) in enumerate(factors)
        if f isa SymbolicMatrix && f.special_type == H_type
            push!(H_indices, i)
        end
    end

    n_H = length(H_indices)
    if isodd(n_H)
        return 0
    end
    if n_H == 0
        return tr_val(factors)
    end

    dim = measure.dim
    n_factors = length(factors)


    wires = Dict{Int,Any}()
    for k = 1:n_H
        idx = H_indices[k]
        next_h_idx = H_indices[mod1(k+1, n_H)]
        prev_h_idx = H_indices[mod1(k-1, n_H)]


        fwd_consts = SymbolicMatrix[]
        curr = mod1(idx + 1, n_factors)
        while curr != next_h_idx
            push!(fwd_consts, factors[curr])
            curr = mod1(curr + 1, n_factors)
        end

        # Backward segment (reversed with adjoint)
        bwd_consts = SymbolicMatrix[]
        curr = mod1(idx - 1, n_factors)
        while curr != prev_h_idx
            push!(bwd_consts, factors[curr]') # Adjoint because we go backwards
            curr = mod1(curr - 1, n_factors)
        end

        wires[idx] = (
            next_h_idx,
            isempty(fwd_consts) ? nothing : fwd_consts,
            prev_h_idx,
            isempty(bwd_consts) ? nothing : bwd_consts,
        )
    end

    partitions = get_pair_partitions(n_H)
    total_val = 0

    for pi in partitions
        # GOE: <H_ij H_kl> = delta_ik delta_jl + delta_il delta_jk
        # Sum over all 2^(n_H/2) contraction-type choices

        # Pre-build partner lookup: m -> (partner_m, pair_idx)
        partner_lookup = Dict{Int,Tuple{Int,Int}}()
        for (p_idx, (u, v)) in enumerate(pi)
            partner_lookup[u] = (v, p_idx)
            partner_lookup[v] = (u, p_idx)
        end

        choice_combinations = collect(Iterators.product(fill([1, 2], n_H ÷ 2)...))

        for choices in choice_combinations

            visited = falses(n_H)
            current_partition_traces = []

            visited_ports = falses(n_H, 2)

            for start_m = 1:n_H
                for start_port in [1, 2]
                    if !visited_ports[start_m, start_port]
                        curr_trace_factors = SymbolicMatrix[]
                        curr_m = start_m
                        curr_port = start_port

                        while !visited_ports[curr_m, curr_port]
                            visited_ports[curr_m, curr_port] = true

                            if curr_port == 2
                                # Exit Port 2 (Output), take forward wire
                                dest_idx, fwd, prev_idx, bwd = wires[H_indices[curr_m]]
                                if fwd !== nothing
                                    append!(curr_trace_factors, fwd)
                                end

                                # Land at Port 1 (Input) of some H
                                landed_m = 1
                                while H_indices[landed_m] != dest_idx
                                    landed_m += 1
                                end
                                visited_ports[landed_m, 1] = true

                                # Use Wick contraction jump
                                partner_m, pair_idx = partner_lookup[landed_m]

                                if choices[pair_idx] == 2 # delta_il delta_jk (P1 -> P2)
                                    curr_m = partner_m
                                    curr_port = 2
                                else # delta_ik delta_jl (P1 -> P1)
                                    curr_m = partner_m
                                    curr_port = 1
                                end
                            else
                                # Exit Port 1 (Input), take backward wire
                                dest_idx, fwd, prev_idx, bwd = wires[H_indices[curr_m]]
                                if bwd !== nothing
                                    append!(curr_trace_factors, bwd)
                                end

                                # Land at Port 2 (Output) of some H
                                landed_m = 1
                                while H_indices[landed_m] != prev_idx
                                    landed_m += 1
                                end
                                visited_ports[landed_m, 2] = true

                                partner_m, pair_idx = partner_lookup[landed_m]

                                if choices[pair_idx] == 2 # delta_il delta_jk (P2 -> P1)
                                    curr_m = partner_m
                                    curr_port = 1
                                else # delta_ik delta_jl (P2 -> P2)
                                    curr_m = partner_m
                                    curr_port = 2
                                end
                            end
                        end

                        if isempty(curr_trace_factors)
                            push!(current_partition_traces, dim)
                        else
                            push!(current_partition_traces, tr_val(curr_trace_factors))
                        end
                    end
                end
            end



            term_val =
                isempty(current_partition_traces) ? 1 : prod(current_partition_traces)
            total_val += term_val
        end
    end
    return total_val
end


function fallback_integrate(t::LazyTrace, measure::GSEMeasure)
    # GSE-GOE duality: <Tr(H^k)>_GSE(d) = (-1)^(k/2+1) <Tr(H^k)>_GOE(-d)
    factors = t.factors
    H_type = :GSE
    n_H = count(f -> f isa SymbolicMatrix && f.special_type == H_type, factors)

    if isodd(n_H)
        return 0
    end

    goe_m = GOEMeasure(measure.dim, MetadataMatcher(H_type))
    goe_res = integrate(t, goe_m)

    # Substitute d -> -d
    dim = measure.dim
    res_subbed = Symbolics.substitute(goe_res, Dict(dim => -dim))


    final_sign = ((-1)^(n_H ÷ 2 + 1))
    return final_sign * res_subbed
end

function fallback_integrate(t::LazyTrace, measure::GinUEMeasure)
    dim = measure.dim
    if isempty(t.cycles)
        return t.prefactor
    end

    matcher = measure.matcher === nothing ? MetadataMatcher(:GinUE) : measure.matcher
    G_type = (matcher isa MetadataMatcher) ? matcher.type_tag : :GinUE


    all_factors = Any[]
    cycle_ranges = UnitRange{Int}[]
    total_factors = 0
    for cycle in t.cycles
        start_idx = total_factors + 1
        append!(all_factors, cycle)
        total_factors += length(cycle)
        push!(cycle_ranges, start_idx:total_factors)
    end


    G_indices = Int[]
    G_bar_indices = Int[]
    for (i, f) in enumerate(all_factors)
        if f isa SymbolicMatrix && f.special_type == G_type
            if f.is_adj
                push!(G_bar_indices, i)
            else
                push!(G_indices, i)
            end
        end
    end

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
                        error("Connectivity error in GinUE integration")
                    end

                    if mat_segment !== nothing
                        append!(curr_trace_factors, mat_segment)
                    end

                    # Find partner via permutation
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

    all_factors = Any[]
    cycle_ranges = UnitRange{Int}[]
    total_factors = 0
    for cycle in t.cycles
        start_idx = total_factors + 1
        append!(all_factors, cycle)
        total_factors += length(cycle)
        push!(cycle_ranges, start_idx:total_factors)
    end

    G_indices = Int[]
    for (i, f) in enumerate(all_factors)
        if f isa SymbolicMatrix && f.special_type == G_type
            push!(G_indices, i)
        end
    end

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
        # Track which ports of which G matrices we have visited
        # visited_ports[m_index, port_index] where port_index 1=Input, 2=Output
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

                        # Exit current node via curr_port
                        dest_factor_idx = 0
                        mat_segment = nothing
                        if curr_port == 2
                            # From Port 2 (Output), follow forward wire to next Port 1
                            dest_factor_idx, mat_segment = wires[curr_factor_idx]
                            landed_port = 1
                        else
                            # From Port 1 (Input), follow backward wire to previous Port 2
                            dest_factor_idx, mat_segment = reverse_wires[curr_factor_idx]
                            landed_port = 2
                        end

                        if mat_segment !== nothing
                            append!(curr_trace_factors, mat_segment)
                        end

                        landed_m = pos_map[dest_factor_idx]
                        f_landed = all_factors[dest_factor_idx]

                        visited_ports[landed_m, landed_port] = true

                        # Jump to partner via Wick contraction.
                        # Contraction matches indices: row matches row, col matches col.
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

                        # Determine partner port via index matching
                        f_partner = all_factors[G_indices[partner_m]]
                        # If same is_adj, row is at same port index. If different, they swap.
                        if f_landed.is_adj == f_partner.is_adj
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
    G_type = :GinSE
    all_factors = vcat(t.cycles...)
    n_G = count(f -> f isa SymbolicMatrix && f.special_type == G_type, all_factors)

    if isodd(n_G)
        return 0
    end

    ginoe_m = GinOEMeasure(measure.dim, MetadataMatcher(G_type))
    ginoe_res = integrate(t, ginoe_m)
    dim = measure.dim
    res_subbed = Symbolics.substitute(ginoe_res, Dict(dim => -dim))
    final_sign = ((-1)^(n_G ÷ 2 + 1))
    return final_sign * res_subbed
end

# Generate asymptotic methods for all Gaussian measure types
for (T_measure, d_ctor) in [
    (GUEMeasure, dGUE), (GOEMeasure, dGOE), (GSEMeasure, dGSE),
    (GinUEMeasure, dGinUE), (GinOEMeasure, dGinOE), (GinSEMeasure, dGinSE),
]
    @eval function asymptotic(expr, measure::$T_measure, order = 1)
        d = measure.dim
        if d isa Symbolics.Num || !(d isa Integer)
            exact_res = integrate(expr, measure)
            return _expand_asymptotic(exact_res, d, order)
        end

        d_asymp = Symbolics.variable(:d_asymp)
        m_sym = $d_ctor(d_asymp)
        exact_res = integrate(expr, m_sym)
        return _expand_asymptotic(exact_res, d_asymp, order)
    end
end
