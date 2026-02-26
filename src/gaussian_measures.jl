# Gaussian Random Matrix measures (GUE, GOE, GSE)

struct GUEMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end
GUEMeasure(dim) = GUEMeasure(dim, nothing)

struct GOEMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end
GOEMeasure(dim) = GOEMeasure(dim, nothing)

struct GSEMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end
GSEMeasure(dim) = GSEMeasure(dim, nothing)

struct GinUEMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end
GinUEMeasure(dim) = GinUEMeasure(dim, nothing)

struct GinOEMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end
GinOEMeasure(dim) = GinOEMeasure(dim, nothing)

struct GinSEMeasure{D,M} <: AbstractMeasure
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

# Integration rules for each ensemble

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
    dim = measure.dim
    if isempty(t.cycles)
        return t.prefactor
    end

    matcher = measure.matcher === nothing ? MetadataMatcher(:GUE) : measure.matcher
    H_type = (matcher isa MetadataMatcher) ? matcher.type_tag : :GUE

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

    H_indices = Int[]
    for (i, f) in enumerate(all_factors)
        if f isa SymbolicMatrix && f.special_type == H_type
            push!(H_indices, i)
        end
    end

    n_H = length(H_indices)
    if isodd(n_H)
        return 0
    end
    
    all_slots = sort(H_indices)
    if isempty(all_slots)
        return _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    end

    # Build Wires and evaluate constants
    wires, reverse_wires = _build_wires(H_indices, Int[], cycle_ranges, all_factors)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    partitions = get_pair_partitions(n_H)
    total_val = 0

    h_map = Dict(idx => m for (m, idx) in enumerate(H_indices))

    for pi in partitions
        # GUE Contraction: E[H_ij H_kl] = delta_il delta_jk
        # This means Output of H_u connects to Input of H_v and vice-versa.
        
        partner_map = Dict{Int,Int}()
        for (u, v) in pi
            partner_map[u] = v
            partner_map[v] = u
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

                        # 1. Wick Jump
                        # In GUE, Output (2) always goes to Input (1) of partner
                        # and Input (1) always goes to Output (2) of partner.
                        m = h_map[s]
                        partner_m = partner_map[m]
                        s = H_indices[partner_m]
                        p = (p == 1 ? 2 : 1)
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
                        push!(current_partition_traces, dim)
                    else
                        push!(current_partition_traces, tr_val(curr_factors))
                    end
                end
            end
        end

        total_val += isempty(current_partition_traces) ? 1 : prod(current_partition_traces)
    end

    return constant_part * total_val
end



function fallback_integrate(t::LazyTrace, measure::GOEMeasure)
    dim = measure.dim
    if isempty(t.cycles)
        return t.prefactor
    end

    matcher = measure.matcher === nothing ? MetadataMatcher(:GOE) : measure.matcher
    H_type = (matcher isa MetadataMatcher) ? matcher.type_tag : :GOE

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

    H_indices = Int[]
    for (i, f) in enumerate(all_factors)
        if f isa SymbolicMatrix && f.special_type == H_type
            push!(H_indices, i)
        end
    end

    n_H = length(H_indices)
    if isodd(n_H)
        return 0
    end
    
    all_slots = sort(H_indices)
    if isempty(all_slots)
        return _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    end

    # Build Wires and evaluate constants
    wires, reverse_wires = _build_wires(H_indices, Int[], cycle_ranges, all_factors)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    partitions = get_pair_partitions(n_H)
    total_val = 0

    h_map = Dict(idx => m for (m, idx) in enumerate(H_indices))

    for pi in partitions
        # GOE Contraction: E[H_ij H_kl] = delta_il delta_jk + delta_ik delta_jl
        # 2nd term swaps ports for the partner comparison
        
        partner_lookup = Dict{Int,Int}()
        for (u, v) in pi
            partner_lookup[u] = v
            partner_lookup[v] = u
        end

        choice_combinations = collect(Iterators.product(fill([1, 2], n_H ÷ 2)...))

        for choices in choice_combinations
            # Each pair has a choice of contraction term
            # choices[pair_idx] == 1: delta_il delta_jk (Port 1 <-> Port 2)
            # choices[pair_idx] == 2: delta_ik delta_jl (Port 1 <-> Port 1, Port 2 <-> Port 2)
            
            # Map m_index to its choice based on pair index
            pair_choices = Dict{Int, Int}()
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

                            # 1. Wick Jump
                            m = h_map[s]
                            partner_m = partner_lookup[m]
                            choice = pair_choices[m]
                            
                            s = H_indices[partner_m]
                            if choice == 1
                                # Standard GUE-like jump: 1->2, 2->1
                                p = (p == 1 ? 2 : 1)
                            else
                                # Transposed jump: 1->1, 2->2 (delta_ik delta_jl)
                                p = p 
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
                            push!(current_partition_traces, dim)
                        else
                            push!(current_partition_traces, tr_val(curr_factors))
                        end
                    end
                end
            end

            total_val += isempty(current_partition_traces) ? 1 : prod(current_partition_traces)
        end
    end

    return constant_part * total_val
end


function fallback_integrate(t::LazyTrace, measure::GSEMeasure)
    # GSE-GOE duality: <Tr(H^k)>_GSE(d) = (-1)^(k/2+1) <Tr(H^k)>_GOE(-d)
    all_factors = vcat(t.cycles...)
    H_type = :GSE
    n_H = count(f -> f isa SymbolicMatrix && f.special_type == H_type, all_factors)

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
                        # If orientations match, row is at same port index. If different, they swap.
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
