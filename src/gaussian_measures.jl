# Gaussian Random Matrix measures (GUE, GOE, GSE)

struct GUEMeasure{M,D}
    H::M
    dim::D
end

struct GOEMeasure{M,D}
    H::M
    dim::D
end

struct GSEMeasure{M,D}
    H::M
    dim::D
end

struct GinUEMeasure{M,D}
    G::M
    dim::D
end

struct GinOEMeasure{M,D}
    G::M
    dim::D
end

struct GinSEMeasure{M,D}
    G::M
    dim::D
end

"""
    dGUE(H, dim)

Defines the measure for the **Gaussian Unitary Ensemble (GUE)**.
H is a complex Hermitian matrix (H = H^\\dagger).

Expectation values are given by **Wick's Theorem** (Isserlis' Theorem) with the contraction:
```math
\\langle H_{ij} \\bar{H}_{kl} \\rangle = \\delta_{il} \\delta_{jk}
```
This implies \\langle \\text{Tr}(H^2) \\rangle = d^2.

Reference:
- Mehta, M. L. (2004). *Random Matrices*.
"""
dGUE(H, dim) = GUEMeasure(H, dim)

"""
    dGOE(H, dim)

Defines the measure for the **Gaussian Orthogonal Ensemble (GOE)**.
H is a real symmetric matrix (H = H^T).

The Wick contraction is:
```math
\\langle H_{ij} H_{kl} \\rangle = \\delta_{ik} \\delta_{jl} + \\delta_{il} \\delta_{jk}
```
This implies \\langle \\text{Tr}(H^2) \\rangle = d^2 + d.

Reference:
- Mehta, M. L. (2004). *Random Matrices*.
"""
dGOE(H, dim) = GOEMeasure(H, dim)

"""
    dGSE(H, dim)

Defines the measure for the **Gaussian Symplectic Ensemble (GSE)**.
H is a Hermitian quaternionic self-dual matrix. Dimension d must be even.

The integration uses the property:
```math
\\langle \\text{Tr}(H^k) \\rangle_{GSE}(d) = (-1)^{\\frac{k}{2} + 1} \\langle \\text{Tr}(H^k) \\rangle_{GOE}(-d)
```
This implies \\langle \\text{Tr}(H^2) \\rangle = d^2 - d.

Reference:
- Mehta, M. L. (2004). *Random Matrices*.
"""
function dGSE(H, dim)
    if dim isa Integer && isodd(dim)
        throw(ArgumentError("GSE dimension must be even, got $dim"))
    end
    return GSEMeasure(H, dim)
end

"""
    dGinUE(G, dim)

Defines the measure for the **Complex Ginibre Ensemble (GinUE)**.
G is a general complex matrix with i.i.d. complex Gaussian entries.

Contractions are given by:
```math
\\langle G_{ij} \\bar{G}_{kl} \\rangle = \\delta_{ik} \\delta_{jl}
```
All other contractions (including \$\\langle G_{ij} G_{kl} \\rangle\$) are zero.
"""
dGinUE(G, dim) = GinUEMeasure(G, dim)

"""
    dGinOE(G, dim)

Defines the measure for the **Real Ginibre Ensemble (GinOE)**.
G is a general real matrix with i.i.d. real Gaussian entries.

The Wick contraction is:
```math
\\langle G_{ij} G_{kl} \\rangle = \\delta_{ik} \\delta_{jl}
```
"""
dGinOE(G, dim) = GinOEMeasure(G, dim)

"""
    dGinSE(G, dim)

Defines the measure for the **Symplectic Ginibre Ensemble (GinSE)**.
G is a quaternionic matrix. Dimension d must be even.

Uses duality relation with GinOE.
"""
function dGinSE(G, dim)
    if dim isa Integer && isodd(dim)
        throw(ArgumentError("GinSE dimension must be even, got $dim"))
    end
    return GinSEMeasure(G, dim)
end

function _setup_gaussian_subs(H_sym, ensemble_type)
    subs_dict = Dict{Any,Any}()
    H_atomic_lookup = Dict{Any,Tuple}()

    if H_sym isa AbstractArray
        for i = 1:size(H_sym, 1)
            for j = 1:size(H_sym, 2)
                h_ij_num = _safe_Num(H_sym[i, j])
                h_ij_un = Symbolics.unwrap(h_ij_num)
                h_atomic = Symbolics.variable(:H_atomic, i, j)

                H_atomic_lookup[Symbolics.unwrap(h_atomic)] = (i, j)
                subs_dict[h_ij_un] = h_atomic

                if ensemble_type in (:GUE, :GSE)
                    # Hermitian: conj(H_{ij}) = H_{ji}
                    hb_atomic = Symbolics.variable(:H_bar_atomic, i, j)
                    H_atomic_lookup[Symbolics.unwrap(hb_atomic)] = (j, i)
                    subs_dict[Symbolics.unwrap(conj(h_ij_un))] = hb_atomic
                    subs_dict[Symbolics.unwrap(Base.conj(h_ij_un))] = hb_atomic
                elseif ensemble_type == :GOE
                    # Real symmetric: conj(H_{ij}) = H_{ij}
                    subs_dict[Symbolics.unwrap(conj(h_ij_un))] = h_atomic
                    subs_dict[Symbolics.unwrap(Base.conj(h_ij_un))] = h_atomic
                elseif ensemble_type == :GinUE
                    # Non-Hermitian complex: conj(G_{ij}) = G_bar_{ij}
                    gb_atomic = Symbolics.variable(:G_bar_atomic, i, j)
                    H_atomic_lookup[Symbolics.unwrap(gb_atomic)] = (i, j, :conj)
                    subs_dict[Symbolics.unwrap(conj(h_ij_un))] = gb_atomic
                    subs_dict[Symbolics.unwrap(Base.conj(h_ij_un))] = gb_atomic
                elseif ensemble_type in (:GinOE, :GinSE)
                    # GinOE: RealEntries, conj(G_{ij}) = G_{ij}
                    # GinSE: handled via GOE/GSE-like logic or simply Real mapping if we use Wick.
                    subs_dict[Symbolics.unwrap(conj(h_ij_un))] = h_atomic
                    subs_dict[Symbolics.unwrap(Base.conj(h_ij_un))] = h_atomic
                end
            end
        end
    end
    return subs_dict, H_atomic_lookup
end

"""
    integrate(expr, measure::GUEMeasure)
"""
function integrate(expr::AbstractArray, measure::GUEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::GOEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::GSEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::GinUEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::GinOEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::GinSEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function IntU.measure_info(measure::GUEMeasure)
    subs_dict, H_atomic_lookup = _setup_gaussian_subs(measure.H, :GUE)
    matcher = LookupMatcher(H_atomic_lookup, Dict{Any,Tuple}())
    return (subs_dict, matcher, measure.dim, :GUE)
end

function IntU.measure_info(measure::GinUEMeasure)
    subs_dict, H_atomic_lookup = _setup_gaussian_subs(measure.G, :GinUE)
    matcher = LookupMatcher(H_atomic_lookup, Dict{Any,Tuple}())
    return (subs_dict, matcher, measure.dim, :GinUE)
end

function IntU.measure_info(measure::GinOEMeasure)
    subs_dict, H_atomic_lookup = _setup_gaussian_subs(measure.G, :GinOE)
    matcher = LookupMatcher(H_atomic_lookup, Dict{Any,Tuple}())
    return (subs_dict, matcher, measure.dim, :GinOE)
end

function IntU.measure_info(measure::GinSEMeasure)
    subs_dict, H_atomic_lookup = _setup_gaussian_subs(measure.G, :GinSE)
    matcher = LookupMatcher(H_atomic_lookup, Dict{Any,Tuple}())
    return (subs_dict, matcher, measure.dim, :GinSE)
end

function fallback_integrate(t::LazyTrace, measure::GUEMeasure)
    factors = t.factors
    H_name = measure.H isa SymbolicMatrix ? measure.H.name : :H

    # Identify indices of H factors
    H_indices = Int[]
    for (i, f) in enumerate(factors)
        if f.name == H_name
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

    # Build wires between H slots
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
        # Each pairing (u, v) in pi connects H_indices[u] and H_indices[v]
        # In GUE, <H_ij H_kl> = delta_il delta_jk.
        # This means we connect output(u) to input(v) AND input(u) to output(v).
        # This is equivalent to connecting tracing through the cycle.

        # We can map this to a permutation in S_{n_H}.
        # Pairing (u, v) corresponds to a swap (u v) in terms of connections.
        perm_map = Dict{Int,Int}()
        for (u, v) in pi
            perm_map[u] = v
            perm_map[v] = u
        end

        # Count cycles in the resulting connection graph
        visited = falses(n_H)
        current_partition_traces = []

        for start_m = 1:n_H
            if !visited[start_m]
                curr_trace_factors = SymbolicMatrix[]
                curr_m = start_m
                while !visited[curr_m]
                    visited[curr_m] = true
                    # Traverse from H_curr_m to its paired partner
                    paired_m = perm_map[curr_m]
                    # The wire starts AFTER paired_m and leads to some other H
                    # <H_ij H_kl> means output of u connects to input of v.
                    # Output of u is the wire starting at H_indices[u].
                    # Input of v is the wire ending at H_indices[v].

                    # So start at current slot, take the wire, land at next H.
                    # Then take paired partner of that H.

                    # Wire from curr_m:
                    dest_factor_idx, mat_segment = wires[H_indices[curr_m]]
                    if mat_segment !== nothing
                        append!(curr_trace_factors, mat_segment)
                    end

                    # We landed at H with factor index `dest_factor_idx`.
                    # Let's find its m-index.
                    next_m = 1
                    while H_indices[next_m] != dest_factor_idx
                        ;
                        next_m += 1;
                    end

                    # Now we are at next_m, but Wick contraction says we jump to its partner!
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

function IntU.measure_info(measure::GOEMeasure)
    subs_dict, H_atomic_lookup = _setup_gaussian_subs(measure.H, :GOE)
    matcher = LookupMatcher(H_atomic_lookup, Dict{Any,Tuple}())
    return (subs_dict, matcher, measure.dim, :GOE)
end

function fallback_integrate(t::LazyTrace, measure::GOEMeasure)
    factors = t.factors
    H_name = measure.H isa SymbolicMatrix ? measure.H.name : :H

    H_indices = Int[]
    for (i, f) in enumerate(factors)
        if f.name == H_name
            push!(H_indices, i)
        end
    end

    n_H = length(H_indices)
    if isodd(n_H)
        ;
        return 0;
    end
    if n_H == 0
        ;
        return tr_val(factors);
    end

    dim = measure.dim
    n_factors = length(factors)

    # Pre-build wires
    wires = Dict{Int,Any}() # index -> (dest_idx, segment, rev_segment)
    # rev_segment would be adjoints but H is symmetric so same.
    for k = 1:n_H
        idx = H_indices[k]
        next_h_idx = H_indices[mod1(k+1, n_H)]
        prev_h_idx = H_indices[mod1(k-1, n_H)]

        # Forward segment
        fwd_consts = SymbolicMatrix[]
        curr = mod1(idx + 1, n_factors)
        while curr != next_h_idx
            push!(fwd_consts, factors[curr])
            curr = mod1(curr + 1, n_factors)
        end

        # Backward segment (for GOE swaps)
        # Traverses from idx backwards to prev_h_idx
        # This is used if we enter H_idx from its "output" port and want to go to its "input" port.
        # But for symmetric traces, it's just the reverse order of factors.
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
        # Each pairing contributes 2^(n_H/2) terms? No.
        # <H_ij H_kl> = delta_ik delta_jl + delta_il delta_jk.
        # We must sum over all 2^(n_H/2) choices of contraction types.

        choice_combinations = collect(Iterators.product(fill([1, 2], n_H ÷ 2)...))

        for choices in choice_combinations
            # choices is a vector of 1 or 2 for each pair in pi
            visited = falses(n_H)
            current_partition_traces = []

            # We track visited (m, port)
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
                                    ;
                                    landed_m += 1;
                                end
                                visited_ports[landed_m, 1] = true

                                # Use Wick contraction jump
                                pair_idx = 0;
                                partner_m = 0
                                for (p_idx, (u, v)) in enumerate(pi)
                                    if u == landed_m
                                        ;
                                        pair_idx = p_idx;
                                        partner_m = v;
                                        break;
                                    end
                                    if v == landed_m
                                        ;
                                        pair_idx = p_idx;
                                        partner_m = u;
                                        break;
                                    end
                                end

                                if choices[pair_idx] == 2 # delta_il delta_jk (P1 -> P2)
                                    curr_m = partner_m;
                                    curr_port = 2
                                else # delta_ik delta_jl (P1 -> P1)
                                    curr_m = partner_m;
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
                                    ;
                                    landed_m += 1;
                                end
                                visited_ports[landed_m, 2] = true

                                pair_idx = 0;
                                partner_m = 0
                                for (p_idx, (u, v)) in enumerate(pi)
                                    if u == landed_m
                                        ;
                                        pair_idx = p_idx;
                                        partner_m = v;
                                        break;
                                    end
                                    if v == landed_m
                                        ;
                                        pair_idx = p_idx;
                                        partner_m = u;
                                        break;
                                    end
                                end

                                if choices[pair_idx] == 2 # delta_il delta_jk (P2 -> P1)
                                    curr_m = partner_m;
                                    curr_port = 1
                                else # delta_ik delta_jl (P2 -> P2)
                                    curr_m = partner_m;
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

            # Each (partition, choice-set) pair contributes one term: prod of cycle traces.

            # The result for one choice set is prod(current_partition_traces)
            term_val =
                isempty(current_partition_traces) ? 1 : prod(current_partition_traces)
            total_val += term_val
        end
    end
    return total_val
end

function IntU.measure_info(measure::GSEMeasure)
    subs_dict, H_atomic_lookup = _setup_gaussian_subs(measure.H, :GSE)
    matcher = LookupMatcher(H_atomic_lookup, Dict{Any,Tuple}())
    return (subs_dict, matcher, measure.dim, :GSE)
end

function fallback_integrate(t::LazyTrace, measure::GSEMeasure)
    # The moment of GSE is related to GOE by duality:
    # <Tr(H^k)>_GSE(d) = (-1)^(k/2 + 1) * <Tr(H^k)>_GOE(-d)
    # For multiple traces, it's more complex, but for a single LazyTrace:
    factors = t.factors
    H_name = measure.H isa SymbolicMatrix ? measure.H.name : :H
    n_H = count(f -> f.name == H_name, factors)

    if isodd(n_H)
        ;
        return 0;
    end

    # 1. Integrate as GOE
    goe_res = integrate(t, dGOE(measure.H, measure.dim))

    # 2. Substitution d -> -d
    dim = measure.dim
    # We substitute both instances if dim is Num or just handle it
    res_subbed = Symbolics.substitute(goe_res, Dict(dim => -dim))

    # 3. Apply overall sign
    final_sign = ((-1)^(n_H ÷ 2 + 1))
    return final_sign * res_subbed
end

function fallback_integrate(t::LazyTrace, measure::GinUEMeasure)
    dim = measure.dim
    if isempty(t.cycles)
        return t.prefactor
    end

    G_name = measure.G isa SymbolicMatrix ? measure.G.name : :G

    # 1. Collect all factors and find Gaussian slots
    all_factors = SymbolicMatrix[]
    cycle_ranges = UnitRange{Int}[]
    total_factors = 0
    for cycle in t.cycles
        start_idx = total_factors + 1
        append!(all_factors, cycle)
        total_factors += length(cycle)
        push!(cycle_ranges, start_idx:total_factors)
    end

    # GinUE: <G_ij G_bar_kl> = delta_ik delta_jl
    G_indices = Int[]
    G_bar_indices = Int[]
    for (i, f) in enumerate(all_factors)
        if f.name == G_name
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

    wires = _build_wires(G_indices, G_bar_indices, cycle_ranges, all_factors)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    perms = permutations(1:n_G_bar)
    total_val = 0
    pos_map = Dict(idx => i for (i, idx) in enumerate(all_slots))
    
    # Map from G index to its position in all_slots
    # Map from G_bar index to its position in all_slots
    
    for p in perms
        visited = falses(length(all_slots))
        current_partition_traces = []

        for start_pos = 1:length(all_slots)
            if !visited[start_pos]
                curr_trace_factors = SymbolicMatrix[]
                curr_pos = start_pos
                while !visited[curr_pos]
                    visited[curr_pos] = true
                    curr_factor_idx = all_slots[curr_pos]
                    
                    dest_factor_idx, mat_segment = wires[curr_factor_idx]
                    if mat_segment !== nothing
                        append!(curr_trace_factors, mat_segment)
                    end
                    
                    # Landed at dest_factor_idx. Find its partner.
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

    G_name = measure.G isa SymbolicMatrix ? measure.G.name : :G
    
    all_factors = SymbolicMatrix[]
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
        if f.name == G_name
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
    
    wires = _build_wires(G_indices, Int[], cycle_ranges, all_factors)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    
    partitions = get_pair_partitions(n_G)
    total_val = 0
    pos_map = Dict(idx => i for (i, idx) in enumerate(G_indices))
    
    for pi in partitions
        visited = falses(n_G)
        current_partition_traces = []
        
        for start_m = 1:n_G
            if !visited[start_m]
                curr_trace_factors = SymbolicMatrix[]
                curr_m = start_m
                while !visited[curr_m]
                    visited[curr_m] = true
                    dest_factor_idx, mat_segment = wires[G_indices[curr_m]]
                    if mat_segment !== nothing
                        append!(curr_trace_factors, mat_segment)
                    end
                    landed_m = pos_map[dest_factor_idx]
                    
                    # Find Wick partner
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
                    curr_m = partner_m
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

function fallback_integrate(t::LazyTrace, measure::GinSEMeasure)
    G_name = measure.G isa SymbolicMatrix ? measure.G.name : :G
    all_factors = vcat(t.cycles...)
    n_G = count(f -> f.name == G_name, all_factors)

    if isodd(n_G)
        return 0
    end

    ginoe_res = integrate(t, dGinOE(measure.G, measure.dim))
    dim = measure.dim
    res_subbed = Symbolics.substitute(ginoe_res, Dict(dim => -dim))
    final_sign = ((-1)^(n_G ÷ 2 + 1))
    return final_sign * res_subbed
end

"""
    asymptotic(expr, measure::GUEMeasure, order=1)
"""
function asymptotic(expr, measure::GUEMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dGUE(measure.H, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end

"""
    asymptotic(expr, measure::GOEMeasure, order=1)
"""
function asymptotic(expr, measure::GOEMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dGOE(measure.H, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end

"""
    asymptotic(expr, measure::GSEMeasure, order=1)
"""
function asymptotic(expr, measure::GSEMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dGSE(measure.H, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end

"""
    asymptotic(expr, measure::GinUEMeasure, order=1)
"""
function asymptotic(expr, measure::GinUEMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dGinUE(measure.G, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end

"""
    asymptotic(expr, measure::GinOEMeasure, order=1)
"""
function asymptotic(expr, measure::GinOEMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dGinOE(measure.G, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end

"""
    asymptotic(expr, measure::GinSEMeasure, order=1)
"""
function asymptotic(expr, measure::GinSEMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dGinSE(measure.G, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end
