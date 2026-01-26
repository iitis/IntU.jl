# Gaussian Random Matrix measures (GUE, GOE, GSE)

struct GUEMeasure{M, D}
    H::M
    dim::D
end

struct GOEMeasure{M, D}
    H::M
    dim::D
end

"""
    dGUE(H, dim)

Define the measure for the Gaussian Unitary Ensemble (GUE).
`H` is the symbolic matrix representing the Hermitian Gaussian random matrix.
`dim` is the dimension (symbolic or integer).

Expectation values are defined by Wick's theorem with the contraction:
`< H_{ij} H_{kl} > = delta_{il} * delta_{jk}`

This normalization corresponds to `< Tr(H^2) > = dim^2`.
Note: H is Hermitian, so `conj(H_{ij})` is treated as `H_{ji}`.
"""
dGUE(H, dim) = GUEMeasure(H, dim)

"""
    dGOE(H, dim)

Define the measure for the Gaussian Orthogonal Ensemble (GOE).
`H` is the symbolic matrix representing the real symmetric Gaussian random matrix.
"""
dGOE(H, dim) = GOEMeasure(H, dim)

"""
    integrate(expr, measure::GUEMeasure)
"""
function integrate(expr::AbstractArray, measure::GUEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::GOEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr, measure::GUEMeasure)
    H_sym = measure.H
    dim = measure.dim
    
    subs_dict = Dict{Any, Any}()
    H_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    
    if H_sym isa AbstractArray
        for i in 1:size(H_sym, 1)
            for j in 1:size(H_sym, 2)
                h_ij_num = _safe_Num(H_sym[i,j])
                h_ij_un = Symbolics.unwrap(h_ij_num)
                h_atomic = Symbolics.variable(:H_atomic, i, j)
                
                H_atomic_lookup[Symbolics.unwrap(h_atomic)] = (i, j)
                
                subs_dict[h_ij_un] = h_atomic
                
                # Handle conjugates
                # conj(H_{ij}) = H_{ji}
                # We map conj(H) to H_bar_atomic which we will identify later as H_{ji}
                # Actually, simpler: map conj(H_{ij}) to H_atomic(j, i)
                # But we can't easily create a variable "pointing" to another variable's indices
                # implicitly without separate lookup.
                
                # So we map conj(h_{ij}) to a new atomic variable h_bar_{ij}
                # and in the lookup we store (j, i) for it!
                
                hb_atomic = Symbolics.variable(:H_bar_atomic, i, j)
                # Note: For H_{ij}, the conjugate is H_{ji}.
                # So if we see conj(H_{ij}), we treat it as H at indices (j, i).
                H_atomic_lookup[Symbolics.unwrap(hb_atomic)] = (j, i)
                
                subs_dict[Symbolics.unwrap(conj(h_ij_un))] = hb_atomic
                subs_dict[Symbolics.unwrap(Base.conj(h_ij_un))] = hb_atomic
            end
        end
    end

    # For GUE, we don't separate U and U_bar. All are H.
    # We pass empty U_bar_lookup because we mapped everything to H_atomic_lookup (or compatible)
    return _integrate_core(expr, dim, subs_dict, H_atomic_lookup, Dict(), :GUE)
end

function integrate(t::LazyTrace, measure::GUEMeasure)
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
    wires = Dict{Int, Any}()
    for k in 1:n_H
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
        perm_map = Dict{Int, Int}()
        for (u, v) in pi
            perm_map[u] = v
            perm_map[v] = u
        end
        
        # Count cycles in the resulting connection graph
        visited = falses(n_H)
        current_partition_traces = []
        
        for start_m in 1:n_H
            if !visited[start_m]
                curr_trace_factors = SymbolicMatrix[]
                curr_m = start_m
                while !visited[curr_m]
                    visited[curr_m] = true
                    # Traverse from H_curr_m to its paired partner
                    paired_m = perm_map[curr_m]
                    # The wire starts AFTER paired_m and leads to some other H
                    # but wait...
                    # <H_ij H_kl> means output of u connects to input of v.
                    # Output of u is the wire starting at H_indices[u].
                    # Input of v is the wire ending at H_indices[v].
                    
                    # So start at current slot, take the wire, land at next H.
                    # Then take paired partner of that H.
                    
                    # Wire from curr_m:
                    dest_factor_idx, mat_segment = wires[H_indices[curr_m]]
                    if mat_segment !== nothing append!(curr_trace_factors, mat_segment) end
                    
                    # We landed at H with factor index `dest_factor_idx`.
                    # Let's find its m-index.
                    next_m = 1
                    while H_indices[next_m] != dest_factor_idx; next_m += 1; end
                    
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

function integrate(expr, measure::GOEMeasure)
    H_sym = measure.H
    dim = measure.dim
    
    subs_dict = Dict{Any, Any}()
    # Reusing H_atomic_lookup for GOE variables
    H_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    
    if H_sym isa AbstractArray
        for i in 1:size(H_sym, 1)
            for j in 1:size(H_sym, 2)
                h_ij_num = _safe_Num(H_sym[i,j])
                h_ij_un = Symbolics.unwrap(h_ij_num)
                h_atomic = Symbolics.variable(:H_atomic, i, j)
                
                H_atomic_lookup[Symbolics.unwrap(h_atomic)] = (i, j)
                
                subs_dict[h_ij_un] = h_atomic
                
                # Handle conjugates for GOE (Real Symmetric)
                # conj(H_{ij}) = H_{ij} = H_{ji}
                # But to maintain structure, we map conj to same atomic variable.
                # However, logic in core might handle conjugate "slots".
                # Actually, simpler: map conj(h) -> h_atomic.
                # AND map H[j,i] to same variables if i != j?
                # Usually user defines H as symmetric matrix?
                # If H is just generic symbolic array, we impose symmetry via the measure logic.
                
                subs_dict[Symbolics.unwrap(conj(h_ij_un))] = h_atomic
                subs_dict[Symbolics.unwrap(Base.conj(h_ij_un))] = h_atomic
            end
        end
    end

    # For GOE, measure type :GOE
    return _integrate_core(expr, dim, subs_dict, H_atomic_lookup, Dict(), :GOE)
end

function integrate(t::LazyTrace, measure::GOEMeasure)
    factors = t.factors
    H_name = measure.H isa SymbolicMatrix ? measure.H.name : :H
    
    H_indices = Int[]
    for (i, f) in enumerate(factors)
        if f.name == H_name
            push!(H_indices, i)
        end
    end
    
    n_H = length(H_indices)
    if isodd(n_H); return 0; end
    if n_H == 0; return tr_val(factors); end
    
    dim = measure.dim
    n_factors = length(factors)
    
    # Pre-build wires
    wires = Dict{Int, Any}() # index -> (dest_idx, segment, rev_segment)
    # rev_segment would be adjoints but H is symmetric so same.
    for k in 1:n_H
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
        
        wires[idx] = (next_h_idx, isempty(fwd_consts) ? nothing : fwd_consts, 
                      prev_h_idx, isempty(bwd_consts) ? nothing : bwd_consts)
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
            
            for start_m in 1:n_H
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
                                if fwd !== nothing append!(curr_trace_factors, fwd) end
                                
                                # Land at Port 1 (Input) of some H
                                landed_m = 1
                                while H_indices[landed_m] != dest_idx; landed_m += 1; end
                                visited_ports[landed_m, 1] = true
                                
                                # Use Wick contraction jump
                                pair_idx = 0; partner_m = 0
                                for (p_idx, (u, v)) in enumerate(pi)
                                    if u == landed_m; pair_idx = p_idx; partner_m = v; break; end
                                    if v == landed_m; pair_idx = p_idx; partner_m = u; break; end
                                end
                                
                                if choices[pair_idx] == 2 # delta_il delta_jk (P1 -> P2)
                                    curr_m = partner_m; curr_port = 2
                                else # delta_ik delta_jl (P1 -> P1)
                                    curr_m = partner_m; curr_port = 1
                                end
                            else
                                # Exit Port 1 (Input), take backward wire
                                dest_idx, fwd, prev_idx, bwd = wires[H_indices[curr_m]]
                                if bwd !== nothing append!(curr_trace_factors, bwd) end
                                
                                # Land at Port 2 (Output) of some H
                                landed_m = 1
                                while H_indices[landed_m] != prev_idx; landed_m += 1; end
                                visited_ports[landed_m, 2] = true
                                
                                pair_idx = 0; partner_m = 0
                                for (p_idx, (u, v)) in enumerate(pi)
                                    if u == landed_m; pair_idx = p_idx; partner_m = v; break; end
                                    if v == landed_m; pair_idx = p_idx; partner_m = u; break; end
                                end
                                
                                if choices[pair_idx] == 2 # delta_il delta_jk (P2 -> P1)
                                    curr_m = partner_m; curr_port = 1
                                else # delta_ik delta_jl (P2 -> P2)
                                    curr_m = partner_m; curr_port = 2
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
            
            # Since each cycle is counted twice (once per direction for real symmetric edges? 
            # No, GOE trace graph is directed by the trace cycle but Wick edges are undirected.
            # Actually, standard Wick logic for GOE traces:
            # Each choices set corresponds to a single term.
            # But the way I count cycles (starting from any port) might double count?
            # Yes, if I start at Port 2 and go fwd, versus start at Port 1 and go bwd.
            # In GOE, the "cycle" is undirected? No, tr(factors) is symmetric if matrices are symmetric.
            # But SymbolicMatrix' is explicitly there.
            
            # Re-eval: In GUE, we had directed cycles. In GOE, we have undirected cycles?
            # But we are in a trace, which has a direction.
            # If we go backwards, we get tr(factors'). 
            # If factors are real symmetric, tr(factors') = tr(factors).
            
            # To avoid double counting, we observe that each port is connected to exactly one other port.
            # So the graph is a collection of disjoint cycles.
            # My loop finds each cycle.
            # If I sum over ports, I visit each cycle of length L exactly L times?
            # No, if I visit each node in the cycle once.
            
            # The result for one choice set is prod(current_partition_traces)
            term_val = isempty(current_partition_traces) ? 1 : prod(current_partition_traces)
            total_val += term_val
        end
    end
    return total_val
end

"""
    asymptotic(expr, measure::GUEMeasure, order=1)
"""
function asymptotic(expr, measure::GUEMeasure, order=1)
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
function asymptotic(expr, measure::GOEMeasure, order=1)
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
