# Haar measure integration

# Dummy type to represent the measure
struct HaarMeasure{M, D}
    U::M
    dim::D
end
"""
    dU(U, dim)

Define the Haar measure for the Unitary group U(d).
`U` is the symbolic matrix representing the unitary, and `dim` is the dimension (symbolic or integer).
"""
dU(U, dim) = HaarMeasure(U, dim)

"""
    integrate(expr, measure::HaarMeasure)
"""
function integrate(expr::AbstractArray, measure::HaarMeasure)
    return map(e -> integrate(e, measure), expr)
end

function fallback_integrate(expr, measure::HaarMeasure)
    U_sym = measure.U
    dim = measure.dim
    
    # Substitute Re(U) and Im(U)
    subs_dict = Dict{Any, Any}()
    U_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    U_bar_lookup = Dict{Any, Tuple{Int, Int}}()
    

    if U_sym isa AbstractArray
        for i in 1:size(U_sym, 1)
            for j in 1:size(U_sym, 2)
                u_ij_num = _safe_Num(U_sym[i,j])
                u_ij_un = Symbolics.unwrap(u_ij_num)
                u_atomic = Symbolics.variable(:U_atomic, i, j)
                u_bar_atomic = Symbolics.variable(:U_bar_atomic, i, j)
                
                U_atomic_lookup[Symbolics.unwrap(u_atomic)] = (i, j)
                U_bar_lookup[Symbolics.unwrap(u_bar_atomic)] = (i, j)
                
                subs_dict[u_ij_un] = u_atomic
                
                c_ij_un = Symbolics.unwrap(conj(u_ij_num))
                subs_dict[c_ij_un] = u_bar_atomic
                
                bc_ij_un = Symbolics.unwrap(Base.conj(u_ij_num))
                subs_dict[bc_ij_un] = u_bar_atomic
            end
        end
    end

    return _robust_real_num(_integrate_core(expr, dim, subs_dict, U_atomic_lookup, U_bar_lookup))
end

"""
    asymptotic(expr, measure::HaarMeasure, order=1)

Returns the series expansion of the integral in powers of `1/d`.
"""
function asymptotic(expr, measure::HaarMeasure, order=1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end
    
    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dU(measure.U, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end

"""
    integrate(t::LazyTrace, measure::HaarMeasure)

Integrate a single trace of matrices over the Haar measure.
Uses the graphical Weingarten calculus.
"""
function fallback_integrate(t::LazyTrace, measure::HaarMeasure)
    # 1. Identify U and U_dag instances
    U_indices = Int[]
    U_bar_indices = Int[]
    
    factors = t.factors
    n_factors = length(factors)
    
    for (i, f) in enumerate(factors)
        if f.special_type == :U
            push!(U_indices, i)
        elseif f.special_type == :U_dag
            push!(U_bar_indices, i)
        end
    end
    
    n_U = length(U_indices)
    n_U_bar = length(U_bar_indices)
    
    if n_U != n_U_bar
        return 0
    end
    
    dim = measure.dim
    
    if n_U == 0
        if isempty(factors)
            return dim
        end
        return tr_val(factors)
    end
    
    
    # 2. Build Wires
    wires = Dict{Int, Any}()
    all_slots = sort([U_indices; U_bar_indices])
    n_slots = length(all_slots)
    
    for k in 1:n_slots
        start_idx = all_slots[k]
        end_idx = all_slots[mod1(k+1, n_slots)]
        
        consts = SymbolicMatrix[]
        curr = mod1(start_idx + 1, n_factors)
        while curr != end_idx
            push!(consts, factors[curr])
            curr = mod1(curr + 1, n_factors)
        end
        wires[start_idx] = (end_idx, isempty(consts) ? nothing : consts)
    end
    
    u_map = Dict{Int, Int}()
    for (m, idx) in enumerate(U_indices) u_map[idx] = m end
    ub_map = Dict{Int, Int}()
    for (m, idx) in enumerate(U_bar_indices) ub_map[idx] = m end
    
    permutations = collect(Combinatorics.permutations(1:n_U))
    total_val = 0
    
    for sigma in permutations
        for tau in permutations
             
            inv_tau = invperm(tau)
            P = [sigma[inv_tau[i]] for i in 1:n_U]
            cycle_type = get_cycle_type(P)
            wg_val = weingarten(cycle_type, dim)
            
            if _symbolic_isequal(wg_val, 0)
                continue
            end
            
            visited_U = falses(n_U)
            visited_Ub = falses(n_U)
            current_term_traces = []
            
            for start_m in 1:n_U
                if !visited_U[start_m]
                    curr_trace_factors = SymbolicMatrix[]
                    curr_type = 1 
                    curr_idx = start_m
                    
                    while true
                        if curr_type == 1
                            if visited_U[curr_idx] break end
                            visited_U[curr_idx] = true
                            next_ub_m = sigma[curr_idx]
                            start_factor_idx = U_bar_indices[next_ub_m]
                            dest_factor_idx, mat_segment = wires[start_factor_idx]
                            if mat_segment !== nothing append!(curr_trace_factors, mat_segment) end
                            if haskey(u_map, dest_factor_idx)
                                curr_type = 1; curr_idx = u_map[dest_factor_idx]
                            else
                                curr_type = 2; curr_idx = ub_map[dest_factor_idx]
                            end
                        else
                            if visited_Ub[curr_idx] break end
                            visited_Ub[curr_idx] = true
                            next_u_m = inv_tau[curr_idx]
                            start_factor_idx = U_indices[next_u_m]
                            dest_factor_idx, mat_segment = wires[start_factor_idx]
                            if mat_segment !== nothing append!(curr_trace_factors, mat_segment) end
                            if haskey(u_map, dest_factor_idx)
                                curr_type = 1; curr_idx = u_map[dest_factor_idx]
                            else
                                curr_type = 2; curr_idx = ub_map[dest_factor_idx]
                            end
                        end
                    end
                    if isempty(curr_trace_factors)
                        push!(current_term_traces, dim)
                    else
                        push!(current_term_traces, tr_val(curr_trace_factors))
                    end
                end
            end
            
            # Check Ub cycles (if any isolated ones exist - theoretically shouldn't for connected trace)
             for start_m in 1:n_U_bar
                 if !visited_Ub[start_m]
                    curr_trace_factors = SymbolicMatrix[]
                    curr_type = 2
                    curr_idx = start_m
                    while true
                        if curr_type == 1
                            if visited_U[curr_idx] break end
                            visited_U[curr_idx] = true
                            next_ub_m = sigma[curr_idx]
                            start_factor_idx = U_bar_indices[next_ub_m]
                            dest_factor_idx, mat_segment = wires[start_factor_idx]
                            if mat_segment !== nothing append!(curr_trace_factors, mat_segment) end
                            if haskey(u_map, dest_factor_idx)
                                curr_type = 1; curr_idx = u_map[dest_factor_idx]
                            else
                                curr_type = 2; curr_idx = ub_map[dest_factor_idx]
                            end
                        else
                            if visited_Ub[curr_idx] break end
                            visited_Ub[curr_idx] = true
                            next_u_m = inv_tau[curr_idx]
                            start_factor_idx = U_indices[next_u_m]
                            dest_factor_idx, mat_segment = wires[start_factor_idx]
                            if mat_segment !== nothing append!(curr_trace_factors, mat_segment) end
                            if haskey(u_map, dest_factor_idx)
                                curr_type = 1; curr_idx = u_map[dest_factor_idx]
                            else
                                curr_type = 2; curr_idx = ub_map[dest_factor_idx]
                            end
                        end
                    end
                     if isempty(curr_trace_factors)
                        push!(current_term_traces, dim)
                    else
                         push!(current_term_traces, tr_val(curr_trace_factors))
                    end
                 end
            end
            
            term_prod = isempty(current_term_traces) ? 1 : prod(current_term_traces)
            total_val += term_prod * wg_val
        end
    end
    
    return total_val
end
