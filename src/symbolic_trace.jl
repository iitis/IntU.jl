# src/symbolic_trace.jl

"""
    SymbolicMatrix(name::Symbol)
    SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol)

A wrapper associated with a symbolic name to represent a matrix in a coordinate-free way.
`special_type` can be `:U`, `:U_dag`, or `:Constant`.
"""
struct SymbolicMatrix
    name::Symbol
    is_adj::Bool
    special_type::Symbol  # :U, :U_dag, :Constant
end

SymbolicMatrix(name::Symbol) = SymbolicMatrix(name, false, :Constant)

import Base: *, adjoint, show

function adjoint(A::SymbolicMatrix)
    new_type = A.special_type
    if A.special_type == :U
        new_type = :U_dag
    elseif A.special_type == :U_dag
        new_type = :U
    end
    # For constants, new_type remains :Constant (or we could have :Constant_dag, but :Constant handles both usually)
    return SymbolicMatrix(A.name, !A.is_adj, new_type)
end

function show(io::IO, A::SymbolicMatrix)
    print(io, A.name)
    if A.is_adj
        print(io, "'")
    end
end

"""
    LazyTrace

Represents the trace of a sequence of SymbolicMatrices.
"""
struct LazyTrace
    factors::Vector{SymbolicMatrix}
end

function *(A::SymbolicMatrix, B::SymbolicMatrix)
    return [A, B]
end

function *(A::Vector{SymbolicMatrix}, B::SymbolicMatrix)
    return push!(copy(A), B)
end

function *(A::SymbolicMatrix, B::Vector{SymbolicMatrix})
    return vcat([A], B)
end

"""
    tr_lazy(product)

Create a LazyTrace from a product of SymbolicMatrices (or a single one).
"""
function tr_lazy(product::AbstractVector)
    return LazyTrace(collect(SymbolicMatrix, product))
end

function tr_lazy(product::SymbolicMatrix)
    return LazyTrace([product])
end

function show(io::IO, t::LazyTrace)
    print(io, "tr(")
    for (i, f) in enumerate(t.factors)
        print(io, f)
        if i < length(t.factors)
            print(io, " * ")
        end
    end
    print(io, ")")
end




"""
    tr_val(p::SymbolicMatrixProduct)

Symbolic function representing the trace of a matrix product.
"""
function tr_val(factors::Vector{SymbolicMatrix})
    # Create a symbolic variable representing this trace
    # This avoids issues with Term multiplication and simplification
    if isempty(factors)
        return 1 # Should handle dim separately, but trace of empty is not passed here usually
    end
    
    # Construct a nice string representation
    # e.g. "tr(A * B)"
    s_parts = String[]
    for (i, f) in enumerate(factors)
        push!(s_parts, string(f))
    end
    name = "tr_val(" * join(s_parts, "*") * ")"
    return Symbolics.variable(Symbol(name); T=Real)
end
# Symbolics metadata might go here if needed.
# Actually, we rely on Term wrapping it.

"""
    integrate(t::LazyTrace, measure::HaarMeasure)

Integrate a single trace of matrices over the Haar measure.
Uses the graphical Weingarten calculus.
"""
function integrate(t::LazyTrace, measure::HaarMeasure)
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
            
            if isequal(wg_val, 0)
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
