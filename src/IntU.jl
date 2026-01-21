module IntU

using Symbolics
using Combinatorics
using LinearAlgebra
using SymbolicUtils

include("Weingarten.jl")
using .Weingarten

export integrate, dU, integrate_indices

# Dummy type to represent the measure
struct HaarMeasure{T, N, D}
    U::AbstractArray{T, N}
    dim::D
end
dU(U::AbstractArray{T,N}, dim) where {T,N} = HaarMeasure{T,N,typeof(dim)}(U, dim)

"""
    integrate(expr, measure::HaarMeasure)
"""
function integrate(expr::AbstractArray, measure::HaarMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr, measure::HaarMeasure)
    if expr isa Complex
        val_re = integrate(real(expr), measure)
        val_im = integrate(imag(expr), measure)
        return _robust_real(val_re + im * val_im)
    end
    
    U_sym = measure.U
    dim = measure.dim
    
    IM_dummy = Symbolics.variable(:IM_dummy)
    
    # Expand FIRST to reveal Re(u) and Im(u) terms from products like U*conj(U)
    expr = Symbolics.expand(expr)

    # Predicate to check if term is a Sum
    function is_add(t)
        Symbolics.iscall(t) && Symbolics.operation(t) == (+)
    end

    # Rewrite rules to ensure abs(z)^2 becomes (z * conj(z)) or real^2 + imag^2
    r_abs_sq = @rule abs(~x)^2 => (~x) * conj(~x)
    r_abs = @rule abs(~x) => hypot(real(~x), imag(~x))
    
    # Rewrite hypot(x,y) -> ((x + i*y)(x - i*y))^(1//2) ONLY if x, y are Sums
    r_hypot_sum = @rule hypot(~x::is_add, ~y::is_add) => ((~x + im*~y) * (~x - im*~y))^(1//2)
    # Fallback/Default for non-sums:
    r_hypot_default = @rule hypot(~x, ~y) => ((~x)^2 + (~y)^2)^(1//2)
    
    expr_unwrapped = Symbolics.unwrap(expr)
    
    # Apply rewrites
    chain = SymbolicUtils.Chain([r_abs_sq, r_abs, r_hypot_sum, r_hypot_default])
    expr_rewritten = SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(chain))(expr_unwrapped)
    
    # Manual power fixing function
    function fix_powers(t)
        if Symbolics.iscall(t)
            op = Symbolics.operation(t)
            if op == (^)
                args = Symbolics.arguments(t)
                base = args[1]
                expon = args[2]
                
                if Symbolics.iscall(base) && Symbolics.operation(base) == (^)
                    base_args = Symbolics.arguments(base)
                    inner_base = base_args[1]
                    inner_expon = base_args[2]
                    
                    new_expon = inner_expon * expon
                    if isinteger(new_expon)
                        return inner_base^Int(new_expon)
                    else
                        return inner_base^new_expon
                    end
                end
                
                if expon isa Rational && isinteger(expon)
                     return base^Int(expon)
                end
            end
        end
        return t
    end
    
    expr_rewritten = SymbolicUtils.Postwalk(fix_powers)(expr_rewritten)
    expr = Num(expr_rewritten)
    
    # Expand again
    expr = Symbolics.expand(expr)

    # Substitute Re(U) and Im(U)
    subs_dict = Dict{Any, Any}()
    U_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    U_bar_lookup = Dict{Any, Tuple{Int, Int}}()
    
    if U_sym isa AbstractArray
        for i in 1:size(U_sym, 1)
            for j in 1:size(U_sym, 2)
                u_ij = U_sym[i,j]
                u_atomic = Symbolics.variable(:U_atomic, i, j)
                u_bar_atomic = Symbolics.variable(:U_bar_atomic, i, j)
                
                U_atomic_lookup[Symbolics.unwrap(u_atomic)] = (i, j)
                U_bar_lookup[Symbolics.unwrap(u_bar_atomic)] = (i, j)
                
                subs_dict[u_ij] = u_atomic
                subs_dict[conj(u_ij)] = u_bar_atomic
                subs_dict[real(u_ij)] = (1//2) * (u_atomic + u_bar_atomic)
                subs_dict[imag(u_ij)] = (-1//2) * IM_dummy * (u_atomic - u_bar_atomic)
            end
        end
    end
    
    # Substitute
    expr_subbed = Symbolics.substitute(expr, subs_dict)
    
    # Expand
    expanded_expr = Symbolics.expand(expr_subbed)

    # Substitute IM_dummy back to im
    # Use -1 for IM_dummy^2 first to simplify, then IM_dummy to im
    expanded_expr = Symbolics.substitute(expanded_expr, Dict(IM_dummy^2 => -1))
    expanded_expr = Symbolics.substitute(expanded_expr, Dict(IM_dummy => im))

    # Helper to traverse product
    function process_term_wrapped(term)
        return process_term(term, U_atomic_lookup, U_bar_lookup, dim)
    end

    # Local helper to sum integrals over terms
    function integrate_num_expr(ex)
        ex_un = Symbolics.unwrap(ex)
        if Symbolics.iscall(ex_un) && Symbolics.operation(ex_un) == (+)
            terms = Symbolics.arguments(ex_un)
            return sum(process_term_wrapped, terms)
        end
        return process_term_wrapped(ex)
    end

    # Handle Complex{Num} result from substitution
    final_res = integrate_num_expr(expanded_expr)
    return _robust_real(final_res)
end

function _robust_real(x)
    if x isa AbstractArray
        return map(_robust_real, x)
    end
    
    # Standardize to underlying Julia type if it's a constant
    unwrapped = Symbolics.unwrap(x)
    
    if unwrapped isa Complex
        return iszero(imag(unwrapped)) ? real(unwrapped) : unwrapped
    end
    
    # If it's a Num and symbolically real, return the real part
    try
        r = Symbolics.simplify(Symbolics.real(x))
        i = Symbolics.simplify(Symbolics.imag(x))
        if isequal(Symbolics.unwrap(i), 0)
            return Symbolics.unwrap(r)
        end
    catch
    end
    
    return unwrapped
end

function _iszero(x)
    x_un = Symbolics.unwrap(x)
    if x_un isa Number
        return iszero(x_un)
    end
    return isequal(x_un, 0)
end

function process_term(term, U_atomic_lookup, U_bar_lookup, dim)
    term = Symbolics.unwrap(term)
    
    if Symbolics.iscall(term)
        op = Symbolics.operation(term)
        if op == real
            return process_term(Symbolics.arguments(term)[1], U_atomic_lookup, U_bar_lookup, dim)
        elseif op == imag
            return 0
        elseif op == abs
            return process_term(Symbolics.arguments(term)[1], U_atomic_lookup, U_bar_lookup, dim)
        elseif op == (+)
            return sum(t -> process_term(t, U_atomic_lookup, U_bar_lookup, dim), Symbolics.arguments(term))
        end
    end

    coeff = 1 // 1
    u_indices = Vector{Tuple{Int, Int}}()
    u_bar_indices = Vector{Tuple{Int, Int}}()
    
    function traverse(t)
        t_unwrapped = Symbolics.unwrap(t)
        
        if haskey(U_atomic_lookup, t_unwrapped)
            push!(u_indices, U_atomic_lookup[t_unwrapped])
            return
        elseif haskey(U_bar_lookup, t_unwrapped)
            push!(u_bar_indices, U_bar_lookup[t_unwrapped])
            return
        end
        
        # Robust fallback
        for (k, v) in U_atomic_lookup
            if isequal(k, t_unwrapped) || string(k) == string(t_unwrapped)
                push!(u_indices, v)
                return
            end
        end
        for (k, v) in U_bar_lookup
            if isequal(k, t_unwrapped) || string(k) == string(t_unwrapped)
                push!(u_bar_indices, v)
                return
            end
        end
        
        if t isa Number && !(t isa Num) && !(t isa Complex{Num})
            coeff *= t
            return
        end
        
        if Symbolics.iscall(t_unwrapped)
            op = Symbolics.operation(t_unwrapped)
            args = Symbolics.arguments(t_unwrapped)
            
            if op == (*)
                for arg in args; traverse(arg); end
                return
            elseif op == (^)
                base = args[1]
                p_val = Symbolics.unwrap(args[2])
                p = try parse(Int, string(p_val)) catch; nothing end
                if p isa Integer
                    for _ in 1:p; traverse(base); end
                    return
                end
            elseif op == real || op == Base.real
                traverse(args[1]); return
            elseif op == imag || op == Base.imag
                coeff = 0; return
            elseif op == abs || op == Base.abs
                traverse(args[1]); return
            elseif op == conj || op == Base.conj
                 traverse(args[1]); return
            end
        end
        
        coeff *= t
    end
    traverse(term)
    
    n = length(u_indices)
    if n != length(u_bar_indices)
        return 0
    end
    if n == 0
        return coeff
    end
    
    return coeff * integrate_indices(u_indices, u_bar_indices, dim)
end


function integrate_indices(U_idxs::Vector{Tuple{Int, Int}}, U_bar_idxs::Vector{Tuple{Int, Int}}, dim)
    n = length(U_idxs)
    I = [x[1] for x in U_idxs]
    J = [x[2] for x in U_idxs]
    I_bar = [x[1] for x in U_bar_idxs]
    J_bar = [x[2] for x in U_bar_idxs]
    
    valid = get_matching_permutations(I, I_bar)
    valid_taus = get_matching_permutations(J, J_bar)
    
    if isempty(valid) || isempty(valid_taus)
        return 0
    end
    
    total = 0 // 1
    for sigma in valid
        for tau in valid_taus
            inv_tau = invperm(tau)
            P = [sigma[inv_tau[i]] for i in 1:n]
            # Cycle type of sigma * tau^-1
            cycle_type = get_cycle_type(P)
            wg_val = weingarten(cycle_type, dim)
            total += wg_val
        end
    end
    return total
end

function get_matching_permutations(target::Vector{Int}, source::Vector{Int})
    n = length(target)
    if n != length(source)
        return Vector{Vector{Int}}()
    end

    # Group indices of source and target by value
    source_groups = Dict{Int, Vector{Int}}()
    for (idx, val) in enumerate(source)
        push!(get!(source_groups, val, Int[]), idx)
    end
    
    target_groups = Dict{Int, Vector{Int}}()
    for (idx, val) in enumerate(target)
        push!(get!(target_groups, val, Int[]), idx)
    end

    # Check if value sets and counts match
    if length(source_groups) != length(target_groups)
        return Vector{Vector{Int}}()
    end
    for (val, t_idxs) in target_groups
        if !haskey(source_groups, val) || length(source_groups[val]) != length(t_idxs)
            return Vector{Vector{Int}}()
        end
    end

    # Generate permutations for each value group
    vals = collect(keys(target_groups))
    group_perms = [collect(permutations(source_groups[v])) for v in vals]
    
    res = Vector{Vector{Int}}()
    for combined_p in Iterators.product(group_perms...)
        p = zeros(Int, n)
        for (v_idx, v_perms) in enumerate(combined_p)
            t_idxs = target_groups[vals[v_idx]]
            for (i, t_idx) in enumerate(t_idxs)
                p[t_idx] = v_perms[i]
            end
        end
        push!(res, p)
    end
    
    return res
end

function get_cycle_type(p::Vector{Int})
    visited = falses(length(p))
    lengths = Int[]
    for i in 1:length(p)
        if !visited[i]
            curr = i
            len = 0
            while !visited[curr]
                visited[curr] = true
                curr = p[curr]
                len += 1
            end
            push!(lengths, len)
        end
    end
    sort!(lengths, rev=true)
    return lengths
end

end # module
