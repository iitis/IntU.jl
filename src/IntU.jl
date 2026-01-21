module IntU

using Symbolics
using Combinatorics
using LinearAlgebra
using SymbolicUtils

include("Weingarten.jl")
using .Weingarten

export integrate, dU

# Dummy type to represent the measure
struct HaarMeasure{T, N, D}
    U::AbstractArray{T, N}
    dim::D
end
dU(U::AbstractArray{T,N}, dim) where {T,N} = HaarMeasure{T,N,typeof(dim)}(U, dim)

"""
    integrate(expr, measure::HaarMeasure)
"""
function integrate(expr, measure::HaarMeasure)
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
    # This is necessary because Symbolics.expand doesn't always handle abs(sum)^2 well.
    r_abs_sq = @rule abs(~x)^2 => (~x) * conj(~x)
    r_abs = @rule abs(~x) => hypot(real(~x), imag(~x))
    
    # Rewrite hypot(x,y) -> ((x + i*y)(x - i*y))^(1//2) ONLY if x, y are Sums
    # This optimization is crucial for Tr(U) (Sums) to avoid expansion explosion.
    r_hypot_sum = @rule hypot(~x::is_add, ~y::is_add) => ((~x + im*~y) * (~x - im*~y))^(1//2)
    # Fallback/Default for non-sums:
    r_hypot_default = @rule hypot(~x, ~y) => ((~x)^2 + (~y)^2)^(1//2)
    
    expr_unwrapped = Symbolics.unwrap(expr)
    
    # Apply rewrites in sequence: abs^2 -> z*conj(z), then remaining abs -> hypot, then hypot -> pow
    chain = SymbolicUtils.Chain([r_abs_sq, r_abs, r_hypot_sum, r_hypot_default])
    expr_rewritten = SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(chain))(expr_unwrapped)
    
    # Manual power fixing function to handle nested and rational powers robustly
    function fix_powers(t)
        if Symbolics.iscall(t)
            op = Symbolics.operation(t)
            if op == (^)
                args = Symbolics.arguments(t)
                base = args[1]
                expon = args[2]
                
                # Case 1: Nested power ((x^p)^q) -> x^(p*q)
                # Note: Postwalk visits children first, so base is already processed
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
                
                # Case 2: Rational integer power (x^(2//1)) -> x^2
                if expon isa Rational && isinteger(expon)
                     return base^Int(expon)
                end
            end
        end
        return t
    end
    
    expr_rewritten = SymbolicUtils.Postwalk(fix_powers)(expr_rewritten)
    
    expr = Num(expr_rewritten)
    
    # Expand again to distribute the square
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
    
    # Also handle IM_dummy substitution immediately if possible
    expr_subbed = Symbolics.substitute(expr_subbed, Dict(IM_dummy => im))
    
    # Expand
    expanded_expr = Symbolics.expand(expr_subbed)

    # Substitute IM_dummy back to im (which makes it Complex{Num} potentially)
    expanded_expr = Symbolics.substitute(expanded_expr, Dict(IM_dummy^2 => -1))
    expanded_expr = Symbolics.substitute(expanded_expr, Dict(IM_dummy => im))

    # Helper to traverse product
    function process_term(term)
        term = Symbolics.unwrap(term)
        
        # Handle wrappers introduced by Symbolics (e.g. real(Sum), imag(Sum), abs(Sum))
        if Symbolics.iscall(term)
            op = Symbolics.operation(term)
            if op == real
                return process_term(Symbolics.arguments(term)[1])
            elseif op == imag
                return 0
            elseif op == abs
                # abs on sum might be tricky, but on monomials we unwrap it
                return process_term(Symbolics.arguments(term)[1])
            elseif op == (+)
                return sum(process_term, Symbolics.arguments(term))
            end
        end

        coeff = 1
        u_indices = Vector{Tuple{Int, Int}}()
        u_bar_indices = Vector{Tuple{Int, Int}}()
        
        function traverse(t)
            t_unwrapped = Symbolics.unwrap(t)
            
            # Check atomic lookups with unwrapped variable
            if haskey(U_atomic_lookup, t_unwrapped)
                push!(u_indices, U_atomic_lookup[t_unwrapped])
                return
            elseif haskey(U_bar_lookup, t_unwrapped)
                push!(u_bar_indices, U_bar_lookup[t_unwrapped])
                return
            end
            
            # Robust fallback using isequal and string comparison
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
            
            # Literal numbers
            if t isa Number && !(t isa Num) && !(t isa Complex{Num})
                coeff *= t
                return
            end
            
            if Symbolics.iscall(t_unwrapped)
                op = Symbolics.operation(t_unwrapped)
                args = Symbolics.arguments(t_unwrapped)
                
                op_str = string(op)
                
                if op_str == "*"
                    for arg in args
                        traverse(arg)
                    end
                    return
                elseif op_str == "^"
                    base = args[1]
                    p_val = Symbolics.unwrap(args[2])
                    p = try parse(Int, string(p_val)) catch; nothing end
                    if p isa Integer
                        for _ in 1:p
                            traverse(base)
                        end
                        return
                    end
                elseif op_str == "real" || op_str == "Base.real"
                    traverse(args[1])
                    return
                elseif op_str == "imag" || op_str == "Base.imag"
                    coeff = 0
                    return
                elseif op_str == "abs" || op_str == "Base.abs"
                    traverse(args[1])
                    return
                elseif op_str == "conj" || op_str == "Base.conj"
                     traverse(args[1])
                     return
                end
            end
            
            coeff *= t
        end
        traverse(term)
        
        n = length(u_indices)
        if n != length(u_bar_indices)
            return 0 # Orthogonality
        end
        if n == 0
            return coeff
        end
        
        return coeff * integrate_indices(u_indices, u_bar_indices, dim)
    end

    # Local helper to sum integrals over terms
    function integrate_num_expr(ex)
        ex_un = Symbolics.unwrap(ex)
        if Symbolics.iscall(ex_un) && Symbolics.operation(ex_un) == (+)
            terms = Symbolics.arguments(ex_un)
            return sum(process_term, terms)
        end
        return process_term(ex)
    end

    # Handle Complex{Num} by splitting
    if expanded_expr isa Complex
        re_part = real(expanded_expr)
        im_part = imag(expanded_expr)
        
        val_re = integrate_num_expr(re_part)
        val_im = integrate_num_expr(im_part)
        return val_re + im * val_im
    end

    return integrate_num_expr(expanded_expr)
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
    perms = permutations(1:n)
    valid = Vector{Vector{Int}}()
    for p in perms
        match = true
        for i in 1:n
            if target[i] != source[p[i]]
                match = false
                break
            end
        end
        if match; push!(valid, p); end
    end
    return valid
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
