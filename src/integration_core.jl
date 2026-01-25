# Core integration logic

# Resolve Num(::Complex) ambiguity
Symbolics.Num(z::Complex) = Complex(Num(real(z)), Num(imag(z)))
Base.isequal(a::Complex{Num}, b::Num) = isequal(real(a), b) && isequal(imag(a), 0)
Base.isequal(a::Num, b::Complex{Num}) = isequal(b, a)

function _integrate_core(expr, dim, subs_dict, U_atomic_lookup, U_bar_lookup)
    if expr isa Complex
        val_re = _integrate_core(real(expr), dim, subs_dict, U_atomic_lookup, U_bar_lookup)
        val_im = _integrate_core(imag(expr), dim, subs_dict, U_atomic_lookup, U_bar_lookup)
        return _robust_real(val_re + im * val_im)
    end

    expr_un = Symbolics.unwrap(expr)
    if expr_un isa Number && (expr_un isa Real || expr_un isa Complex)
        return expr_un
    end
    
    # Check if unwrapping revealed a complex object
    if expr_un isa Complex
        return _integrate_core(expr_un, dim, subs_dict, U_atomic_lookup, U_bar_lookup)
    end

    # Wrap in Num if not already, then expand
    expr_num = _safe_Num(expr_un)
    try
        expr_num = Symbolics.expand(expr_num)
    catch
    end

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
    
    # Rewrite complex(x, y) calls to x + im*y so expanding works
    r_complex = @rule complex(~x, ~y) => ~x + im*~y
    r_complex_base = @rule Base.complex(~x, ~y) => ~x + im*~y

    expr_unwrapped = Symbolics.unwrap(expr)
    
    # Apply rewrites
    chain = SymbolicUtils.Chain([r_abs_sq, r_abs, r_hypot_sum, r_hypot_default, r_complex, r_complex_base])
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
    if expr_rewritten isa Complex
        return _integrate_core(expr_rewritten, dim, subs_dict, U_atomic_lookup, U_bar_lookup)
    end
    expr_num = _safe_Num(expr_rewritten)
    
    # Expand again
    try
        expr_num = Symbolics.expand(_safe_Num(Symbolics.unwrap(expr_num)))
    catch
    end

    # Substitute using SymbolicUtils but fallback to Postwalk on TypeError
    function robust_substitute(ex, dict)
        try
            return Symbolics.substitute(Symbolics.wrap(ex), dict)
        catch
            # Truly brute force fallback using only symbols if possible
            return SymbolicUtils.Postwalk(x -> x isa Symbolics.Num ? get(dict, x, x) : x)(ex)
        end
    end

    expr_subbed = if expr_num isa Complex
        robust_substitute(Symbolics.unwrap(real(expr_num)), subs_dict) + im * robust_substitute(Symbolics.unwrap(imag(expr_num)), subs_dict)
    else
        robust_substitute(Symbolics.unwrap(expr_num), subs_dict)
    end
    
    # Expand
    expanded_expr = try
        Symbolics.expand(_safe_Num(Symbolics.unwrap(expr_subbed)))
    catch
        _safe_Num(expr_subbed)
    end
    
    # Apply rewrites again to catch complex(...) introduced by substitution
    expanded_expr = SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(chain))(Symbolics.unwrap(expanded_expr))
    # Safe wrap and expand again to distribute
    try
        expanded_expr = Symbolics.expand(_safe_Num(expanded_expr))
    catch e
        println("DEBUG: Expansion failed: ", e)
        expanded_expr = _safe_Num(expanded_expr)
    end


    # Helper to traverse product
    function process_term_wrapped(term)
        return process_term(term, U_atomic_lookup, U_bar_lookup, dim)
    end

    # Local helper to sum integrals over terms
    function integrate_num_expr(ex)
        ex_un = Symbolics.unwrap(ex)
        if ex_un isa Complex
             return integrate_num_expr(real(ex_un)) + im * integrate_num_expr(imag(ex_un))
        end
        if ex_un isa Number && (ex_un isa Real)
             return ex_un
        end
        if Symbolics.iscall(ex_un)
            op = Symbolics.operation(ex_un)
            if op == (+)
                terms = Symbolics.arguments(ex_un)
                return sum(process_term_wrapped, terms)
            end
            if op == complex || op == Base.complex
                args = Symbolics.arguments(ex_un)
                return integrate_num_expr(args[1]) + im * integrate_num_expr(args[2])
            end
        end
        return process_term_wrapped(ex)
    end

    # Handle result
    final_res = integrate_num_expr(expanded_expr)
    return _robust_real(final_res)
end

function _safe_Num(x)
    if x isa Num || x isa Complex{Num}
        return x
    end
    # Check if it's a raw symbolic object from Symbolics/SymbolicUtils
    if !(x isa Number) && !(x isa AbstractArray)
        # If it's not a standard number/array, it's likely a symbolic object.
        # We can try to wrap it in Num.
        return Num(x)
    end
    x_un = Symbolics.unwrap(x)
    if x_un isa Complex
        return Complex(Num(real(x_un)), Num(imag(x_un)))
    end
    return Num(x_un)
end

function _robust_real(x)
    if x isa AbstractArray
        return map(_robust_real, x)
    end
    
    x_un = Symbolics.unwrap(x)

    # If it's already a plain Complex, check if it's real
    if x_un isa Complex
        if _iszero(imag(x_un))
            return _robust_real(real(x_un))
        end
        return Complex(_robust_real(real(x_un)), _robust_real(imag(x_un)))
    end

    # If it's a plain Real, we're done
    if x_un isa Real
        return x_un
    end

    # Try symbolic simplification to see if it's real
    try
        nx = _safe_Num(x_un)
        ix = Symbolics.simplify(imag(nx))
        if _iszero(ix)
            rx = Symbolics.simplify(real(nx))
            ux = Symbolics.unwrap(rx)
            # If unwrapped part is a real number, return it
            if ux isa Real
                return ux
            end
            return ux
        end
    catch
    end
    
    return x_un
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
        
        t_str = string(t_unwrapped)
        
        # Exact match or string match for robustness
        found = false
        for (k, v) in U_atomic_lookup
            if isequal(k, t_unwrapped) || string(k) == t_str
                push!(u_indices, v)
                found = true; break
            end
        end
        found && return
        
        for (k, v) in U_bar_lookup
            if isequal(k, t_unwrapped) || string(k) == t_str
                push!(u_bar_indices, v)
                found = true; break
            end
        end
        found && return
        
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
    
    val = integrate_indices(u_indices, u_bar_indices, dim)
    if isequal(val, 0)
        return 0
    end
    return coeff * val
end


"""
    integrate_indices(U_idxs, U_bar_idxs, dim)

Low-level integration function using Weingarten calculus.
`U_idxs` is a list of (i, j) tuples for each U_{i,j} in the product.
`U_bar_idxs` is a list of (k, l) tuples for each \bar{U}_{k,l}.
"""
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

# Overloads for symbolic types to make them "just work" with tr, det etc.
function LinearAlgebra.tr(A::Symbolics.Arr{T, 2}) where T
    return sum(A[i,i] for i in 1:size(A,1))
end

"""
    _poly_degree(p, d)

Helper to get the degree of a polynomial `p` in variable `d`.
"""
function _poly_degree(p, d)
    p_un = Symbolics.unwrap(p)
    if isequal(p_un, d) return 1 end
    if !(p_un isa SymbolicUtils.BasicSymbolic) return 0 end
    if Symbolics.iscall(p_un)
        op = Symbolics.operation(p_un)
        args = Symbolics.arguments(p_un)
        if op == (+)
            return maximum(a -> _poly_degree(a, d), args, init=0)
        elseif op == (*)
            return sum(a -> _poly_degree(a, d), args)
        elseif op == (^)
            # Handle power if exponent is an integer
            base = args[1]
            expon = Symbolics.unwrap(args[2])
            if expon isa Integer
                return _poly_degree(base, d) * expon
            elseif expon isa Rational && isinteger(expon)
                return _poly_degree(base, d) * Int(expon)
            end
        end
    end
    return 0
end


"""
    _expand_asymptotic(ex, d, order)

Helper to expand a rational function of `d` in powers of `1/d`.
Uses Taylor expansion in ε = 1/d around ε = 0.
"""
function _expand_asymptotic(ex, d, order)
    ex_un = Symbolics.unwrap(ex)
    if ex_un isa AbstractArray
        return map(e -> _expand_asymptotic(e, d, order), ex_un)
    end
    if ex_un isa Complex
        re_part = _expand_asymptotic(real(ex_un), d, order)
        im_part = _expand_asymptotic(imag(ex_un), d, order)
        return re_part + im * im_part
    end
    if !(ex_un isa Symbolics.Num) && !(ex_un isa SymbolicUtils.BasicSymbolic)
        return ex_un
    end
    
    # We want expansion in 1/d. Let eps = 1/d.
    # Ensure it is a single fraction
    ex_sim = SymbolicUtils.simplify_fractions(ex_un)
    num = Symbolics.numerator(ex_sim)
    den = Symbolics.denominator(ex_sim)
    
    ϵ = Symbolics.variable(:ϵ)
    
    # Get degrees to clear denominators
    n = _poly_degree(num, d)
    m = _poly_degree(den, d)
    max_deg = max(n, m)
    
    # Substitute d -> 1/ϵ and clear denominators by multiplying by ϵ^max_deg
    # P_eps(ϵ) = num(1/ϵ) * ϵ^max_deg
    # Q_eps(ϵ) = den(1/ϵ) * ϵ^max_deg
    p_eps = Symbolics.simplify(Symbolics.substitute(num, Dict(d => 1/ϵ)) * ϵ^max_deg)
    q_eps = Symbolics.simplify(Symbolics.substitute(den, Dict(d => 1/ϵ)) * ϵ^max_deg)
    
    # Now we have a well-behaved rational function f_eps = p_eps / q_eps
    f_eps = p_eps / q_eps
    
    total = 0 // 1
    curr_deriv = f_eps
    
    for k in 0:order
        # Evaluate at ϵ = 0. Substitute into p_eps and q_eps separately to avoid 0/0 if possible,
        # but p_eps/q_eps should be fine after simplification.
        val = try
            # If q_eps(0) is 0, then f_eps is singular at 0. 
            # This would mean the integral has positive powers of d.
            # We handle this by using limit or just substitute.
            Symbolics.substitute(curr_deriv, Dict(ϵ => 0))
        catch
            Symbolics.substitute(Symbolics.simplify(curr_deriv), Dict(ϵ => 0))
        end
        
        if !isequal(val, 0)
            term = (val // factorial(k)) * (1/d)^k
            total += term
        end
        
        if k < order
            curr_deriv = Symbolics.derivative(curr_deriv, ϵ)
        end
    end
    
    return Symbolics.simplify(total)
end
