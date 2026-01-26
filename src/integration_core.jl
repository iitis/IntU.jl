# Core integration logic

# Helper to wrap complex numbers in Num safely
_to_Num(z::Complex) = Complex(Num(real(z)), Num(imag(z)))
_to_Num(z) = Num(z)
Base.isequal(a::Complex{Num}, b::Num) = isequal(real(a), b) && isequal(imag(a), 0)
Base.isequal(a::Num, b::Complex{Num}) = isequal(b, a)

function _integrate_core(expr, dim, subs_dict, U_atomic_lookup, U_bar_lookup, measure_type=:U)
    if expr isa Complex
        val_re = _integrate_core(real(expr), dim, subs_dict, U_atomic_lookup, U_bar_lookup, measure_type)
        val_im = _integrate_core(imag(expr), dim, subs_dict, U_atomic_lookup, U_bar_lookup, measure_type)
        return _robust_real(val_re + im * val_im)
    end

    expr_un = Symbolics.unwrap(expr)
    if expr_un isa Number && (expr_un isa Real || expr_un isa Complex)
        return expr_un
    end
    
    # Check if unwrapping revealed a complex object
    if expr_un isa Complex
        return _integrate_core(expr_un, dim, subs_dict, U_atomic_lookup, U_bar_lookup, measure_type)
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
    
    r_real = @rule real(~x) => (1//2) * (~x + conj(~x))
    r_real_base = @rule Base.real(~x) => (1//2) * (~x + conj(~x))
    r_imag = @rule imag(~x) => (1//(2im)) * (~x - conj(~x))
    r_imag_base = @rule Base.imag(~x) => (1//(2im)) * (~x - conj(~x))
    
    # Rewrite hypot(x,y) -> ((x + i*y)(x - i*y))^(1//2) ONLY if x, y are Sums
    r_hypot_sum = @rule hypot(~x::is_add, ~y::is_add) => ((~x + im*~y) * (~x - im*~y))^(1//2)
    # Fallback/Default for non-sums:
    r_hypot_default = @rule hypot(~x, ~y) => ((~x)^2 + (~y)^2)^(1//2)
    
    # Rewrite complex(x, y) calls to x + im*y so expanding works
    r_complex = @rule complex(~x, ~y) => ~x + im*~y
    r_complex_base = @rule Base.complex(~x, ~y) => ~x + im*~y

    expr_unwrapped = Symbolics.unwrap(expr)
    
    # Apply rewrites
    chain = SymbolicUtils.Chain([r_abs_sq, r_abs, r_real, r_real_base, r_imag, r_imag_base, r_hypot_sum, r_hypot_default, r_complex, r_complex_base])
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
        return _integrate_core(expr_rewritten, dim, subs_dict, U_atomic_lookup, U_bar_lookup, measure_type)
    end
    expr_num = _safe_Num(expr_rewritten)
    
    # Expand again
    try
        expr_num = Symbolics.expand(_safe_Num(Symbolics.unwrap(expr_num)))
    catch
    end

    # Substitute using Symbolics.substitute with a robust Postwalk fallback
    function robust_substitute(ex, dict)
        # 1. Try high-level substitute (works for Num and mapped sub-expressions)
        try
            res = Symbolics.substitute(Symbolics.wrap(ex), dict)
            # If nothing changed, try deeper dive
            if isequal(res, Symbolics.wrap(ex))
                 throw(error("No change"))
            end
            return res
        catch
            # 2. Fallback to manual Postwalk traversal
            # We unwrap both the expression and the dictionary keys for raw comparison
            unwrapped_dict = Dict(Symbolics.unwrap(k) => Symbolics.unwrap(v) for (k, v) in dict)
            
            p_res = SymbolicUtils.Postwalk(x -> begin
                u = Symbolics.unwrap(x)
                if haskey(unwrapped_dict, u)
                    return unwrapped_dict[u]
                end
                return x
            end)(Symbolics.unwrap(ex))
            
            return Symbolics.wrap(p_res)
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
        return process_term(term, U_atomic_lookup, U_bar_lookup, dim, measure_type)
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
    return _to_Num(x_un)
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

function process_term(term, U_atomic_lookup, U_bar_lookup, dim, measure_type=:U)
    term = Symbolics.unwrap(term)
    
    if Symbolics.iscall(term)
        op = Symbolics.operation(term)
        if op == (+)
            return sum(t -> process_term(t, U_atomic_lookup, U_bar_lookup, dim, measure_type), Symbolics.arguments(term))
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
                # Check if it's a conjugate entry for Symplectic
                if length(v) == 3 && v[3] == :conj
                    if measure_type == :Sp
                        if !(dim isa Integer)
                             error("Symplectic integration with conjugates requires integer dimension")
                        end
                        i_idx, j_idx, _ = v
                        m_dim = div(dim, 2)
                        
                        di, si = (i_idx <= m_dim ? (i_idx + m_dim, 1) : (i_idx - m_dim, -1))
                        dj, sj = (j_idx <= m_dim ? (j_idx + m_dim, 1) : (j_idx - m_dim, -1))
                        
                        coeff *= (si * sj)
                        push!(u_indices, (di, dj))
                    else
                        error("Conjugate flag :conj found for non-Symplectic measure")
                    end
                else
                    push!(u_indices, v)
                end
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
            elseif op == conj || op == Base.conj
                inner = Symbolics.unwrap(args[1])
                inner_str = string(inner)
                # conj(U) -> U_bar
                for (k, v) in U_atomic_lookup
                    if isequal(k, inner) || string(k) == inner_str
                        if measure_type == :Sp
                             # Perform same Sp duality rewrite here!
                             if !(dim isa Integer)
                                  error("Symplectic integration with conjugates requires integer dimension")
                             end
                             i_i, j_j = v
                             m_m = div(dim, 2)
                             di, si = (i_i <= m_m ? (i_i + m_m, 1) : (i_i - m_m, -1))
                             dj, sj = (j_j <= m_m ? (j_j + m_m, 1) : (j_j - m_m, -1))
                             coeff *= (si * sj)
                             push!(u_indices, (di, dj))
                        else
                             push!(u_bar_indices, v)
                        end
                        return
                    end
                end
                # conj(U_bar) -> U
                for (k, v) in U_bar_lookup
                    if isequal(k, inner) || string(k) == inner_str
                        push!(u_indices, v)
                        return
                    end
                end
            end
        end
        
        coeff *= t
    end
    traverse(term)
    
    n_u = length(u_indices)
    n_bar = length(u_bar_indices)

    if measure_type == :U
        if n_u != n_bar
            return 0
        end
        if n_u == 0
            return coeff
        end
        
        val = integrate_indices(u_indices, u_bar_indices, dim)
        if isequal(val, 0)
            return 0
        end
        return coeff * val
        
    elseif measure_type == :O
        # Orthogonal measure: O combined with O_bar (which is same as O)
        # We combine both lists
        all_indices = [u_indices; u_bar_indices]
        n_total = length(all_indices)
        
        if n_total % 2 != 0
            return 0
        end
        if n_total == 0
            return coeff
        end
        
        val = integrate_indices_orthogonal(all_indices, dim)
        if isequal(val, 0)
            return 0
        end
        return coeff * val
        
    elseif measure_type == :Sp
        # Symplectic measure
        # We combine both lists (u_indices and u_bar_indices are treated same for real/symplectic generally,
        # but technically S_bar = J S J^T ? No S is unitary symplectic.
        # S_ij is treated as variable.
        # Formula uses S_{i_1 j_1} ... S_{i_2k j_2k}.
        all_indices = [u_indices; u_bar_indices]
        n_total = length(all_indices)
        
        if n_total % 2 != 0
            return 0
        end
        if n_total == 0
            return coeff
        end
        
        val = integrate_indices_symplectic(all_indices, dim)
        if isequal(val, 0)
            return 0
        end
        return coeff * val
    else
        error("Unknown measure type: $measure_type")
    end
end


"""
    integrate_indices(U_idxs, U_bar_idxs, dim)

Low-level integration function using Weingarten calculus (Unitary).
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

"""
    integrate_indices_orthogonal(indices, dim)
    
Low-level integration function using Orthogonal Weingarten calculus.
Indices are a list of (i, j) for O_{ij}.
Formula: sum_{pi, sigma in PairPartitions} delta_pi(i) * delta_sigma(j) * Wg(pi, sigma)
"""
function integrate_indices_orthogonal(indices::Vector{Tuple{Int, Int}}, dim)
    n = length(indices) # This is 2k
    if n % 2 != 0; return 0; end
    
    I = [x[1] for x in indices]
    J = [x[2] for x in indices]
    
    valid_pi = get_matching_pair_partitions_filtered(I)
    valid_sigma = get_matching_pair_partitions_filtered(J)
    
    if isempty(valid_pi) || isempty(valid_sigma)
        return 0
    end
    
    total = 0 // 1
    for pi in valid_pi
        for sigma in valid_sigma
            # Wg(pi, sigma)
            val = weingarten_orthogonal_val(pi, sigma, dim)
            total += val
        end
    end
    return total
end

"""
    get_matching_pair_partitions_filtered(indices)

Returns all pair partitions of 1..2k such that for every pair (a,b), indices[a] == indices[b].
"""
function get_matching_pair_partitions_filtered(indices::Vector{Int})
    n = length(indices)
    
    # Validation: each value must appear an even number of times.
    counts = Dict{Int, Vector{Int}}()
    for (pos, val) in enumerate(indices)
        push!(get!(counts, val, Int[]), pos)
    end
    
    for (val, pos_list) in counts
        if length(pos_list) % 2 != 0
            return Vector{Vector{Tuple{Int, Int}}}()
        end
    end
    
    # For each group of positions, generate all pair partitions
    # Then take Cartesian product
    
    # We can reuse get_pair_partitions but we need to map indices back to positions
    
    group_partitions = Vector{Vector{Vector{Tuple{Int, Int}}}}()
    
    for (val, pos_list) in counts
        k_local = length(pos_list)
        # Generate partitions of 1..k_local
        local_parts_idx = get_pair_partitions(k_local)
        
        # Map back to real positions
        real_parts = Vector{Vector{Tuple{Int, Int}}}()
        for p in local_parts_idx
            mapped = [(pos_list[a], pos_list[b]) for (a,b) in p]
            push!(real_parts, mapped)
        end
        push!(group_partitions, real_parts)
    end
    
    # Combine
    res = Vector{Vector{Tuple{Int, Int}}}()
    for combined in Iterators.product(group_partitions...)
        # flattened list of pairs
        full_part = Vector{Tuple{Int, Int}}()
        for part in combined
            append!(full_part, part)
        end
        # Sort internal pairs and list for consistency (optional but good for debugging)
        push!(res, full_part)
    end
    
    return res
end


"""
    integrate_indices_symplectic(indices, dim)

Low-level integration function using Symplectic Weingarten calculus.
Formula: sum_{pi, sigma} J_pi(i) * J_sigma(j) * Wg^Sp(pi, sigma)
"""
function integrate_indices_symplectic(indices::Vector{Tuple{Int, Int}}, dim)
    n = length(indices)
    k = div(n, 2)
    partitions = get_pair_partitions(n)
    
    I = [x[1] for x in indices]
    J_idx = [x[2] for x in indices]
    
    # Check if dim is clearly integer for simplified J evaluation
    dim_int = dim isa Integer ? dim : nothing
    
    total = 0 // 1
    
    # Optimization: Filter partitions that yield non-zero J-contraction first?
    # J_{uv} is non-zero only if |u - v| = d/2 (appropriately signed).
    # This requires looking at the actual index values.
    
    # Pre-calculate contractions
    # For each partition pi, calculate J_pi(I).
    # If I contains symbolic variables, this might return a symbolic expressions.
    # Currently we might assume indices are integers.
    
    pi_contractions = Dict{Any, Any}()
    sigma_contractions = Dict{Any, Any}()
    
    for pi in partitions
        val_I = compute_symplectic_contraction(pi, I, dim)
        if !isequal(val_I, 0)
             pi_contractions[pi] = val_I
        end
        
        val_J = compute_symplectic_contraction(pi, J_idx, dim)
        if !isequal(val_J, 0)
             sigma_contractions[pi] = val_J
        end
    end
    
    if isempty(pi_contractions) || isempty(sigma_contractions)
        return 0
    end
    
    for (pi, val_pi) in pi_contractions
        for (sigma, val_sigma) in sigma_contractions
            wg = weingarten_symplectic_val(pi, sigma, dim)
            total += val_pi * val_sigma * wg
        end
    end
    
    return total
end

"""
    compute_symplectic_contraction(partition, indices, dim)

Computes prod_{(u, v) in partition} J(indices[u], indices[v]).
J = [0 I; -I 0].
"""
function compute_symplectic_contraction(partition, indices, dim)
    val = 1
    
    # If dim is unknown, we can't evaluate J numerically easily unless indices are 
    # constructed such that we know.
    # But usually dim is symbolic `d`.
    # AND indices are something like 1, 2, ...
    # We need to know `d` to know if index i and j are symplectic pairs.
    # If `d` is symbolic, we can't determine this range.
    
    # However, usually when integrating, `dim` matches the matrix size indices.
    # If the user defines `dSp(S, d)` and `S` is 2x2. Then `d=2`.
    # If `S` is symbolic size, indices are typically symbolic too?
    # For now, we assume numeric indices and try to infer symplectic structure.
    # Or we assume d is the numeric size involved.
    
    # Fallback: if d is symbolic, we check if we can convert it to number if indices are numbers.
    # If indices are huge integers, we assume standard basis.
    
    # CRITICAL: We need `d` to define J.
    # If d is symbolic, we can try to look at max index?
    # No, that's unsafe.
    
    is_d_numeric = (dim isa Integer)
    
    for (u, v) in partition
        idx_u = indices[u]
        idx_v = indices[v]
        
        j_val = symplectic_form(idx_u, idx_v, dim)
        if isequal(j_val, 0)
            return 0
        end
        val *= j_val
    end
    return val
end

function symplectic_form(i, j, dim)
    # J matrix evaluation
    # returns 1 if j = i + m, -1 if j = i - m, 0 otherwise
    # where m = dim / 2.
    
    # If i, j are symbolic, return symbolic representation?
    # For now, simplistic implementation for numeric indices.
    
    if !(i isa Integer) || !(j isa Integer)
        # Symbolic indices not fully supported yet for contraction
        return 0 
    end
    
    # We require dim to be known integer for this check 
    # Or we return 0 if we can't check.
    if !(dim isa Integer)
         # Try to resolve or error?
         # If the user provided d symbolic but indices 1,2, then we can likely assume d=2 if max index is 2?
         # Unsafe.
         # For now, return 0 or error.
         # Actually, better to error if we can't compute.
         # But usually we want to return 0 for non-matching.
         # Let's try to assume small dimension if indices are small?
         # Re-eval plan: "Support numeric indices". 
         # We need dim to be numeric.
         return 0 
         # TODO: Support J(i, j) symbolic object.
    end
    
    m = div(dim, 2)
    
    # 1 <= i, j <= 2m
    if i < 1 || i > 2*m || j < 1 || j > 2*m
        return 0
    end
    
    if j == i + m
        return 1
    elseif j == i - m
        return -1
    else
        return 0
    end
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
    # Use Symbolics' built-in degree function for maximum reliability
    try
        # Symbolics.degree(poly, var)
        # We need to make sure both are wrapped
        deg = Symbolics.degree(Symbolics.wrap(p), Symbolics.wrap(d))
        return Int(Symbolics.unwrap(deg))
    catch
        # Fallback to very basic manual checking if degree fails
        p_un = Symbolics.unwrap(p)
        d_un = Symbolics.unwrap(d)
        if isequal(p_un, d_un) return 1 end
        return 0
    end
end


"""
    asymptotic(ex, d, order=1)

Generic asymptotic expansion of a rational function `ex` in powers of `1/d`.
"""
function asymptotic(ex, d, order=1)
    return _expand_asymptotic(ex, d, order)
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
    if !(ex_un isa Symbolics.Num) && !(ex_un isa SymbolicUtils.BasicSymbolic) && !(ex_un isa Number)
        return ex_un
    end
    
    d_un = Symbolics.unwrap(d)
    d_num = Symbolics.wrap(d_un)
    
    # We want expansion in 1/d. Let eps = 1/d.
    # To handle rational functions correctly:
    # Get numerator and denominator after simplification to ensure rational form
    ex_sim = Symbolics.simplify(Symbolics.wrap(ex_un))
    
    # Handle the case where simplify result is still complex (e.g. constant complex or complex num)
    if ex_sim isa Complex
        re_part = _expand_asymptotic(real(ex_sim), d, order)
        im_part = _expand_asymptotic(imag(ex_sim), d, order)
        return re_part + im * im_part
    end
    
    num = Symbolics.numerator(ex_sim)
    den = Symbolics.denominator(ex_sim)
    
    # Get degrees to clear internal denominators when d -> 1/eps
    n = _poly_degree(num, d_un)
    m = _poly_degree(den, d_un)
    
    ϵ = Symbolics.variable(:ϵ)
    ϵ_un = Symbolics.unwrap(ϵ)
    
    # If n > m, we have positive powers of d (pole at infinity).
    # If m >= n, we have rational part starting from (1/d)^(m-n).
    
    # f_analytic = ex(1/eps) * eps^(n-m)
    # Use expand=true to ensure internal 1/ϵ terms are cleared immediately
    p_eps = Symbolics.simplify(Symbolics.substitute(num, Dict(d_num => 1/ϵ)) * ϵ^n; expand=true)
    q_eps = Symbolics.simplify(Symbolics.substitute(den, Dict(d_num => 1/ϵ)) * ϵ^m; expand=true)
    
    # Now (p_eps / q_eps) is analytic at 0.
    f_analytic = Symbolics.simplify(p_eps / q_eps; expand=true)
    
    total = Num(0)
    curr_deriv = f_analytic
    # ex(1/eps) = f_analytic(eps) * eps^(m-n)
    diff = m - n
    
    # We need to expand until we reach requested order in 1/d
    # power = k + diff. We want power <= order.
    # So k <= order - diff.
    # But if diff is negative (polynomial), we want to capture all terms.
    max_k = order - diff
    
    for k in 0:max_k
        # Robust evaluation at ϵ = 0
        # Aggressively simplify to clear denominators before substitution
        val = try
            v = Symbolics.substitute(curr_deriv, Dict(ϵ => 0))
            vu = Symbolics.unwrap(v)
            if isnan(vu) || isinf(vu)
                # If NaN/Inf, try more aggressive simplification
                curr_sim = Symbolics.simplify(curr_deriv; expand=true)
                Symbolics.substitute(curr_sim, Dict(ϵ => 0))
            else
                v
            end
        catch
            # Fallback to direct replacement in expression (Postwalk-like recursion)
            try
                curr_sim = Symbolics.simplify(curr_deriv; expand=true)
                Symbolics.substitute(curr_sim, Dict(ϵ => 0))
            catch
                # Last resort: very basic search and replace
                Symbolics.wrap(SymbolicUtils.Postwalk(x -> isequal(x, ϵ_un) ? 0 : x)(Symbolics.unwrap(curr_deriv)))
            end
        end
        
        if !_iszero(val)
            power = k + diff
            term = (val * (1//factorial(k))) * (1/d_num)^power
            total += term
        end
        
        if k < max_k
            curr_deriv = Symbolics.derivative(curr_deriv, ϵ)
            # Occasional simplify to manage expression size
            if k % 2 == 0
                 try
                     curr_deriv = Symbolics.simplify(curr_deriv; expand=true)
                 catch
                     # If simplify fails, continue with unsimplified expression
                 end
            end
        end
    end
    
    return Symbolics.simplify(total)
end

