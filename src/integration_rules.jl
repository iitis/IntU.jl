
"""
    _dispatch_paired(u, ub, d, integrate_fn)

Helper for rules where U and Ū indices must pair up (|u| == |ub|).
Used by :U, :COE, :GinUE, etc.
"""
function _dispatch_paired(u, ub, d, integrate_fn)
    d = _ensure_symbolic_dim(d)
    length(u) != length(ub) ? 0 : (length(u) == 0 ? 1 : integrate_fn(u, ub, d))
end

"""
    _dispatch_merged(u, ub, d, integrate_fn)

Helper for rules where U and Ū indices are merged into one list
and the total count must be even. Used by :O, :GUE, :GOE, :GSE, etc.
"""
function _dispatch_merged(u, ub, d, integrate_fn)
    d = _ensure_symbolic_dim(d)
    all_indices = [u; ub]
    length(all_indices) % 2 != 0 ? 0 :
    (length(all_indices) == 0 ? 1 : integrate_fn(all_indices, d))
end

INTEGRATION_RULES[:U] = (u, ub, d, mt) -> _dispatch_paired(u, ub, d, integrate_indices)

INTEGRATION_RULES[:V] =
    (u, ub, d, mt) -> begin
        k = mt isa Tuple ? mt[2] : d

        function is_out(idx)
            j = idx[2]
            j === nothing && return false
            ju = Symbolics.unwrap(j)
            ku = Symbolics.unwrap(k)
            if ju isa Integer && ku isa Integer
                return ju > ku
            end
            diff = Symbolics.simplify(ju - ku)
            diff_u = Symbolics.unwrap(diff)
            if diff_u isa Number && real(diff_u) > 0
                return true
            end
            return false
        end

        if any(is_out, u) || any(is_out, ub)
            return 0
        end
        d = _ensure_symbolic_dim(d)
        length(u) != length(ub) ? 0 : (length(u) == 0 ? 1 : integrate_indices(u, ub, d))
    end
INTEGRATION_RULES[:psi] =
    (u, ub, d, mt) -> begin
        if any(x -> !_symbolic_isequal(x[2], 1), u) || any(x -> !_symbolic_isequal(x[2], 1), ub)
            return 0
        end
        d = _ensure_symbolic_dim(d)
        length(u) != length(ub) ? 0 : (length(u) == 0 ? 1 : integrate_indices(u, ub, d))
    end

INTEGRATION_RULES[:O] =
    (u, ub, d, mt) -> _dispatch_merged(u, ub, d, integrate_indices_orthogonal)

INTEGRATION_RULES[:Sp] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        d_un = Symbolics.unwrap(d)
        if d_un isa Number && isinteger(d_un)
            d_int = Int(d_un)
            if isodd(d_int)
                return 0
            end
            m = div(d_int, 2)
            sign_factor = 1
            converted_ub = Vector{Tuple{Any,Any}}()

            for (p, q) in ub
                pk = p <= m ? p + m : p - m
                qk = q <= m ? q + m : q - m
                if (p <= m && q > m) || (p > m && q <= m)
                    sign_factor *= -1
                end
                push!(converted_ub, (pk, qk))
            end

            all_indices = [u; converted_ub]
            length(all_indices) % 2 != 0 ? 0 :
            (
                length(all_indices) == 0 ? 1 :
                sign_factor * integrate_indices_symplectic(all_indices, d)
            )
        elseif d_un isa Number
            throw(
                ArgumentError(
                    "Symplectic integration requires integer dimension d, got $d_un",
                ),
            )
        else
            m = d / 2
            sign_factor = 1
            converted_ub = Vector{Tuple{Any,Any}}()

            for (p, q) in ub
                pk = p + m
                qk = q + m
                push!(converted_ub, (pk, qk))
            end

            all_indices = [u; converted_ub]
            length(all_indices) % 2 != 0 ? 0 :
            (
                length(all_indices) == 0 ? 1 :
                sign_factor * integrate_indices_symplectic(all_indices, d)
            )
        end
    end

INTEGRATION_RULES[:GUE] =
    (u, ub, d, mt) -> _dispatch_merged(u, ub, d, integrate_indices_gue)

INTEGRATION_RULES[:GOE] =
    (u, ub, d, mt) -> _dispatch_merged(u, ub, d, integrate_indices_goe)

INTEGRATION_RULES[:GSE] =
    (u, ub, d, mt) -> _dispatch_merged(u, ub, d, integrate_indices_gse)

INTEGRATION_RULES[:COE] =
    (u, ub, d, mt) -> _dispatch_paired(u, ub, d, integrate_indices_coe)

INTEGRATION_RULES[:CSE] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        d_un = Symbolics.unwrap(d)
        length(u) != length(ub) ? 0 :
        (
            begin
                phys_dim = mt isa Tuple ? mt[2] : d
                length(u) == 0 ? 1 : integrate_indices_cse(u, ub, d, phys_dim)
            end
        )
    end

INTEGRATION_RULES[:Perm] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        all_indices = [u; ub]
        length(all_indices) == 0 ? 1 : integrate_indices_permutation(all_indices, d)
    end

INTEGRATION_RULES[:CPerm] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        all_indices = [u; ub]
        length(all_indices) == 0 ? 1 :
        integrate_indices_centered_permutation(all_indices, d)
    end

INTEGRATION_RULES[:Design] =
    (u, ub, d, mt) -> begin
        _, t_val = mt
        if length(u) != length(ub) || length(u) > t_val
            length(u) != length(ub) ? 0 :
            throw(
                ArgumentError(
                    "Integrand degree ($(length(u))) exceeds design order t=$t_val",
                ),
            )
        end
        length(u) == 0 ? 1 : integrate_indices(u, ub, d)
    end

INTEGRATION_RULES[:DiagUnitary] =
    (u, ub, d, mt) -> begin
        if length(u) != length(ub) || any(x -> x[1] != x[2], u) || any(x -> x[1] != x[2], ub)
            return 0
        end
        length(u) == 0 ? 1 : integrate_indices_diagonal(u, ub, d)
    end

INTEGRATION_RULES[:GinUE] =
    (u, ub, d, mt) -> _dispatch_paired(u, ub, d, integrate_indices_ginue)

INTEGRATION_RULES[:GinOE] =
    (u, ub, d, mt) -> _dispatch_merged(u, ub, d, integrate_indices_ginoe)

INTEGRATION_RULES[:GinSE] =
    (u, ub, d, mt) -> _dispatch_merged(u, ub, d, integrate_indices_ginse)

function _unitary_single_entry_moment(n::Int, dim)
    n == 0 && return one(Num)

    if dim isa Integer
        denom = one(BigInt)
        for j = 0:(n-1)
            denom *= (BigInt(dim) + j)
        end
        return factorial(big(n)) // denom
    end

    denom = one(Num)
    for j = 0:(n-1)
        denom *= (dim + j)
    end
    return factorial(big(n)) / denom
end

function _count_values(values)
    counts = Dict{Any,Int}()
    for v in values
        counts[v] = get(counts, v, 0) + 1
    end
    return counts
end

function _unitary_row_col_abs_moment(counts::Dict{Any,Int}, dim)
    n = sum(values(counts))
    n == 0 && return one(Num)

    num = one(BigInt)
    for c in values(counts)
        num *= factorial(big(c))
    end

    if dim isa Integer
        denom = one(BigInt)
        for j = 0:(n-1)
            denom *= (BigInt(dim) + j)
        end
        return num // denom
    end

    denom = one(Num)
    for j = 0:(n-1)
        denom *= (dim + j)
    end
    return num / denom
end

"""
    integrate_indices(U_idxs, U_bar_idxs, dim)

Low-level integration function using Weingarten calculus (Unitary).
"""
function integrate_indices(U_idxs::Vector{<:Tuple}, U_bar_idxs::Vector{<:Tuple}, dim)
    n = length(U_idxs)

    if n > 0
        u0 = U_idxs[1]
        ub0 = U_bar_idxs[1]
        if all(x -> x == u0, U_idxs) && all(x -> x == ub0, U_bar_idxs)
            return u0 == ub0 ? _unitary_single_entry_moment(n, dim) : 0
        end

        I = [x[1] for x in U_idxs]
        J = [x[2] for x in U_idxs]
        I_bar = [x[1] for x in U_bar_idxs]
        J_bar = [x[2] for x in U_bar_idxs]

        if all(x -> x == I[1], I) && all(x -> x == I_bar[1], I_bar)
            I[1] == I_bar[1] || return 0
            J_counts = _count_values(J)
            J_bar_counts = _count_values(J_bar)
            return J_counts == J_bar_counts ? _unitary_row_col_abs_moment(J_counts, dim) : 0
        end

        if all(x -> x == J[1], J) && all(x -> x == J_bar[1], J_bar)
            J[1] == J_bar[1] || return 0
            I_counts = _count_values(I)
            I_bar_counts = _count_values(I_bar)
            return I_counts == I_bar_counts ? _unitary_row_col_abs_moment(I_counts, dim) : 0
        end
    end

    I = [x[1] for x in U_idxs]
    J = [x[2] for x in U_idxs]
    I_bar = [x[1] for x in U_bar_idxs]
    J_bar = [x[2] for x in U_bar_idxs]

    valid = get_matching_permutations(I, I_bar)
    valid_taus = get_matching_permutations(J, J_bar)

    if isempty(valid) || isempty(valid_taus)
        return 0
    end

    cycle_counts = Dict{Vector{Int},Int}()
    n_v = length(valid)
    n_vt = length(valid_taus)
    if n_v * n_vt > 100 # Only if there are enough terms to justify overhead
        p = Progress(n_v * n_vt; dt = 10.0, desc = "Calculating cycle types... ")
        for sigma in valid
            for tau in valid_taus
                inv_tau = invperm(tau)
                P = [sigma[inv_tau[i]] for i = 1:n]
                cycle_type = get_cycle_type(P)
                cycle_counts[cycle_type] = get(cycle_counts, cycle_type, 0) + 1
                next!(p)
            end
        end
    else
        for sigma in valid
            for tau in valid_taus
                inv_tau = invperm(tau)
                P = [sigma[inv_tau[i]] for i = 1:n]
                cycle_type = get_cycle_type(P)
                cycle_counts[cycle_type] = get(cycle_counts, cycle_type, 0) + 1
            end
        end
    end

    if dim isa Integer
        total = zero(Rational{BigInt})
        for (ct, count) in cycle_counts
            total += count * weingarten(ct, dim)
        end
        return total
    else
        total_num = zero(Num)
        total_den = one(Num)

        for (ct, count) in cycle_counts
            wnum, wden = IntegrateUnitary.weingarten_raw(ct, dim)
            if isequal(total_den, 1)
                total_num = count * wnum
                total_den = wden
            else
                if isequal(total_den, wden)
                    total_num += count * wnum
                else
                    total_num = total_num * wden + count * wnum * total_den
                    total_den = total_den * wden
                end
            end
        end

        return total_num / total_den
    end
end

function get_matching_permutations(target::AbstractVector, source::AbstractVector)
    n = length(target)
    if n != length(source)
        return Vector{Vector{Int}}()
    end


    source_groups = Dict{Any,Vector{Int}}()
    for (idx, val) in enumerate(source)
        push!(get!(source_groups, val, Int[]), idx)
    end

    target_groups = Dict{Any,Vector{Int}}()
    for (idx, val) in enumerate(target)
        push!(get!(target_groups, val, Int[]), idx)
    end


    if length(source_groups) != length(target_groups)
        return Vector{Vector{Int}}()
    end
    for (val, t_idxs) in target_groups
        if !haskey(source_groups, val) || length(source_groups[val]) != length(t_idxs)
            return Vector{Vector{Int}}()
        end
    end


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
    n = length(p)
    if n > 64
        visited = falses(n)
        lengths = Int[]
        for i = 1:n
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
        sort!(lengths, rev = true)
        return lengths
    end

    visited = UInt64(0)
    lengths = Int[]
    sizehint!(lengths, n)
    for i = 1:n
        if (visited & (UInt64(1) << (i - 1))) == 0
            curr = i
            len = 0
            while (visited & (UInt64(1) << (curr - 1))) == 0
                visited |= (UInt64(1) << (curr - 1))
                curr = p[curr]
                len += 1
            end
            push!(lengths, len)
        end
    end
    sort!(lengths, rev = true)
    return lengths
end

"""
    integrate_indices_orthogonal(indices, dim)
    
Low-level integration function using Orthogonal Weingarten calculus.
Indices are a list of (i, j) for O_{ij}.
Formula: sum_{pi, sigma in PairPartitions} delta_pi(i) * delta_sigma(j) * Wg(pi, sigma)
"""

"""
    _orthogonal_row_sum_denom(dim, k)

Compute `prod_{j=0}^{k-1} (dim + 2j)`, the denominator for the orthogonal
row-sum shortcut when all row or column indices are uniform.
"""
function _orthogonal_row_sum_denom(dim, k)
    denom = one(Rational{BigInt})
    for j = 0:(k-1)
        denom *= (dim + 2*j)
    end
    return denom
end

function integrate_indices_orthogonal(indices::AbstractVector, dim)
    n = length(indices) # This is 2k
    if n % 2 != 0
        return 0
    end
    k = n ÷ 2

    I = [x[1] for x in indices]
    J = [x[2] for x in indices]

    is_I_uniform = all(x -> x == I[1], I)
    is_J_uniform = all(x -> x == J[1], J)

    if is_I_uniform || is_J_uniform
        denom = _orthogonal_row_sum_denom(dim, k)
        if is_I_uniform && is_J_uniform
            num = 1
            for i = 1:k
                num *= (2*i-1)
            end
            return num / denom
        elseif is_I_uniform
            valid_sigma = get_matching_pair_partitions_filtered(J)
            return length(valid_sigma) / denom
        else # is_J_uniform
            valid_pi = get_matching_pair_partitions_filtered(I)
            return length(valid_pi) / denom
        end
    end

    valid_pi = get_matching_pair_partitions_filtered(I)
    valid_sigma = get_matching_pair_partitions_filtered(J)

    if isempty(valid_pi) || isempty(valid_sigma)
        return 0
    end

    pi_counts = Dict{Any,Int}()
    for pi in valid_pi
        c_pi = canonicalize_pair_partition(pi)
        pi_counts[c_pi] = get(pi_counts, c_pi, 0) + 1
    end

    sigma_counts = Dict{Any,Int}()
    for sigma in valid_sigma
        c_sigma = canonicalize_pair_partition(sigma)
        sigma_counts[c_sigma] = get(sigma_counts, c_sigma, 0) + 1
    end

    type_counts = Dict{Vector{Int},BigInt}()
    for (c_pi, count_pi) in pi_counts
        for (c_sigma, count_sigma) in sigma_counts
            ct = get_full_cycle_type(c_pi, c_sigma)
            type_counts[ct] =
                get(type_counts, ct, zero(BigInt)) + BigInt(count_pi) * BigInt(count_sigma)
        end
    end

    if !(Symbolics.unwrap(dim) isa Number)
        w, type_to_idx, _ = get_weingarten_reduced_data(k, dim, return_rationals = true)

        # Exact rational summation to bypass Symbolics.simplify bugs
        local_total_rat = nothing
        for (ct, count) in type_counts
            term = _rational_mul(w[type_to_idx[ct]], count)
            if local_total_rat === nothing
                local_total_rat = term
            else
                local_total_rat = _rational_add(local_total_rat, term)
            end
        end

        if local_total_rat === nothing
            return Num(0)
        end
        return from_rational(local_total_rat)
    end

    w, type_to_idx, _ = get_weingarten_reduced_data(k, dim)
    T_val = eltype(w)
    total = zero(T_val)
    for (ct, count) in type_counts
        val = w[type_to_idx[ct]]
        total += count * val
    end
    return total
end

"""
    get_matching_pair_partitions_filtered(indices)

Returns all pair partitions of 1..2k such that for every pair (a,b), indices[a] == indices[b].
"""
function get_matching_pair_partitions_filtered(indices::Vector{Int})
    n = length(indices)


    counts = Dict{Int,Vector{Int}}()
    for (pos, val) in enumerate(indices)
        push!(get!(counts, val, Int[]), pos)
    end

    for (val, pos_list) in counts
        if length(pos_list) % 2 != 0
            return Vector{Vector{Tuple{Int,Int}}}()
        end
    end
    group_partitions = Vector{Vector{Vector{Tuple{Int,Int}}}}()

    for (val, pos_list) in counts
        k_local = length(pos_list)
        local_parts_idx = get_pair_partitions(k_local)

        real_parts = Vector{Vector{Tuple{Int,Int}}}()
        for p in local_parts_idx
            mapped = [(pos_list[a], pos_list[b]) for (a, b) in p]
            push!(real_parts, mapped)
        end
        push!(group_partitions, real_parts)
    end


    res = Vector{Vector{Tuple{Int,Int}}}()
    for combined in Iterators.product(group_partitions...)
        full_part = Vector{Tuple{Int,Int}}()
        for part in combined
            append!(full_part, part)
        end
        push!(res, full_part)
    end

    return res
end


"""
    integrate_indices_symplectic(indices, dim)

Low-level integration function using Symplectic Weingarten calculus.
Formula: sum_{pi, sigma} J_pi(i) * J_sigma(j) * Wg^Sp(pi, sigma)
"""
function integrate_indices_symplectic(indices::AbstractVector, dim)
    n = length(indices)
    if n % 2 != 0
        return 0
    end
    k = n ÷ 2
    partitions = get_pair_partitions(n)

    I = [x[1] for x in indices]
    J = [x[2] for x in indices]

    pi_contractions = Dict{Any,Any}()
    for pi in partitions
        val_I = compute_symplectic_contraction(pi, I, dim)
        if !_iszero(val_I)
            c_pi = canonicalize_pair_partition(pi)
            pi_contractions[c_pi] = get(pi_contractions, c_pi, 0) + val_I
        end
    end

    if isempty(pi_contractions)
        return 0
    end

    sigma_contractions = Dict{Any,Any}()
    for sigma in partitions
        val_J = compute_symplectic_contraction(sigma, J, dim)
        if !_iszero(val_J)
            c_sigma = canonicalize_pair_partition(sigma)
            sigma_contractions[c_sigma] = get(sigma_contractions, c_sigma, 0) + val_J
        end
    end

    if isempty(sigma_contractions)
        return 0
    end

    w, type_to_idx, _ =
        get_weingarten_reduced_data(k, -dim, return_rationals = !(dim isa Integer))

    if !(Symbolics.unwrap(dim) isa Number)
        # Use exact rational polynomial arithmetic to avoid Symbolics.simplify bugs
        local_total_rat = nothing
        for (c_pi, val_pi) in pi_contractions
            for (c_sigma, val_sigma) in sigma_contractions
                ct = get_full_cycle_type(c_pi, c_sigma)
                loops = length(ct)
                sign = ((-1)^loops)
                vp = Symbolics.unwrap(val_pi)
                vs = Symbolics.unwrap(val_sigma)
                vp_int = vp isa Number ? Int(vp) : Int(Symbolics.value(vp))
                vs_int = vs isa Number ? Int(vs) : Int(Symbolics.value(vs))
                int_coeff = BigInt(vp_int * vs_int * sign)
                term = _rational_mul(w[type_to_idx[ct]], int_coeff)
                if local_total_rat === nothing
                    local_total_rat = term
                else
                    local_total_rat = _rational_add(local_total_rat, term)
                end
            end
        end

        if local_total_rat === nothing
            return Num(0)
        end
        return from_rational(local_total_rat)
    end

    T = eltype(w)
    total = zero(T)

    n_pi = length(pi_contractions)
    n_sigma = length(sigma_contractions)
    if n_pi * n_sigma > 100
        p = Progress(
            n_pi * n_sigma;
            dt = 10.0,
            desc = "Calculating symplectic integrals... ",
        )
        for (c_pi, val_pi) in pi_contractions
            for (c_sigma, val_sigma) in sigma_contractions
                ct = get_full_cycle_type(c_pi, c_sigma)
                loops = length(ct)
                wg = ((-1)^loops) * w[type_to_idx[ct]]
                total += val_pi * val_sigma * wg
                next!(p)
            end
        end
    else
        for (c_pi, val_pi) in pi_contractions
            for (c_sigma, val_sigma) in sigma_contractions
                ct = get_full_cycle_type(c_pi, c_sigma)
                loops = length(ct)
                wg = ((-1)^loops) * w[type_to_idx[ct]]
                total += val_pi * val_sigma * wg
            end
        end
    end

    return total
end

"""
    integrate_indices_gue(indices, dim)
    
Low-level integration function using Wick's theorem for GUE.
Formula: sum_{pi in PairPartitions} prod_{(u, v) in pi} delta(i_u, j_v) * delta(j_u, i_v)
"""
function integrate_indices_gue(indices::AbstractVector, dim)
    n = length(indices) # Must be even
    partitions = get_pair_partitions(n)

    total = 0 // 1

    for pi in partitions
        term_val = 1
        possible = true

        for (u, v) in pi
            (i_u, j_u) = indices[u]
            (i_v, j_v) = indices[v]

            if !_symbolic_isequal(i_u, j_v) || !_symbolic_isequal(j_u, i_v)
                possible = false
                break
            end
        end

        if possible
            total += 1
        end
    end

    return total
end

function integrate_indices_gue(indices::Vector{Any}, dim)
    return integrate_indices_gue(Vector{Tuple{Any,Any}}(indices), dim)
end

"""
    integrate_indices_goe(indices, dim)
    
Low-level integration function using Wick's theorem for GOE.
Formula: sum_{pi in PairPartitions} prod_{(u, v) in pi} (delta(i_u, k_v)*delta(j_u, l_v) + delta(i_u, l_v)*delta(j_u, k_v))
where pair is H_{i_u j_u} and H_{k_v l_v} (indices re-labeled for clarity).
"""
function integrate_indices_goe(indices::AbstractVector, dim)
    n = length(indices) # Must be even
    partitions = get_pair_partitions(n)

    total = 0 // 1

    for pi in partitions
        term_val = 1
        possible = true

        for (u, v) in pi
            (i1, j1) = indices[u]
            (i2, j2) = indices[v]

            val_pair = 0 // 1

            match1 = _symbolic_isequal(i1, i2) && _symbolic_isequal(j1, j2)
            match2 = _symbolic_isequal(i1, j2) && _symbolic_isequal(j1, i2)

            if match1
                val_pair += 1
            end
            if match2
                val_pair += 1
            end

            if val_pair == 0
                possible = false
                break
            end
            term_val *= val_pair
        end

        if possible
            total += term_val
        end
    end

    return total
end


function _get_J(i, j, d)
    return symplectic_form(i, j, d)
end

function integrate_indices_gse(indices::AbstractVector, dim)
    n = length(indices)
    partitions = get_pair_partitions(n)
    total = 0 // 1

    for pi in partitions
        choice_combinations = collect(Iterators.product(fill([1, 2], n ÷ 2)...))
        for choices in choice_combinations
            term_val = 1 // 1
            possible = true
            n_type2 = 0
            for (p_idx, (u, v)) in enumerate(pi)
                (a, b) = indices[u]
                (c, d) = indices[v]

                choice = choices[p_idx]
                if choice == 1
                    if _symbolic_isequal(a, d) && _symbolic_isequal(b, c)
                        term_val *= 1
                    else
                        possible = false;
                        break
                    end
                else
                    n_type2 += 1
                    jac = _get_J(a, c, dim)
                    jbd = _get_J(b, d, dim)
                    if jac == 0 || jbd == 0
                        possible = false;
                        break
                    end
                    term_val *= (jac * jbd)
                end
            end

            if possible
                total += ((-1)^n_type2) * term_val
            end
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
    val = Num(1)

    for (u, v) in partition
        idx_u = indices[u]
        idx_v = indices[v]

        j_val = symplectic_form(idx_u, idx_v, dim)
        if _symbolic_isequal(j_val, 0)
            return 0
        end
        val *= j_val
    end
    return val
end

function symplectic_form(i, j, dim)
    m = dim / 2
    if _iszero(j - (i + m))
        return 1
    elseif _iszero(j - (i - m))
        return -1
    else
        return 0
    end
end

"""
    integrate_indices_ginue(u_indices, ub_indices, dim)

Low-level integration for Complex Ginibre Ensemble.
Formula: sum_{sigma in S_n} prod_m delta(i_m, ib_sigma(m)) * delta(j_m, jb_sigma(m))
"""
function integrate_indices_ginue(u_indices::AbstractVector, ub_indices::AbstractVector, dim)
    n = length(u_indices)
    if n != length(ub_indices)
        return 0
    end
    perms = permutations(1:n)
    total = 0
    for p in perms
        possible = true
        for m = 1:n
            (i, j) = u_indices[m]
            (ib, jb) = ub_indices[p[m]]
            if !_symbolic_isequal(i, ib) || !_symbolic_isequal(j, jb)
                possible = false
                break
            end
        end
        if possible
            total += 1
        end
    end
    return total
end

"""
    integrate_indices_ginoe(indices, dim)

Low-level integration for Real Ginibre Ensemble.
Formula: sum_{pi in PairPartitions} prod_{(u, v) in pi} delta(i_u, i_v) * delta(j_u, j_v)
"""
function integrate_indices_ginoe(indices::AbstractVector, dim)
    n = length(indices)
    if isodd(n)
        return 0
    end

    partitions = get_pair_partitions(n)
    total = 0
    for pi in partitions
        possible = true
        for (u, v) in pi
            (i1, j1) = indices[u]
            (i2, j2) = indices[v]
            if !_symbolic_isequal(i1, i2) || !_symbolic_isequal(j1, j2)
                possible = false
                break
            end
        end
        if possible
            total += 1
        end
    end
    return total
end

"""
    integrate_indices_ginse(indices, dim)

Low-level integration for Symplectic Ginibre Ensemble.
Uses duality with GinOE.
"""
function integrate_indices_ginse(indices::AbstractVector, dim)
    n = length(indices)
    if isodd(n)
        return 0
    end

    res_oe = integrate_indices_ginoe(indices, dim)

    final_sign = ((-1)^(n ÷ 2 + 1))
    return final_sign * res_oe
end

"""
    integrate_indices_diagonal(U_idxs, U_bar_idxs, dim)

Low-level integration function for the Diagonal Unitary group (Torus).
The integral is non-zero (equal to 1) only if the multisets of indices 
of U and U_bar are identical.
Indices are already checked to be diagonal (i == j) in process_term.
"""
function integrate_indices_diagonal(U_idxs::AbstractVector, U_bar_idxs::AbstractVector, dim)
    n = length(U_idxs)
    if n != length(U_bar_idxs)
        return 0
    end

    rows_u = [x[1] for x in U_idxs]
    rows_ubar = [x[1] for x in U_bar_idxs]
    sort!(rows_u)
    sort!(rows_ubar)

    if rows_u == rows_ubar
        return 1
    else
        return 0
    end
end
