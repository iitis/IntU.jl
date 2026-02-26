# Weingarten calculus functions

"""
    conjugate_partition(part::Vector{Int})

Returns the conjugate partition \\lambda' of a partition \\lambda.
The conjugate partition is obtained by transposing the Young diagram of \\lambda.
Mathematically, \\lambda'_i = \\text{card}\\{j : \\lambda_j \\ge i\\}.
"""
function conjugate_partition(part::Vector{Int})
    isempty(part) && return Int[]
    # Number of columns is the first element (largest part)
    # The length of the i-th column is the number of parts >= i
    cols = part[1]
    [count(>=(i), part) for i = 1:cols]
end

"""
    character_at_id(part::Vector{Int})

Calculates the character of the symmetric group S_n at the identity element for the 
irreducible representation (irrep) corresponding to the partition `part` (\\\\lambda).
This is equivalent to the dimension f^\\\\lambda of the irrep.

The dimension is given by the **Hook Length Formula**:
```math
f^\\\\lambda = \\\\frac{n!}{\\\\prod_{(i,j) \\\\in \\\\lambda} h_{\\\\lambda}(i,j)}
```
where h_\\\\lambda(i,j) is the hook length of the cell (i,j) in the Young diagram of \\\\lambda.

Reference:
- Frame, J. S., Robinson, G. de B., & Thrall, R. M. (1954). The hook graphs of the symmetric group.
"""
function character_at_id(part::Vector{Int})
    n = sum(part)
    conj_part = conjugate_partition(part)

    # Hook length formula
    denom = 1
    for i = 1:length(part)
        for j = 1:part[i]
            term = part[i] - j + conj_part[j] - i + 1
            denom *= term
        end
    end

    return factorial(n) // denom
end

"""
    murnaghan_nakayama(lambda, mu)

Computes the character \\chi^\\lambda(\\mu) of the symmetric group S_n for the irrep 
\\lambda at the conjugacy class with cycle type \\mu using the **Murnaghan-Nakayama rule**.

The rule states:
```math
\\chi^\\lambda(\\mu) = \\sum_{T \\in RIM(\\lambda, \\mu)} (-1)^{\\text{ht}(T)}
```
where the sum is over "rim hook" tableaux of shape \\lambda and content \\mu.

Reference:
- Murnaghan, F. D. (1937). The characters of the symmetric group.
- Nakayama, T. (1940). On some finite group of substitutions.
"""
function murnaghan_nakayama(lambda::Vector{Int}, mu::Vector{Int})
    isempty(mu) && return isempty(lambda) ? 1 : 0
    isempty(lambda) && return 0

    # Check if partitions have same weight
    n = sum(lambda)
    if n != sum(mu)
        return 0
    end

    return calculate_character(lambda, mu)
end


"""
    irrep_dimension(part::Vector{Int}, d)

Computes the dimension s_\\lambda(1^d) of the irreducible representation of the unitary 
group U(d) (or the Schur polynomial at ones) corresponding to the partition \\lambda.

The dimension is given by the **Hook-Content Formula**:
```math
\\text{dim}_d(\\lambda) = \\prod_{(i,j) \\in \\lambda} \\frac{d + j - i}{h_\\lambda(i,j)}
```
where h_\\lambda(i,j) is the hook length and j-i is the content of the cell (i,j).

This implementation supports symbolic dimension d, returning a rational function in d.

Reference:
- Stanley, R. P. (1999). *Enumerative Combinatorics*, Vol. 2.
"""
@memoize function irrep_dimension(part::Vector{Int}, d)
    conj_part = conjugate_partition(part)

    if d isa Integer
        prod_val = one(Rational{BigInt})
        for i = 1:length(part)
            for j = 1:part[i]
                hook_length = part[i] - i + conj_part[j] - j + 1
                term = d + j - i
                prod_val *= Rational{BigInt}(term, hook_length)
            end
        end
        return prod_val
    else
        num = one(Num)
        den = one(BigInt)
        for i = 1:length(part)
            for j = 1:part[i]
                hook_length = part[i] - i + conj_part[j] - j + 1
                num *= (d + j - i)
                den *= hook_length
            end
        end
        # Keep rational coefficient separate to avoid float pollution
        return (1 // den) * num
    end
end

function get_binary_partition(part::Vector{Int})
    rev_part = reverse(part)
    prepended = [0; rev_part]
    diffs = diff(prepended)

    res = Int[]
    for d in diffs
        for _ = 1:d
            push!(res, 1)
        end
        push!(res, 0)
    end
    return res
end

# MN Algorithm
function mn_inner(R::Vector{Int}, m::Vector{Int}, t::Int)
    if t > length(m)
        return 1
    end

    target_len = m[t]

    c = 0
    s = 1
    limit = min(target_len, length(R))
    for j = 1:(limit-1)
        if R[j] == 0
            s = -s
        end
    end

    len_R = length(R)
    for i = 1:(len_R-target_len)
        if R[i] != R[i+target_len-1]
            s = -s
        end

        if R[i] == 1 && R[i+target_len] == 0
            R[i] = 0
            R[i+target_len] = 1
            c += s * mn_inner(R, m, t + 1)
            R[i] = 1
            R[i+target_len] = 0
        end
    end

    return c
end

@memoize function calculate_character(lambda::Vector{Int}, mu::Vector{Int})
    R = get_binary_partition(lambda)
    return mn_inner(R, mu, 1)
end


"""
    weingarten(partition_type::Vector{Int}, d)

Computes the **Unitary Weingarten function** \\text{Wg}(\\sigma, d) where \\sigma 
is a permutation with cycle type given by `partition_type`.

The Weingarten function is defined as the sum over irreducible representations of S_n:
```math
\\text{Wg}(\\sigma, d) = \\frac{1}{(n!)^2} \\sum_{\\lambda \\vdash n, \\ell(\\lambda) \\le d} \\frac{(f^\\lambda)^2 \\chi^\\lambda(\\sigma)}{s_\\lambda(1^d)}
```
where f^\\lambda is the dimension of the S_n irrep, \\chi^\\lambda(\\sigma) is the 
character, and s_\\lambda(1^d) is the dimension of the U(d) irrep.

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitay, orthogonal and symplectic groups. *Communications in Mathematical Physics*.
"""
@memoize function weingarten(partition_type::Vector{Int}, d)
    wnum, wden = weingarten_raw(partition_type, d)
    res = wnum / wden
    if !(d isa Integer)
        return Symbolics.simplify(res)
    end
    return res
end

# Internal version that returns (numerator, denominator) for symbolic d
@memoize function weingarten_raw(partition_type::Vector{Int}, d)
    n = sum(partition_type)
    n_fact = factorial(big(n))
    parts = partitions(n)

    if d isa Integer
        # ... (keep existing Integer logic if needed, but we mostly care about symbolic)
        sum_val = 0 // 1
        for lam in parts
            if length(lam) > d
                continue
            end
            f_lam = character_at_id(lam)
            chi_lam_mu = calculate_character(lam, partition_type)
            dim_lam = irrep_dimension(lam, d)
            coeff = (Rational{BigInt}(f_lam)^2) * chi_lam_mu
            sum_val += coeff / dim_lam
        end
        res = sum_val / (Rational{BigInt}(n_fact)^2)
        return res, one(Num)
    else
        contents_mult = Dict{Int,Int}()
        for lam in parts
            current_mult = Dict{Int,Int}()
            conj_lam = conjugate_partition(lam)
            for i = 1:length(lam)
                for j = 1:lam[i]
                    c = j - i
                    current_mult[c] = get(current_mult, c, 0) + 1
                end
            end
            for (c, m) in current_mult
                contents_mult[c] = max(get(contents_mult, c, 0), m)
            end
        end

        D = one(Num)
        for (c, m) in contents_mult
            D *= (d + c)^m
        end

        total_num = zero(Num)
        for lam in parts
            f_lam = character_at_id(lam)
            chi_lam_mu = calculate_character(lam, partition_type)
            den_lam = one(BigInt)
            poly_lam = one(Num)
            conj_lam = conjugate_partition(lam)
            current_mult = Dict{Int,Int}()
            for i = 1:length(lam)
                for j = 1:lam[i]
                    hook = lam[i] - i + conj_lam[j] - j + 1
                    den_lam *= hook
                    c = j - i
                    poly_lam *= (d + c)
                    current_mult[c] = get(current_mult, c, 0) + 1
                end
            end

            D_div_poly = one(Num)
            for (c, m) in contents_mult
                m_lam = get(current_mult, c, 0)
                diff_m = m - m_lam
                if diff_m > 0
                    D_div_poly *= (d + c)^diff_m
                end
            end

            coeff = (Rational{BigInt}(f_lam)^2) * chi_lam_mu
            # Pre-expand the terms to keep total_num as a flattened polynomial
            term_num = Symbolics.expand(coeff * den_lam * D_div_poly)
            total_num += term_num
        end

        common_den = Symbolics.expand(D * (Rational{BigInt}(n_fact)^2))
        return total_num, common_den
    end
end


# ==============================================================================
# Orthogonal and Symplectic Weingarten Calculus
# ==============================================================================

"""
    get_pair_partitions(n)

Generate all partitions of the set `{1, ..., n}` into `n/2` disjoint pairs. 
`n` must be even. The number of such partitions is given by the double factorial `(n-1)!!`.

These partitions are also known as **perfect matchings** of the complete graph `K_n`.
"""
@memoize function get_pair_partitions(n::Int)
    if n % 2 != 0
        return Vector{Vector{Tuple{Int,Int}}}()
    end
    if n == 0
        return [Vector{Tuple{Int,Int}}()]
    end

    k = n ÷ 2
    res = Vector{Vector{Tuple{Int,Int}}}()
    buf = Vector{Tuple{Int,Int}}(undef, k)
    used = falses(n)

    function generate(depth)
        if depth > k
            push!(res, copy(buf))
            return
        end

        # Find smallest unused index
        first = 0
        for i = 1:n
            if !used[i]
                first = i
                break
            end
        end
        used[first] = true

        # Pair `first` with each subsequent unused index
        for second = (first+1):n
            if !used[second]
                used[second] = true
                buf[depth] = (first, second)
                generate(depth + 1)
                used[second] = false
            end
        end
        used[first] = false
    end

    generate(1)
    return res
end

"""
    get_full_cycle_type(pi, sigma)

Returns the cycle type of the union of two pair partitions as a sorted partition of k.
The union forms cycles of lengths 2l_1, 2l_2, ... where sum l_i = k.
Returns [l_1, l_2, ...] sorted descending.
"""
function get_full_cycle_type(pi::Vector{Tuple{Int,Int}}, sigma::Vector{Tuple{Int,Int}})
    n = 2 * length(pi)
    neighs = Vector{Tuple{Int,Int}}(undef, n)
    for (u, v) in pi
        neighs[u] = (v, 0)
        neighs[v] = (u, 0)
    end
    for (u, v) in sigma
        u_p = neighs[u][1]
        neighs[u] = (u_p, v)
        v_p = neighs[v][1]
        neighs[v] = (v_p, u)
    end

    visited = falses(n)
    cycle_lengths = Int[]
    for i = 1:n
        if !visited[i]
            len = 0
            curr = i
            visited[curr] = true
            while true
                len += 1
                nxt = neighs[curr][1]
                if visited[nxt]
                    nxt = neighs[curr][2]
                    if visited[nxt]
                        break
                    end
                end
                curr = nxt
                visited[curr] = true
            end
            push!(cycle_lengths, len ÷ 2)
        end
    end
    sort!(cycle_lengths, rev = true)
    return cycle_lengths
end

function count_loops(pi::Vector{Tuple{Int,Int}}, sigma::Vector{Tuple{Int,Int}})
    return length(get_full_cycle_type(pi, sigma))
end


"""
    canonicalize_pair_partition(p)

Returns a sorted list of pairs (min, max), sorted by min, to uniquely identify a matching.
"""
function canonicalize_pair_partition(p::Vector{Tuple{Int,Int}})
    sorted_pairs = [(min(u, v), max(u, v)) for (u, v) in p]
    sort!(sorted_pairs, by = x -> x[1])
    return sorted_pairs
end

"""
    _safe_solve(M, rhs)

Solves `M \\ rhs` but falls back to the Moore-Penrose pseudoinverse `pinv` if `M` is singular.
Used to compute Weingarten function for small dimensions where the Gram matrix is singular 
(e.g., dSp(2) integration yielding d = -2).
"""
function _safe_solve(M::AbstractMatrix{T}, rhs) where {T}
    try
        return M \ rhs
    catch e
        if e isa LinearAlgebra.SingularException || (e isa ErrorException && occursin("singular", lowercase(e.msg)))
            # Fallback to pseudoinverse. We convert to Float64 as pinv requires SVD, 
            # which is not natively supported for Rational{BigInt} in Base.
            M_float = Float64.(M)
            rhs_float = Float64.(rhs)
            w_float = pinv(M_float) * rhs_float
            
            # Recover exact rational values (Weingarten equations always yield rationals)
            if T <: Rational
                IntType = T.parameters[1]
                return rationalize.(IntType, w_float, tol=1e-10)
            else
                return rationalize.(w_float, tol=1e-10)
            end
        end
        rethrow(e)
    end
end

"""
    get_weingarten_orthogonal_data(k, d)

Internal function to generate the **Orthogonal Weingarten matrix**. 
The matrix \$G\$ is a Gram matrix of size \$(2k-1)!! \\times (2k-1)!!\$ where 
\$G_{\\pi, \\sigma} = d^{\\ell(\\pi, \\sigma)}\$.
The Weingarten matrix is the inverse of \$G\$.
"""
@memoize function get_weingarten_reduced_data(k::Int, d)
    parts = get_pair_partitions(2*k)
    N = length(parts)
    pi_id = parts[1]

    # Group partitions by cycle type with respect to pi_id
    type_to_parts = Dict{Vector{Int},Vector{Int}}()
    for i = 1:N
        ct = get_full_cycle_type(parts[i], pi_id)
        push!(get!(type_to_parts, ct, Int[]), i)
    end

    cts = collect(keys(type_to_parts))
    sort!(cts, rev = true)
    n_types = length(cts)
    type_to_idx = Dict(ct => i for (i, ct) in enumerate(cts))

    val_sample = d^1
    T = typeof(val_sample)
    if d isa Integer
        T = Rational{BigInt}
    end

    M = zeros(T, n_types, n_types)
    for i = 1:n_types
        pi_lambda = parts[type_to_parts[cts[i]][1]]
        for j = 1:n_types
            sum_g = zero(T)
            for sigma_idx in type_to_parts[cts[j]]
                loops = count_loops(pi_lambda, parts[sigma_idx])
                if d isa Integer
                    sum_g += (Rational{BigInt}(d)^loops)
                else
                    sum_g += d^loops
                end
            end
            M[i, j] = sum_g
        end
    end

    id_type = fill(1, k)
    id_idx = type_to_idx[id_type]
    rhs = zeros(T, n_types)
    rhs[id_idx] = one(T)

    w = _safe_solve(M, rhs)

    return w, type_to_idx, type_to_parts
end

@memoize function get_weingarten_orthogonal_data(k::Int, d)
    parts = get_pair_partitions(2*k)
    N = length(parts)
    canonical_parts = Any[canonicalize_pair_partition(p) for p in parts]
    lookup = Dict{Any,Int}(c => i for (i, c) in enumerate(canonical_parts))

    w, type_to_idx, _ = get_weingarten_reduced_data(k, d)
    T = eltype(w)

    Wg_mat = zeros(T, N, N)
    for i = 1:N
        for j = 1:N
            ct = get_full_cycle_type(parts[i], parts[j])
            Wg_mat[i, j] = w[type_to_idx[ct]]
        end
    end

    return Wg_mat, lookup
end

function _safe_inv(G::AbstractMatrix)
    n = size(G, 1)
    if n == 1
        return map(x -> 1/x, G)
    end
    # For small symbolic matrices OR BigInt matrices, manual inverse might be safer/faster
    # Generic inv() for Rational{BigInt} can be slow or problematic if not optimized.
    # But usually Base.inv works fine for Rational{BigInt}.

    # Check if elements are Numbers but not AbstractFloat (to avoid precision loss)
    if eltype(G) <: Number && !(eltype(G) <: AbstractFloat)
        # For 2x2, simpler
        if n == 2
            a, b = G[1, 1], G[1, 2]
            c, d = G[2, 1], G[2, 2]
            det = a*d - b*c
            return (1//det) * [d -b; -c a]
        end
    end

    if eltype(G) <: Symbolics.Num
        if n == 2
            a, b = G[1, 1], G[1, 2]
            c, d = G[2, 1], G[2, 2]
            det = a*d - b*c
            return (1/det) * [d -b; -c a]
        end
    end
    return inv(G)
end

"""
    weingarten_orthogonal_val(pi, sigma, d)

Returns the **Orthogonal Weingarten function** value \\\\text{Wg}^O(\\\\pi, \\\\sigma, d).

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
function weingarten_orthogonal_val(
    pi::Vector{Tuple{Int,Int}},
    sigma::Vector{Tuple{Int,Int}},
    d,
)
    return weingarten_orthogonal_val_canonical(
        canonicalize_pair_partition(pi),
        canonicalize_pair_partition(sigma),
        d,
    )
end

"""
    weingarten_orthogonal_val_canonical(c_pi, c_sigma, d)

Internal version of `weingarten_orthogonal_val` that assumes arguments are already canonical.
"""
function weingarten_orthogonal_val_canonical(c_pi, c_sigma, d)
    k = length(c_pi)
    w, type_to_idx, _ = get_weingarten_reduced_data(k, d)

    ct = get_full_cycle_type(c_pi, c_sigma)
    return w[type_to_idx[ct]]
end



"""
    weingarten_symplectic_val(pi, sigma, d)

Returns the **Symplectic Weingarten function** value \\\\text{Wg}^{Sp}(\\\\pi, \\\\sigma, d).
Uses the duality relation:
```math
\\\\text{Wg}^{Sp}(\\\\pi, \\\\sigma, d) = (-1)^k \\\\text{Wg}^{O}(\\\\pi, \\\\sigma, -d)
```
where k is the number of pairs.

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
@memoize function weingarten_symplectic_val(pi, sigma, d)
    k = length(pi)
    # Wg^Sp(d) = (-1)^(k + loops) * Wg^O(-d)
    # We can get Wg^O(-d) from the reduced system directly
    w, type_to_idx, _ = get_weingarten_reduced_data(k, -d)

    ct = get_full_cycle_type(pi, sigma)
    loops = length(ct)

    val_ortho = w[type_to_idx[ct]]
    return ((-1)^loops) * val_ortho
end
