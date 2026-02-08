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
        prod_val = 1 // 1
        for i = 1:length(part)
            for j = 1:part[i]
                hook_length = part[i] - i + conj_part[j] - j + 1
                term = d + j - i
                prod_val *= term // hook_length
            end
        end
        return prod_val
    else
        num = 1
        den = 1
        for i = 1:length(part)
            for j = 1:part[i]
                hook_length = part[i] - i + conj_part[j] - j + 1
                num *= (d + j - i)
                den *= hook_length
            end
        end
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
    # Wg(sigma, d) where sigma has cycle type `partition_type`.
    n = sum(partition_type)
    n_fact = factorial(big(n))

    # Iterate over all partitions of n
    parts = partitions(n)

    sum_val = 0 // 1

    for lam in parts
        # If length(lam) > d, s_lambda(1^d) = 0.
        if d isa Integer && length(lam) > d
            continue
        end

        # char_lam(1^n) = dimension f^lambda
        f_lam = character_at_id(lam)

        # char_lam(mu)
        chi_lam_mu = calculate_character(lam, partition_type)

        # s_lam(1^d)
        dim_lam = irrep_dimension(lam, d)

        term = (
            d isa Integer ? (big(f_lam)^2 * chi_lam_mu) // dim_lam :
            (big(f_lam)^2 * chi_lam_mu) / dim_lam
        )
        sum_val += term
    end

    return (d isa Integer ? sum_val // (n_fact^2) : sum_val * (1 // n_fact^2))
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

    # Recursive generation
    res = Vector{Vector{Tuple{Int,Int}}}()

    # helper
    function generate(current_pairs, remaining)
        if isempty(remaining)
            push!(res, current_pairs)
            return
        end

        first = remaining[1]
        # Try pairing `first` with each other element
        for i = 2:length(remaining)
            second = remaining[i]

            new_pairs = copy(current_pairs)
            push!(new_pairs, (first, second))

            new_remaining = copy(remaining)
            deleteat!(new_remaining, [1, i])

            generate(new_pairs, new_remaining)
        end
    end

    generate(Vector{Tuple{Int,Int}}(), collect(1:n))
    return res
end

"""
    count_loops(pi, sigma)

Count the number of loops `ℓ(π, σ)` in the multigraph formed by 
the union of two pair partitions `π` and `σ`.

The graph has `2k` vertices. Since both `π` and `σ` are perfect matchings, 
their union consists of disjoint cycles of even length.
"""
function count_loops(pi::Vector{Tuple{Int,Int}}, sigma::Vector{Tuple{Int,Int}})
    n = 2 * length(pi)
    # Every vertex has exactly two neighbors: one from pi, one from sigma.
    # Use a fixed-size vector for speed.
    neighs = Vector{Tuple{Int,Int}}(undef, n)
    for (u, v) in pi
        neighs[u] = (v, 0)
        neighs[v] = (u, 0)
    end
    for (u, v) in sigma
        # Store the sigma neighbor in the second slot
        u_p = neighs[u][1]
        neighs[u] = (u_p, v)
        v_p = neighs[v][1]
        neighs[v] = (v_p, u)
    end

    visited = falses(n)
    loops = 0
    for i = 1:n
        if !visited[i]
            loops += 1
            curr = i
            visited[curr] = true
            # Trace the cycle using the fact it's a collection of cycles.
            # Alternate between pi and sigma neighbors.
            while true
                # Neighbor in pi
                nxt = neighs[curr][1]
                if visited[nxt]
                    # Check neighbor in sigma
                    nxt = neighs[curr][2]
                    if visited[nxt]
                        break
                    end
                end
                curr = nxt
                visited[curr] = true
            end
        end
    end
    return loops
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
    get_weingarten_orthogonal_data(k, d)

Internal function to generate the **Orthogonal Weingarten matrix**. 
The matrix \$G\$ is a Gram matrix of size \$(2k-1)!! \\times (2k-1)!!\$ where 
\$G_{\\pi, \\sigma} = d^{\\ell(\\pi, \\sigma)}\$.
The Weingarten matrix is the inverse of \$G\$.
"""
@memoize function get_weingarten_orthogonal_data(k::Int, d)
    parts = get_pair_partitions(2*k)
    N = length(parts)

    # Pre-canonicalize all partitions for lookup
    canonical_parts = Any[canonicalize_pair_partition(p) for p in parts]
    lookup = Dict{Any,Int}(c => i for (i, c) in enumerate(canonical_parts))

    # Build Gram matrix
    # Determine type
    val_sample = d^1
    T = typeof(val_sample)
    if d isa Integer
        T = Rational{Int}
    end

    G = zeros(T, N, N)
    for i = 1:N
        for j = 1:N
            loops = count_loops(parts[i], parts[j])
            if d isa Integer
                G[i, j] = (d^loops) // 1
            else
                G[i, j] = d^loops
            end
        end
    end

    # Invert G
    Wg_mat = try
        inv(G)
    catch e
        error("Failed to invert O(d) Gram matrix for k=\$k. Error: \$e")
    end

    return Wg_mat, lookup
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
        d
    )
end

"""
    weingarten_orthogonal_val_canonical(c_pi, c_sigma, d)

Internal version of `weingarten_orthogonal_val` that assumes arguments are already canonical.
"""
function weingarten_orthogonal_val_canonical(c_pi, c_sigma, d)
    k = length(c_pi)
    Wg_mat, lookup = get_weingarten_orthogonal_data(k, d)

    idx_pi = get(lookup, c_pi, nothing)
    idx_sigma = get(lookup, c_sigma, nothing)

    if idx_pi === nothing || idx_sigma === nothing
        error("Partition not found in generated set for k=\$k")
    end

    return Wg_mat[idx_pi, idx_sigma]
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
    # Wg^Sp(d)(pi, sigma) = (-1)^count_loops(pi, sigma) * Wg^O(-d)(pi, sigma)
    # Note: d -> -d substitution.
    val_ortho = weingarten_orthogonal_val(pi, sigma, -d)
    return ((-1)^count_loops(pi, sigma)) * val_ortho
end
