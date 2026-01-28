# Weingarten calculus functions

"""
    conjugate_partition(part::Vector{Int})

Returns the conjugate of a partition.
"""
function conjugate_partition(part::Vector{Int})
    isempty(part) && return Int[]
    # Number of columns is the first element (largest part)
    # The length of the i-th column is the number of parts >= i
    cols = part[1]
    [count(>=(i), part) for i in 1:cols]
end

"""
    character_at_id(part::Vector{Int})

Calculates the character of the symmetric group identity element (dimension of the irrep),
given by the hook length formula.
"""
function character_at_id(part::Vector{Int})
    n = sum(part)
    conj_part = conjugate_partition(part)
    
    # Hook length formula
    denom = 1
    for i in 1:length(part)
        for j in 1:part[i]
            term = part[i] - j + conj_part[j] - i + 1
            denom *= term
        end
    end
    
    return factorial(n) // denom
end

"""
    murnaghan_nakayama(lambda, mu)

Computes the character of the symmetric group for partition `lambda` at class `mu`
using the Murnaghan-Nakayama rule.
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

Computes the dimension of the irreducible representation of U(d) corresponding to partition `part`.
Uses the hook-content formula:
dim(lambda) = Product_{(i,j) in lambda} (d + j - i) / h_{i,j}

This formulation supports symbolic `d`.
"""
@memoize function irrep_dimension(part::Vector{Int}, d)
    conj_part = conjugate_partition(part)
    cols = length(part) > 0 ? part[1] : 0
    
    # We need to iterate over all boxes (i, j) in the Young diagram
    prod_val = 1 // 1
    
    for i in 1:length(part)
        for j in 1:part[i]
            # Hook length h_{i,j} = lambda[i] - i + lambda'[j] - j + 1
            hook_length = part[i] - i + conj_part[j] - j + 1
            
            # Content c_{i,j} = j - i
            # Term = d + c_{i,j} = d + j - i
            term = d + j - i
            
            # Update product
            prod_val *= (d isa Integer ? term // hook_length : term / hook_length)
        end
    end
    
    return prod_val
end

function get_binary_partition(part::Vector{Int})
    rev_part = reverse(part)
    prepended = [0; rev_part]
    diffs = diff(prepended)
    
    res = Int[]
    for d in diffs
        for _ in 1:d
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
    for j in 1:(limit-1) 
        if R[j] == 0
            s = -s
        end
    end
    
    len_R = length(R)
    for i in 1:(len_R - target_len)
        if R[i] != R[i + target_len - 1]
            s = -s
        end
        
        if R[i] == 1 && R[i + target_len] == 0
            R[i] = 0
            R[i + target_len] = 1
            c += s * mn_inner(R, m, t + 1)
            R[i] = 1
            R[i + target_len] = 0
        end
    end
    
    return c
end

@memoize function calculate_character(lambda::Vector{Int}, mu::Vector{Int})
    R = get_binary_partition(lambda)
    return mn_inner(R, mu, 1)
end


@memoize function weingarten(partition_type::Vector{Int}, d)
    # Wg(sigma, d) where sigma has cycle type `partition_type`.
    n = sum(partition_type)
    
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
        
        term = (d isa Integer ? ((f_lam)^2 * chi_lam_mu) // dim_lam : ((f_lam)^2 * chi_lam_mu) / dim_lam)
        sum_val += term
    end
    
    return (d isa Integer ? sum_val // (factorial(n)^2) : sum_val / (factorial(n)^2))
end


# ==============================================================================
# Orthogonal and Symplectic Weingarten Calculus
# ==============================================================================

"""
    get_pair_partitions(n)

Generate all partitions of the set `{1, ..., n}` into pairs. `n` must be even.
Returns a list of partitions. Each partition is a list of pairs (Tuples).
"""
@memoize function get_pair_partitions(n::Int)
    if n % 2 != 0
        return Vector{Vector{Tuple{Int, Int}}}()
    end
    if n == 0
        return [Vector{Tuple{Int, Int}}()]
    end
    
    # Recursive generation
    res = Vector{Vector{Tuple{Int, Int}}}()
    
    # helper
    function generate(current_pairs, remaining)
        if isempty(remaining)
            push!(res, current_pairs)
            return
        end
        
        first = remaining[1]
        # Try pairing `first` with each other element
        for i in 2:length(remaining)
            second = remaining[i]
            
            new_pairs = copy(current_pairs)
            push!(new_pairs, (first, second))
            
            new_remaining = copy(remaining)
            deleteat!(new_remaining, [1, i])
            
            generate(new_pairs, new_remaining)
        end
    end
    
    generate(Vector{Tuple{Int, Int}}(), collect(1:n))
    return res
end

"""
    count_loops(pi, sigma)

Count the number of loops in the graph formed by superimposing two pair partitions `pi` and `sigma`.
The graph has vertices 1..2k. Edges correspond to pairs in `pi` and `sigma`.
Since both are perfect matchings, the union forms a set of disjoint cycles.
"""
function count_loops(pi::Vector{Tuple{Int, Int}}, sigma::Vector{Tuple{Int, Int}})
    n = 2 * length(pi)
    # Every vertex has exactly two neighbors: one from pi, one from sigma.
    # Use a fixed-size vector for speed.
    neighs = Vector{Tuple{Int, Int}}(undef, n)
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
    for i in 1:n
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
function canonicalize_pair_partition(p::Vector{Tuple{Int, Int}})
    sorted_pairs = [Pair(min(u,v), max(u,v)) for (u,v) in p]
    sort!(sorted_pairs, by=x->x.first)
    return sorted_pairs
end

"""
    get_weingarten_orthogonal_data(k, d)

Internal function to generate and memoize the Weingarten matrix and a lookup Map.
"""
@memoize function get_weingarten_orthogonal_data(k::Int, d)
    parts = get_pair_partitions(2*k)
    N = length(parts)
    
    # Pre-canonicalize all partitions for lookup
    canonical_parts = Any[canonicalize_pair_partition(p) for p in parts]
    lookup = Dict{Any, Int}(c => i for (i, c) in enumerate(canonical_parts))
    
    # Build Gram matrix
    # Determine type
    val_sample = d^1
    T = typeof(val_sample)
    if d isa Integer
        T = Rational{Int}
    end
    
    G = zeros(T, N, N)
    for i in 1:N
        for j in 1:N
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
        error("Failed to invert O(d) Gram matrix for k=$k. Error: $e")
    end
    
    return Wg_mat, lookup
end

"""
    weingarten_orthogonal_val(pi, sigma, d)

Returns the value Wg(pi, sigma) for O(d) using O(1) lookup.
"""
function weingarten_orthogonal_val(pi::Vector{Tuple{Int, Int}}, sigma::Vector{Tuple{Int, Int}}, d)
    k = length(pi)
    Wg_mat, lookup = get_weingarten_orthogonal_data(k, d)
    
    idx_pi = get(lookup, canonicalize_pair_partition(pi), nothing)
    idx_sigma = get(lookup, canonicalize_pair_partition(sigma), nothing)
    
    if idx_pi === nothing || idx_sigma === nothing
        error("Partition not found in generated set for k=$k")
    end
    
    return Wg_mat[idx_pi, idx_sigma]
end



"""
    weingarten_symplectic_val(pi, sigma, d)

Returns the value Wg(pi, sigma) for Sp(d).
Uses the relation Wg^Sp(d)(pi, sigma) = (-1)^k * Wg^O(-d)(pi, sigma)
where k is the number of pairs (length of pi or sigma).
"""
@memoize function weingarten_symplectic_val(pi, sigma, d)
    k = length(pi)
    # Wg^Sp(d)(pi, sigma) = (-1)^k * Wg^O(-d)(pi, sigma)
    # Note: d -> -d substitution.
    val_ortho = weingarten_orthogonal_val(pi, sigma, -d)
    return ((-1)^k) * val_ortho
end

