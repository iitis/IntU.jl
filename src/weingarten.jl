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
function get_pair_partitions(n::Int)
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
function count_loops(pi, sigma)
    n = 2 * length(pi)
    adj = Dict{Int, Vector{Int}}()
    for i in 1:n
        adj[i] = Int[]
    end
    
    for (u, v) in pi
        push!(adj[u], v)
        push!(adj[v], u)
    end
    for (u, v) in sigma
        push!(adj[u], v)
        push!(adj[v], u)
    end
    
    # Count connected components
    visited = falses(n)
    loops = 0
    for i in 1:n
        if !visited[i]
            loops += 1
            # BFS/DFS
            q = [i]
            visited[i] = true
            while !isempty(q)
                curr = pop!(q)
                for neighbor in adj[curr]
                    if !visited[neighbor]
                        visited[neighbor] = true
                        push!(q, neighbor)
                    end
                end
            end
        end
    end
    return loops
end

"""
    orthogonal_gram_matrix(k, d)

Compute the Gram matrix for the Orthogonal group O(d) for 2k indices.
Rows and columns are indexed by pair partitions of 2k elements.
G_{pi, sigma} = d^{loops(pi, sigma)}
"""
function orthogonal_gram_matrix(k::Int, d)
    partitions = get_pair_partitions(2*k)
    N = length(partitions)
    
    # Determine type
    val_sample = d^1
    T = typeof(val_sample)
    if d isa Integer
        T = Rational{Int}
    end
    
    G = zeros(T, N, N)
    
    for i in 1:N
        for j in 1:N
            loops = count_loops(partitions[i], partitions[j])
            if d isa Integer
                G[i, j] = (d^loops) // 1
            else
                G[i, j] = d^loops
            end
        end
    end
    return G, partitions
end


"""
    weingarten_orthogonal_matrix(k, d)

Computes the Weingarten matrix (inverse of Gram matrix) for O(d).
Returns (Wg_matrix, partitions).
"""
@memoize function weingarten_orthogonal_matrix(k::Int, d)
    G, parts = orthogonal_gram_matrix(k, d)
    # Invert G
    try
        Wg = inv(G)
        return Wg, parts
    catch e
        # If singular or symbolic issue
        error("Failed to invert O(d) Gram matrix for k=$k. Error: $e")
    end
end

"""
    weingarten_orthogonal_val(pi, sigma, d)

Returns the value Wg(pi, sigma) for O(d).
"""
@memoize function weingarten_orthogonal_val(pi, sigma, d)
    k = length(pi) # pi is list of pairs
    # Note: this is inefficient if called repeatedly for same k but different pi,sigma
    # But we memoize the matrix generation.
    Wg_mat, parts = weingarten_orthogonal_matrix(k, d)
    
    # Find index of pi and sigma
    # This search is O(N) where N=(2k-1)!!, might be slow for large k.
    # We can optimize by canonicalizing representation.
    # For now, simplistic find.
    
    # Canonicalize pairs for comparison: (min, max), sorted by min
    function canonicalize(p)
        sorted_pairs = [Pair(min(u,v), max(u,v)) for (u,v) in p]
        sort!(sorted_pairs, by=x->x.first)
        return sorted_pairs
    end
    
    pi_c = canonicalize(pi)
    sigma_c = canonicalize(sigma)
    
    parts_c = [canonicalize(p) for p in parts]
    
    idx_pi = findfirst(isequal(pi_c), parts_c)
    idx_sigma = findfirst(isequal(sigma_c), parts_c)
    
    if idx_pi === nothing || idx_sigma === nothing
        error("Partition not found in generated set")
    end
    
    return Wg_mat[idx_pi, idx_sigma]
end


"""
    weingarten_symplectic_val(pi, sigma, d)

Returns the value Wg(pi, sigma) for Sp(d).
Uses the relation Wg^Sp(d)(pi, sigma) = (-1)^(k + loops(pi, sigma)) * Wg^O(-d)(pi, sigma)? 
Actually, standard relation: Wg^Sp(d) is related to O(-d).
Using Collins 2006:
Wg_{Sp(d)}(pi, sigma) = (-1)^{l(pi) + l(sigma)} Wg_{O(-d)}(pi, sigma) ?
where l(pi) is number of crossings? No.

Let's use the property that for Sp(d), the Gram matrix is G_{ij} = (-d)^{loops} * (-1)^{sum of crossings is tricky}.
Actually, the most robust way via duality is:
Wg^{Sp(d)}(pi, sigma) = (-1)^{k + loops(pi, sigma)} times [Wg^{O(-d)}(pi, sigma)] NO.

Correct Duality (M. Novak, "Weingarten calculus...", 2012?):
Wg^{Sp(d)}(pi, sigma) = (-1)^{crossings(pi) + crossings(sigma)} Wg^{O(-d)}(pi, sigma) ???

Let's assume the user is happy with just O(d) working perfectly first.
For Sp(d), I will implement a placeholder that warns or tries the O(-d) trick if I can verify it.
Actually, let's implement the `d -> -d` substitution and see.
Most sources suggest replacing d with -d connects O and Sp.
Let's try: Wg^{Sp(d)} = (-1)^{pairing_sign} Wg^{O(-d)}.
The pairing sign depends on the rewriting of J.

For now, I'll export `weingarten_orthogonal_val`.
"""
function weingarten_symplectic_val(pi, sigma, d)
     # Placeholder: Requires handling signs properly.
     # If we assume we mapped pairs correctly with J, maybe it's just O(-d)?
     # TODO: Rigorous check.
     error("Symplectic Weingarten not fully verified yet. Use O(d) with caution.")
end

