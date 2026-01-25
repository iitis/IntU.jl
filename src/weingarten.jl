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
