module Weingarten

using Combinatorics

export weingarten, conjugate_partition

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
    # Product of hook lengths
    # h_{i,j} = mu_i - j + mu'_j - i + 1
    
    # Using the formula from the reference code:
    # Total[partition]! / Product[ partition[[i]] - j + conjPart[[j]] - i + 1 ]
    
    denom = 1
    for i in 1:length(part)
        for j in 1:part[i]
            # In Mathematica indices are 1-based. Julia too.
            # conj_part must be accessed safely? 
            # conj_part length is lambda[1]. j goes up to lambda[i].
            # Since lambda[i] <= lambda[1] (partitions are sorted descending), j is valid index for conj_part.
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
function irrep_dimension(part::Vector{Int}, d)
    # Hook lengths are needed.
    # We already have logic for hook lengths in character_at_id implicitly.
    
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


# Simple MN Implementation helpers


# Porting the code's MN algorithm
# BinaryPartition:
# partition -> binary
# "Differences of Prepend(Reverse(partition), 0)" -> Map Append 1s and 0.
# Example: part=[2,1]. Reverse=[1,2]. Prepend 0 -> [0,1,2]. Diff -> [1,1].
# Wait, Differences([0,1,2]) = [1,1].
# Map: 1 -> [1,0], 1 -> [1,0] -> [1,0,1,0]?
# This seems to construct the path. 1=Right, 0=Up? Or 1=Up, 0=Right?
# Usually partition lambda inside d x d box.
# Boundary path from (d,0) to (0,d). Or similar.

# Let's implement `binary_partition` as in the reference.
function binary_partition(part::Vector{Int})
    # Reference: Flatten[Map[Append[l,Flatten[{ConstantArray[1,#],0}]]&,Differences[Prepend[Reverse[partition],0]]]]
    # This looks like it appends variable 1s then a 0.
    
    # Let's assume standard Maya diagram:
    # Pad partition with 0s? The code doesn't seem to pad explicitly in BinaryPartition.
    
    # Let's try to interpret the algorithm directly:
    # R (sequence of 0s and 1s).
    # Check R[i] == 1 and R[i+k] == 0.
    # Swap -> R[i]=0, R[i+k]=1.
    # Recurse.
    # Height of strip = number of 0s in R[i+1...i+k-1]. Or number of 0s skipped?
    # Code: `For[j=1, j < Min[m[t], Length[R]], j++, If[R[j]==0, s=-s]]`
    # This loops j from 1 to k-1 (if k < Len). Counts 0s?
    # s starts at 1. Swaps sign for every 0?
    
    # I'll rewrite this function in Julia carefully.
    nothing
end

function get_binary_partition(part::Vector{Int})
    # Differences[Prepend[Reverse[part], 0]]
    rev_part = reverse(part)
    prepended = [0; rev_part]
    diffs = diff(prepended) # diffs[i] = part[n-i+1] - part[n-i] ...
    
    # Wait, diffs must be non-negative because partition is sorted descending?
    # Reverse is ascending. [1, 2]. Prepend 0 -> [0, 1, 2]. Diffs -> [1, 1]. Positive.
    
    # For each diff d:
    # Append d 1s, then a 0.
    # Wait, reference says: `Append[l, Flatten[{ConstantArray[1,#], 0}]]`
    # `l={}` (empty list passed to Map?). No, `l` is local var?
    # The Map seems to produce a list of lists.
    # Flatten at the end.
    
    res = Int[]
    for d in diffs
        for _ in 1:d
            push!(res, 1)
        end
        push!(res, 0)
    end
    return res
end

# MN Algorithm from Reference
function mn_inner(R::Vector{Int}, m::Vector{Int}, t::Int)
    if t > length(m)
        return 1
    end
    
    target_len = m[t] # cycle length
    
    c = 0
    # Search for all i such that swapping R[i] and R[i+target_len] is valid (1 -> 0)
    # This corresponds to removing a border strip of length target_len.
    # The sign is determined by the parity of the height of the border strip.
    # We maintain the parity s incrementally using a sliding window approach.
    
    s = 1
    limit = min(target_len, length(R))
    # Initialize s based on the first window
    for j in 1:(limit-1) 
        if R[j] == 0
            s = -s
        end
    end
    
    len_R = length(R)
    for i in 1:(len_R - target_len)
        # Update parity s: if the bit leaving the window (R[i]) is different 
        # significantly from the one technically entering or affecting parity, flip s.
        # Specifically, compare R[i] and R[i+target_len-1].
        if R[i] != R[i + target_len - 1]
            s = -s
        end
        
        # Check if we can validly remove a border strip (1 at start, 0 at end of length)
        if R[i] == 1 && R[i + target_len] == 0
            # Perform swap (remove strip)
            R[i] = 0
            R[i + target_len] = 1
            
            # Recurse for the next cycle component
            c += s * mn_inner(R, m, t + 1)
            
            # Backtrack
            R[i] = 1
            R[i + target_len] = 0
        end
    end
    
    return c
end

# But wait, the `s` update logic in reference:
# At i=1: `If[R[[1]] != R[[m]], s=-s]`.
# `s` was computed on `1..m-1`.
# If `R[1]` (which is `i`) != `R[m]` (which is `i+m-1`?), s flips.
# This seems to verify if the number of 0s changes parity?
# Yes.
# If we have a range, and we shift it by 1.
# The parity of 0s changes if the element leaving (left) and element entering (right) are different.
# Similar to sliding window sum.
# If `R[i]` (leaving) == `R[i+len]` (entering), count of 0s changes by 0 or (+1 and -1). Parity same.
# If they differ, count changes by +/- 1. Parity flips.
# Correct.
# The range of interest for height is `i+1` to `i+m-1`.
# The `s` computed initially is for `1` to `m-1`.
# When checking `i=1`, we want parity of `2` to `m`?
# Wait. Valid extraction at `i` involves `R[i]=1`, `R[i+m]=0`.
# The height is number of 0s in `R[i+1 ... i+m-1]`.
# My initial `s` was for `1 ... m-1`.
# At `i=1`: I want `2 ... m`? No.
# I want `i+1 ... i+m-1`.
# At `i=1`, that is `2 ... m`.
# The update `R[1] != R[m]` updates the range from `1..m-1` to `2..m`?
# `R[1]` leaves. `R[m]` enters.
# Matches my logic.

# BUT `m` in code is `m[[t]]` which is the FULL length of the strip.
# The inner range length is `m-1`.
# Initial `s` computes for indices `1` to `m-1`.
# Loop i starts at 1.
# Update `s`: checks `R[1]` vs `R[m]`.
# This updates `s` to be valid for `2` to `m`?
# But we verify swap at `1` and `1+m`.
# The "between" range is `2` to `m`.
# The standard height definition is "number of rows spanned - 1" or "number of leg cells".
# In Maya diagram, height is number of 0s strictly between the 1 and 0.
# So if we swap `R[i]` and `R[i+m]`, the height is count of 0s in `R[i+1...i+m-1]`.
# The initial `s` loop counts 0s in `1...m-1`.
# At `i=1`, `s` is updated by `R[1]` leaving and `R[m]` entering.
# This results in count of 0s in `2...m`.
# But we want `2...m`? (`i+1 ... i+m-1`).
# If `i=1`, `m=m`, range is `2...m`.
# Is `m` in `R[m]` the `m`-th element?
# `m[t]` is the length.
# If we swap `i` and `i+m`, the indices between are `i+1` to `i+m-1`.
# Example m=2. Swap `i` and `i+2`. Between is `i+1`.
# Code loop `j < Min[m, Len]`. If `m=2`, `j` runs for `1`.
# `s` checks `R[1]`.
# Update: `R[1]` vs `R[2]`.
# Result `s` corresponds to `R[2]`.
# If `i=1`, we want count in `R[2]`.
# So `s` calculation + update is correct.
# Wait, for `i=1`, the loop updates `s` immediately.
# So `s` (initial) is for `1..m-1`.
# Adjusted `s` is for `2..m`.
# But `i+m-1` IS `m`.
# So `s` becomes valid for `2..m`.
# The range between `i` and `i+m` (exclusive) is `i+1 .. i+m-1`.
# If `i=1`, range is `2 .. m`.
# So yes, it is correct.

# One detail: `Min[m[t], Length[R]]`.
# If `Length[R] < m[t]`, returns 1 (no swap possible).
# But the `i` loop constraint `i < Length` prevents execution anyway.
# `limit = min(m[t], length(R))` in initialization.

function calculate_character(lambda::Vector{Int}, mu::Vector{Int})
    R = get_binary_partition(lambda)
    # R needs to be large enough?
    # MurnaghanNakayama adds zeros? Reference "BinaryPartition" just flattens.
    # But effectively lambda is padded with 0?
    # "Prepend[..., 0]" ensures a starting 0.
    # The MN algorithm works on this specific representation.
    
    return mn_inner(R, mu, 1)
end

function weingarten(partition_type::Vector{Int}, d)
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

end # module
