function integrate_indices_centered_permutation(
    all_indices::AbstractVector,
    dim,
)
    # Expand prod(P_u - 1/d)
    # Sum over subsets S of indices:
    # (-1/d)^(k - |S|) * E[prod_{u in S} P_u]

    k = length(all_indices)
    if k >= 8 * sizeof(UInt)
        throw(ArgumentError("Centered permutation integration for k=$k indices is infeasible (2^$k subsets)"))
    end
    total = 0

    # Pre-allocate to avoid repeated allocations in the loop
    subset_indices = Vector{Tuple{Any,Any}}(undef, k)

    for i = 0:(UInt(2)^k-1)
        # count_ones is a fast bitwise operation in Julia
        num_P = count_ones(i)

        # Fill the pre-allocated array and use a view
        idx = 1
        for j = 1:k
            if (i >> (j-1)) & 1 == 1
                subset_indices[idx] = all_indices[j]
                idx += 1
            end
        end

        # subset_indices_view is a view, so no allocation here
        subset_indices_view = @view subset_indices[1:num_P]

        term_val = integrate_indices_permutation(subset_indices_view, dim)
        factor = (-1/dim)^(k - num_P)

        total += term_val * factor
    end
    return total
end
