function integrate_indices_centered_permutation(all_indices::AbstractVector, dim)

    k = length(all_indices)
    if k >= 8 * sizeof(UInt)
        throw(
            ArgumentError(
                "Centered permutation integration for k=$k indices is infeasible (2^$k subsets)",
            ),
        )
    end
    total = 0

    subset_indices = Vector{Tuple{Any,Any}}(undef, k)

    for i = 0:(UInt(2)^k-1)
        num_P = count_ones(i)

        idx = 1
        for j = 1:k
            if (i >> (j-1)) & 1 == 1
                subset_indices[idx] = all_indices[j]
                idx += 1
            end
        end

        subset_indices_view = @view subset_indices[1:num_P]

        term_val = integrate_indices_permutation(subset_indices_view, dim)
        factor = (-1/dim)^(k - num_P)

        total += term_val * factor
    end
    return total
end
