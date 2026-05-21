

"""
    partial_trace(M, dims, subsystem)

Compute the partial trace of matrix `M` over the specified `subsystem`.
`dims` is a tuple/vector of dimensions for each subsystem.
`subsystem` is the index of the subsystem to be TRACED OUT.

Example: For a bipartite system with dims=(2, 2), partial_trace(M, (2, 2), 2)
returns the reduced density matrix of the first subsystem.
"""
function partial_trace(M, dims, subsystem)
    n = length(dims)
    if subsystem < 1 || subsystem > n
        throw(ArgumentError("subsystem index $subsystem out of range 1:$n"))
    end
    total = prod(dims)
    if size(M) != (total, total)
        throw(ArgumentError(
            "Matrix size $(size(M)) does not match product of subsystem dimensions ($total × $total)."
        ))
    end

    target_subs = filter(i -> i != subsystem, 1:n)
    target_dims = dims[target_subs]
    new_dim = prod(target_dims)

    res = similar(M, Union{eltype(M),Num}, new_dim, new_dim)

    function get_strides(d)
        strds = Vector{Int}(undef, length(d))
        s = 1
        for i = length(d):-1:1
            strds[i] = s
            s *= d[i]
        end
        return strds
    end

    full_strides = get_strides(dims)
    target_strides = get_strides(target_dims)

    function to_multi(idx, d, strds)
        m = Vector{Int}(undef, length(d))
        idx -= 1
        for i = 1:length(d)
            m[i] = div(idx, strds[i]) + 1
            idx %= strds[i]
        end
        return m
    end

    function to_flat(m, strds)
        idx = 0
        for i = 1:length(m)
            idx += (m[i] - 1) * strds[i]
        end
        return idx + 1
    end

    traced_dim = dims[subsystem]
    E = eltype(M)
    T = (E <: Number && !(E <: Symbolics.Num)) ? E : Any
    res = Matrix{T}(undef, new_dim, new_dim)

    for i = 1:new_dim
        for j = 1:new_dim
            m_i = to_multi(i, target_dims, target_strides)
            m_j = to_multi(j, target_dims, target_strides)

            val = (E <: Number) ? zero(E) : 0
            for k = 1:traced_dim
                full_m_i = Vector{Int}(undef, n)
                full_m_j = Vector{Int}(undef, n)

                curr_target = 1
                for s = 1:n
                    if s == subsystem
                        full_m_i[s] = k
                        full_m_j[s] = k
                    else
                        full_m_i[s] = m_i[curr_target]
                        full_m_j[s] = m_j[curr_target]
                        curr_target += 1
                    end
                end

                idx_i = to_flat(full_m_i, full_strides)
                idx_j = to_flat(full_m_j, full_strides)
                val += M[idx_i, idx_j]
            end
            res[i, j] = val
        end
    end

    if any(x -> x isa Symbolics.Num || SymbolicUtils.iscall(Symbolics.unwrap(x)), res)
        return map(Symbolics.wrap, res)
    end
    return res
end
