# src/QI.jl

"""
    purity(rho)

Calculate the purity of a density matrix `rho`, defined as Tr(rho^2).
"""
function purity(rho)
    return tr(rho * rho)
end

"""
    average_purity(rho, measure)

Calculate the Haar-average purity of a density matrix `rho` under the given `measure`.
"""
function average_purity(rho, measure)
    return integrate(purity(rho), measure)
end

"""
    fidelity(rho, sigma)

Calculate the fidelity between two density matrices `rho` and `sigma`.
For pure states, this is |<phi|psi>|^2. In this package, we use the 
Hilbert-Schmidt inner product Tr(rho * sigma) as a convenience.
"""
function fidelity(rho, sigma)
    return tr(rho * sigma)
end

"""
    average_fidelity(rho, sigma, measure)

Calculate the Haar-average fidelity between `rho` and `sigma` under the given `measure`.
"""
function average_fidelity(rho, sigma, measure)
    return integrate(fidelity(rho, sigma), measure)
end

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
    # Total dimension should match size(M, 1)
    # We'll implement this using Symbolics-friendly indexing.
    # We can represent M[i1, i2, ..., in; j1, j2, ..., jn]
    # And sum over ik == jk for k == subsystem.
    
    # For now, let's implement bipartite specifically as it's the most common case,
    # or a generic one if possible.
    
    target_subs = filter(i -> i != subsystem, 1:n)
    target_dims = dims[target_subs]
    new_dim = prod(target_dims)
    
    res = Matrix{Complex{Num}}(undef, new_dim, new_dim)
    
    # Helper to calculate strides correctly
    function get_strides(d)
        strds = Vector{Int}(undef, length(d))
        s = 1
        for i in length(d):-1:1
            strds[i] = s
            s *= d[i]
        end
        return strds
    end
    
    full_strides = get_strides(dims)
    target_strides = get_strides(target_dims)
    
    # Helper to convert flat index to multi-index
    function to_multi(idx, d, strds)
        m = Vector{Int}(undef, length(d))
        idx -= 1
        for i in 1:length(d)
            m[i] = div(idx, strds[i]) + 1
            idx %= strds[i]
        end
        return m
    end
    
    # Helper to convert multi-index back to flat
    function to_flat(m, strds)
        idx = 0
        for i in 1:length(m)
            idx += (m[i] - 1) * strds[i]
        end
        return idx + 1
    end

    traced_dim = dims[subsystem]
    
    for i in 1:new_dim
        for j in 1:new_dim
            m_i = to_multi(i, target_dims, target_strides)
            m_j = to_multi(j, target_dims, target_strides)
            
            # Sum over the subsystem index
            val = 0
            for k in 1:traced_dim
                full_m_i = Vector{Int}(undef, n)
                full_m_j = Vector{Int}(undef, n)
                
                # Reconstruct full multi-indices
                curr_target = 1
                for s in 1:n
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
    
    return res
end
