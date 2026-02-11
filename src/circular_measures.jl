"""
    dCUE(U, dim)

Defines the Circular Unitary Ensemble measure for U(d).
This is mathematically equivalent to the Haar measure on U(d).
"""
dCUE(U, dim) = dU(U, dim)

struct COEMeasure{M,D}
    S::M
    dim::D
end

struct CSEMeasure{M,D}
    S::M
    dim::D
end

@doc raw"""
    dCOE(S, dim)

Defines the Circular Orthogonal Ensemble (COE) measure on U(N).
Matrices in COE are symmetric unitary matrices.
They can be modeled as
```math
S = U U^T
```
where \$U \sim \text{Haar}(U(N))\$.
"""
dCOE(S, dim) = COEMeasure(S, dim)

@doc raw"""
    dCSE(S, dim)

Defines the Circular Symplectic Ensemble (CSE) measure on U(2N).
Matrices in CSE are self-dual unitary matrices:
```math
S = S^R = J S^T J^T
```
They can be modeled as
```math
S = U U^R
```
where \$U \sim \text{Haar}(U(2N))\$.
Note: The dimension `dim` corresponds to the size of the matrix, so it must be 2N.
"""
dCSE(S, dim) = CSEMeasure(S, dim)


"""
    integrate(expr, measure::COEMeasure)
"""
function integrate(expr::AbstractArray, measure::COEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::CSEMeasure)
    return map(e -> integrate(e, measure), expr)
end

function IntU.measure_info(measure::COEMeasure)
    S_sym = measure.S
    dim = measure.dim

    subs_dict = Dict{Any,Any}()
    S_atomic_lookup = Dict{Any,Tuple}()

    if S_sym isa AbstractArray
        for i = 1:size(S_sym, 1)
            for j = 1:size(S_sym, 2)
                s_ij_num = _safe_Num(S_sym[i, j])
                s_ij_un = Symbolics.unwrap(s_ij_num)
                s_atomic = Symbolics.variable(:S_atomic, i, j)
                s_bar_atomic = Symbolics.variable(:S_bar_atomic, i, j)

                S_atomic_lookup[Symbolics.unwrap(s_atomic)] = (i, j)
                S_atomic_lookup[Symbolics.unwrap(s_bar_atomic)] = (i, j)

                subs_dict[s_ij_un] = s_atomic
                subs_dict[Symbolics.unwrap(conj(s_ij_num))] = s_bar_atomic
                subs_dict[Symbolics.unwrap(Base.conj(s_ij_num))] = s_bar_atomic
            end
        end
    end

    matcher = LookupMatcher(
        Dict(k => v for (k, v) in S_atomic_lookup if occursin("S_atomic", string(k))),
        Dict(k => v for (k, v) in S_atomic_lookup if occursin("S_bar_atomic", string(k)))
    )

    return (subs_dict, matcher, dim, :COE)
end

function IntU.measure_info(measure::CSEMeasure)
    S_sym = measure.S
    dim = measure.dim

    subs_dict = Dict{Any,Any}()
    S_atomic_lookup = Dict{Any,Tuple}()

    if S_sym isa AbstractArray
        for i = 1:size(S_sym, 1)
            for j = 1:size(S_sym, 2)
                s_ij_num = _safe_Num(S_sym[i, j])
                s_ij_un = Symbolics.unwrap(s_ij_num)

                s_atomic = Symbolics.variable(:S_atomic, i, j)
                s_bar_atomic = Symbolics.variable(:S_bar_atomic, i, j)

                S_atomic_lookup[Symbolics.unwrap(s_atomic)] = (i, j)
                S_atomic_lookup[Symbolics.unwrap(s_bar_atomic)] = (i, j)

                subs_dict[s_ij_un] = s_atomic
                subs_dict[Symbolics.unwrap(conj(s_ij_num))] = s_bar_atomic
                subs_dict[Symbolics.unwrap(Base.conj(s_ij_num))] = s_bar_atomic
            end
        end
    end

    matcher = LookupMatcher(
        Dict(k => v for (k, v) in S_atomic_lookup if occursin("S_atomic", string(k))),
        Dict(k => v for (k, v) in S_atomic_lookup if occursin("S_bar_atomic", string(k)))
    )

    phys_dim = S_sym isa AbstractArray ? size(S_sym, 1) : 0
    return (subs_dict, matcher, dim, (:CSE, phys_dim))
end



@doc raw"""
    integrate_indices_coe(all_indices, dim)

Integration of COE terms by reducing to Haar integration.

## Index Mapping

``S = U U^T``, so each ``S_{ij}`` expands as:
```math
S_{ij} = \sum_k U_{ik} U_{jk}, \qquad \bar{S}_{ij} = \sum_k \bar{U}_{ik} \bar{U}_{jk}
```

Each ``S_{i_k j_k}`` produces two U-type indices ``(i_k, a_k)`` and ``(j_k, a_k)``
with shared dummy column ``a_k``. Similarly each ``\bar{S}_{p_k q_k}`` produces
two ``\bar{U}``-type indices with shared dummy column ``b_k``.

## Weingarten Integration

With ``m`` S-terms and ``m`` ``\bar{S}``-terms, we have ``2m`` U and ``2m``
``\bar{U}`` indices. The Haar integral gives:
```math
\sum_{\sigma, \tau \in S_{2m}} \delta_{\text{rows}}(\sigma) \cdot d^{\#\text{loops}(\tau)} \cdot \text{Wg}(\sigma\tau^{-1}, d)
```

The dummy column summation ``\sum_{a,b} \delta_{\text{cols}}(\tau)`` yields
``d^{\#\text{loops}(\tau)}``, where the number of loops is the number
of connected components in a bipartite graph:
- Vertices: ``a_1 \ldots a_m`` and ``b_1 \ldots b_m``
- Structural edges: ``a_k`` pairs columns ``(2k{-}1, 2k)``; same for ``b_k``
- ``\tau``-edges: ``a_{\lceil r/2 \rceil} \leftrightarrow b_{\lceil \tau(r)/2 \rceil}``

The loop count is computed via union-find on these ``2m`` variable nodes.
"""
function integrate_indices_coe(
    indices::Vector{Tuple{Int,Int}},
    U_bar_indices::Vector{Tuple{Int,Int}},
    dim,
)

    n_s = length(indices)
    n_s_bar = length(U_bar_indices)

    if n_s != n_s_bar
        return 0
    end

    m = n_s 
    n = 2 * m 

    U_rows = Vector{Int}(undef, n)
    for k = 1:m
        i, j = indices[k]
        U_rows[2k-1] = i
        U_rows[2k] = j
    end

    U_bar_rows = Vector{Int}(undef, n)
    for k = 1:m
        p, q = U_bar_indices[k]
        U_bar_rows[2k-1] = p
        U_bar_rows[2k] = q
    end

    valid_sigmas = get_matching_permutations(U_rows, U_bar_rows)
    if isempty(valid_sigmas)
        return 0
    end

    n_fact = try
        factorial(n)
    catch
        BigInt(1)
    end
    
    is_full_group = length(valid_sigmas) == n_fact

    permutations_n = collect(permutations(1:n))

    total = 0 // 1

    if is_full_group
        sum_wg = 0 // 1
        for p_type in partitions(n)
            counts = Dict{Int, Int}()
            for x in p_type
                counts[x] = get(counts, x, 0) + 1
            end
            denom = 1
            for (k, c) in counts
                denom *= (k^c) * factorial(c)
            end
            p_count = factorial(n) // denom
            
            val = weingarten(p_type, dim)
            sum_wg += p_count * val
        end
        
        if _symbolic_isequal(sum_wg, 0)
            return 0
        end

        loop_counts = Dict{Int, Int}()
        for tau in permutations_n
            uf = IntDisjointSets(2 * m)
            for r = 1:n
                u = div(r - 1, 2) + 1
                v_raw = tau[r]
                v = div(v_raw - 1, 2) + 1
                union!(uf, u, m + v)
            end
            loops = num_groups(uf)
            loop_counts[loops] = get(loop_counts, loops, 0) + 1
        end
        
        sum_loops = 0 // 1
        for (loops, count) in loop_counts
            sum_loops += count * (dim isa Integer ? dim : dim)^loops
        end
        
        return sum_wg * sum_loops

    else
        wg_coeffs = Dict{Vector{Int}, Any}()
        
        for tau in permutations_n
            uf = IntDisjointSets(2 * m)
            for r = 1:n
                u = div(r - 1, 2) + 1
                v_raw = tau[r]
                v = div(v_raw - 1, 2) + 1
                union!(uf, u, m + v)
            end
            loops = num_groups(uf)
            weight = (dim isa Integer ? dim : dim)^loops
            
            if _symbolic_isequal(weight, 0)
                continue
            end

            inv_tau = invperm(tau)
            for sigma in valid_sigmas
                P = [sigma[inv_tau[i]] for i = 1:n]
                cycles = get_cycle_type(P)
                wg_coeffs[cycles] = get(wg_coeffs, cycles, 0) + weight
            end
        end

        total = 0 // 1
        for (cycles, coeff) in wg_coeffs
            total += coeff * weingarten(cycles, dim)
        end
        return total
    end
end


@doc raw"""
    integrate_indices_cse(indices, U_bar_indices, dim, phys_dim)

Integration of CSE terms by reducing to Haar integration with symplectic structure.

## Index Mapping

``S = U U^R`` where ``U^R = J U^T J^T`` is the symplectic reverse.
Using the symplectic form ``J = [0\; I; -I\; 0]`` with ``\dim = 2N``:
```math
(U^R)_{aj} = \text{sign}(a)\,\text{sign}(j)\, U_{\text{pair}(j),\,\text{pair}(a)}
```
where ``\text{pair}(k)`` maps ``k \leftrightarrow k \pm N`` and
``\text{sign}(k) = J_{k,\text{pair}(k)}``.

Each ``S_{ij}`` produces two U-type indices with coefficients ``\text{sign}(j)\,\text{sign}(a_k)``,
and each ``\bar{S}_{pq}`` produces two ``\bar{U}``-type indices similarly.

## Signed Loop Count

Like COE, dummy column summation yields a factor per connected component.
However, the symplectic signs introduce parity constraints:
- Each variable ``a_k`` or ``b_k`` contributes ``\text{sign}(\cdot)`` to the product.
- ``\text{sign}(\text{pair}(x)) = -\text{sign}(x)``, so traversing a ``\text{pair}``
  operation flips the sign.
- A component with even total sign-flips contributes ``2N = d``;
  odd total sign-flips gives ``\sum \text{sign}(x) = 0``, killing the term.

The parity is tracked via union-find on these ``2m`` variable nodes.
"""
function integrate_indices_cse(
    indices::Vector{Tuple{Int,Int}},
    U_bar_indices::Vector{Tuple{Int,Int}},
    dim,
    phys_dim,
)
    n_s = length(indices)
    n_s_bar = length(U_bar_indices)

    if n_s != n_s_bar
        return 0
    end

    m = n_s
    n = 2 * m
    half_dim = div(phys_dim, 2)

    U_rows = Vector{Int}(undef, n)
    fixed_sign_coeff = 1

    get_sign(idx) = (idx <= half_dim ? 1 : -1)
    get_pair(idx) = (idx <= half_dim ? idx + half_dim : idx - half_dim)

    for k = 1:m
        i, j = indices[k]
        U_rows[2k-1] = i
        U_rows[2k] = get_pair(j)

        fixed_sign_coeff *= get_sign(j)
    end

    U_bar_rows = Vector{Int}(undef, n)
    for k = 1:m
        p, q = U_bar_indices[k]
        U_bar_rows[2k-1] = p
        U_bar_rows[2k] = get_pair(q)

        fixed_sign_coeff *= get_sign(q)
    end

    valid_sigmas = get_matching_permutations(U_rows, U_bar_rows)

    if isempty(valid_sigmas)
        return 0
    end

    permutations_n = collect(permutations(1:n))
    total_val = 0 // 1

    wg_coeffs = Dict{Vector{Int}, Any}()

    for tau in permutations_n
        possible = true
        uf_parent = collect(1:(2*m))
        uf_parity = zeros(Int, 2*m)

        function find_local(idx)
            if uf_parent[idx] == idx
                return idx, 0
            end
            root, root_parity = find_local(uf_parent[idx])
            uf_parent[idx] = root
            uf_parity[idx] = (uf_parity[idx] + root_parity) % 2
            return root, uf_parity[idx]
        end

        function unite_local(i, j, p)
            root_i, par_i = find_local(i)
            root_j, par_j = find_local(j)
            if root_i != root_j
                uf_parent[root_i] = root_j
                uf_parity[root_i] = (p - par_i - par_j) % 2
                if uf_parity[root_i] < 0
                    uf_parity[root_i] += 2
                end
                return true
            else
                current_rel = (par_i - par_j) % 2
                if current_rel < 0
                    current_rel += 2
                end
                return current_rel == p
            end
        end

        for r = 1:n
            var_idx_L = div(r - 1, 2) + 1
            is_pair_L = ((r - 1) % 2 == 1)
            tr = tau[r]
            var_idx_R = m + div(tr - 1, 2) + 1
            is_pair_R = ((tr - 1) % 2 == 1)
            parity = (is_pair_L != is_pair_R ? 1 : 0)

            if !unite_local(var_idx_L, var_idx_R, Int(parity) % 2)
                possible = false
                break
            end
        end

        if !possible
            continue
        end

        term_weight = 1 // 1
        processed = zeros(Bool, 2*m)
        for i = 1:(2*m)
            if !processed[i]
                root, _ = find_local(i)
                members = Int[]
                for j = i:(2*m)
                    rj, _ = find_local(j)
                    if rj == root
                        push!(members, j)
                        processed[j] = true
                    end
                end
                
                K = length(members)
                dist_sum = sum(find_local(j)[2] for j in members)
                prefactor = (-1)^dist_sum
                
                if K % 2 == 0
                    term_weight *= prefactor * dim
                else
                    term_weight = 0
                    break
                end
            end
        end

        if _symbolic_isequal(term_weight, 0)
            continue
        end

        inv_tau = invperm(tau)
        for sigma in valid_sigmas
            P = [sigma[inv_tau[i]] for i = 1:n]
            cycles = get_cycle_type(P)
            wg_coeffs[cycles] = get(wg_coeffs, cycles, 0) + term_weight
        end
    end

    for (cycles, coeff) in wg_coeffs
        total_val += coeff * weingarten(cycles, dim)
    end
    return fixed_sign_coeff * total_val
end
