# Circular Ensembles (COE, CUE, CSE)

"""
    dCUE(U, dim)

Defines the Circular Unitary Ensemble measure for U(d).
This is mathematically equivalent to the Haar measure on U(d).
"""
dCUE(U, dim) = dU(U, dim)

# Dummy types for new measures
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
where $U \sim \text{Haar}(U(N))$.
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
where $U \sim \text{Haar}(U(2N))$.
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

function fallback_integrate(expr, measure::COEMeasure)
    S_sym = measure.S
    dim = measure.dim

    subs_dict = Dict{Any,Any}()
    S_atomic_lookup = Dict{Any,Tuple{Int,Int}}()

    if S_sym isa AbstractArray
        for i = 1:size(S_sym, 1)
            for j = 1:size(S_sym, 2)
                s_ij_num = _safe_Num(S_sym[i, j])
                s_ij_un = Symbolics.unwrap(s_ij_num)
                # We use :S as a marker
                s_atomic = Symbolics.variable(:S_atomic, i, j)
                s_bar_atomic = Symbolics.variable(:S_bar_atomic, i, j)

                S_atomic_lookup[Symbolics.unwrap(s_atomic)] = (i, j)
                S_atomic_lookup[Symbolics.unwrap(s_bar_atomic)] = (i, j)

                subs_dict[s_ij_un] = s_atomic

                c_ij_un = Symbolics.unwrap(conj(s_ij_num))
                subs_dict[c_ij_un] = s_bar_atomic

                bc_ij_un = Symbolics.unwrap(Base.conj(s_ij_num))
                subs_dict[bc_ij_un] = s_bar_atomic
            end
        end
    end

    # Use standard LookupMatcher. We map everything to :S type internally for process_term
    # But wait, `process_term` uses types :U, :U_bar.
    # We can reuse the same matcher if we map S_atomic -> (:U, i, j) and S_bar_atomic -> (:U_bar, i, j).
    # Then `process_term` will return U_indices and U_bar_indices.
    # We can then pass these to `integrate_indices_coe`.

    # Remap lookup for matcher
    matcher = LookupMatcher(
        Dict(
            k => (v[1], v[2]) for
            (k, v) in S_atomic_lookup if occursin("S_atomic", string(k))
        ),
        Dict(
            k => (v[1], v[2]) for
            (k, v) in S_atomic_lookup if occursin("S_bar_atomic", string(k))
        ),
    )

    return _robust_real_num(_integrate_core(expr, dim, subs_dict, matcher, :COE))
end

function fallback_integrate(expr, measure::CSEMeasure)
    S_sym = measure.S
    dim = measure.dim

    subs_dict = Dict{Any,Any}()
    S_atomic_lookup = Dict{Any,Tuple{Int,Int}}()

    # Similar to COE
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

                c_ij_un = Symbolics.unwrap(conj(s_ij_num))
                subs_dict[c_ij_un] = s_bar_atomic

                bc_ij_un = Symbolics.unwrap(Base.conj(s_ij_num))
                subs_dict[bc_ij_un] = s_bar_atomic
            end
        end
    end

    # Use standard LookupMatcher.
    matcher = LookupMatcher(
        Dict(
            k => (v[1], v[2]) for
            (k, v) in S_atomic_lookup if occursin("S_atomic", string(k))
        ),
        Dict(
            k => (v[1], v[2]) for
            (k, v) in S_atomic_lookup if occursin("S_bar_atomic", string(k))
        ),
    )

    phys_dim = S_sym isa AbstractArray ? size(S_sym, 1) : 0
    return _robust_real_num(
        _integrate_core(expr, dim, subs_dict, matcher, (:CSE, phys_dim)),
    )
end



@doc raw"""
    integrate_indices_coe(all_indices, dim)

Integration of COE terms. 
Mapping:
```math
S_{ij} = \sum_k U_{ik} U_{jk}
```
Each S index pair (i,j) maps to two U index pairs:
(i, k\_new) and (j, k\_new) where k\_new is summed over.

Note: S is symmetric, so distinguishing U and U^T in $U U^T$ matters for index placement.
```math
S_{ij} = (U U^T)_{ij} = \sum_k U_{ik} (U^T)_{kj} = \sum_k U_{ik} U_{jk}
```
This is effectively $U_{ik}$ and $U_{jk}$. Both are "U" type (not conjugate).

```math
\bar{S}_{ij} = \overline{\sum_k U_{ik} U_{jk}} = \sum_k \bar{U}_{ik} \bar{U}_{jk}
```
This is two "U\_bar" type.

So detailed mapping:
- For every occurrence of $S_{ij}$, we generate two geometric U terms: $U_{ia}$ and $U_{ja}$ where 'a' is a fresh summation index.
- For every occurrence of $\bar{S}_{ij}$, we generate two geometric U\_bar terms: $\bar{U}_{ib}$ and $\bar{U}_{jb}$ where 'b' is a fresh summation index.

We then call `integrate_indices` (Haar) on the resulting collection of U and U\_bar indices.
Since we are summing over the dummy indices (a and b), we effectively perform the Weingarten summation.
However, since 'a' and 'b' are dummy indices, summing over them is equivalent to contracting the tensors.
The `integrate_indices` function expects fixed indices.
But wait, we can leave 'a', 'b' as symbolic, get the result (which will contain delta functions involving a, b), and then perform the sum.
Alternatively, we can modify the Weingarten summation to handle the contraction directly.

Let's use the property of Haar integral:
```math
\int U_{i_1 j_1} ... dU = \sum_{\sigma, \tau} \delta_{i, k_\sigma} \delta_{j, l_\tau} Wg
```

Here, our indices are partly fixed (external i, j) and partly dummy (internal summation).
Let's call the dummy indices $d_1, d_2, ...$.
We have $S_{i_1 j_1} ... S_{i_m j_m} \bar{S}_{p_1 q_1} ... \bar{S}_{p_m q_m}$.
Maps to:
$U_{i_1 a_1} U_{j_1 a_1} ... U_{i_m a_m} U_{j_m a_m}$
$\bar{U}_{p_1 b_1} \bar{U}_{q_1 b_1} ... \bar{U}_{p_m b_m} \bar{U}_{q_m b_m}$

Total 2m U terms and 2m U_bar terms.
Let's collect indices.
U indices (row, col):
1. (i_1, a_1)
2. (j_1, a_1)
...
(2m-1). (i_m, a_m)
(2m). (j_m, a_m)

U_bar indices (row, col):
1. (p_1, b_1)
...
(2m). (q_m, b_m)

Integration result is:
```math
\sum_{\sigma, \tau \in S_{2m}} \delta_{rows} \delta_{cols} Wg(\sigma \tau^{-1}, dim)
```

$\delta_{rows}$ connects external indices $i, j, p, q$.
$\delta_{cols}$ connects the dummy indices a, b.
```math
\delta_{cols} = \prod_{r=1}^{2m} \delta(col(U)_r, col(Ubar)_{\tau(r)})
```
Col(U)_r is $a_{\lceil r/2 \rceil}$.
Col(Ubar)_s is $b_{\lceil s/2 \rceil}$.
So we have $\delta(a_{\lceil r/2 \rceil}, b_{\lceil \tau(r)/2 \rceil})$.
We sum over all $a_1...a_m$ and $b_1...b_m$.
Wait, the summation over dummy indices is independent for each 'a' and 'b' if they weren't connected by deltas.
But here we have restrictions.
Actually, we sum over $a_k$ and $b_k$.
The expression is $\sum_{a, b} \delta_{cols} ...$
The $\delta_{cols}$ forces certain a's and b's to be equal.
Specifically, for a fixed $\tau$, the product of deltas identifying a's and b's can be viewed as a graph or partition.
Each connected component in this identification must share the same value.
Since we sum over $a_k, b_k \in 1..dim$, each free connected component contributes a factor of $dim$.
So the sum over a, b of $\delta_{cols}$ is simply $dim^{Loops(\tau)}$.

Let's trace the loops.
We have m "source" nodes $a_1...a_m$ (each appearing twice) and m "target" nodes $b_1...b_m$ (each appearing twice).
Or simpler:
We have 2m positions for U-cols.
Positions 1,2 have value $a_1$.
Positions 3,4 have value $a_2$.
...
Positions 2m-1, 2m have value $a_m$.

We have 2m positions for Ubar-cols.
Positions 1,2 have value $b_1$.
...
Positions 2m-1, 2m have value $b_m$.

$\tau$ maps U-col position $r$ to Ubar-col position $\tau(r)$.
This imposes constraints.
If we view identifying $a_k$ with $b_l$ as an edge, we are counting connected components.
Vertices: $a_1...a_m, b_1...b_m$. Total 2m vertices.
Constraints:
From structure: 'U-col 2k-1' is same node as 'U-col 2k' (both are $a_k$).
From structure: 'Ubar-col 2k-1' is same node as 'Ubar-col 2k' (both are $b_k$).
From $\tau$: 'U-col r' same as 'Ubar-col $\tau(r)$'.

This creates a graph where vertices are the 2m dummy indices (or sets thereof).
Let's formalize.
Graph with 2m vertices $1..2m$ representing the a_k's and b_k's? No.
Graph with $2m$ vertices representing the column slots on LHS (U) and $2m$ vertices on RHS (Ubar).
LHS slots (2k-1, 2k) are connected (same $a_k$).
RHS slots (2k-1, 2k) are connected (same $b_k$).
$\tau$ connects LHS i to RHS $\tau(i)$.
The number of connected components of this graph is the number of independent variables = number of loops.
Since each component contributes a factor of $d$, the contraction value is $d^{\#components}$.
Note: Total vertices = 4m. Edges = m (LHS pairs) + m (RHS pairs) + 2m (Tau).
Wait, we sum over a, b.
So we essentially multiply the Weingarten weight by $d^{\#loops}$.

Algorithm for loops:
Start with disjoint sets of 4m "points"? No, simpler.
We have 2m LHS slots and 2m RHS slots.
Total 4m slots.
Union-find:
1. Unite (LHS 2k-1, LHS 2k) for k=1..m.
2. Unite (RHS 2k-1, RHS 2k) for k=1..m.
3. Unite (LHS r, RHS $\tau(r)$) for r=1..2m.
Count components. But wait, we are uniting slots to see which $a_k$ and $b_l$ are same.
No, we need to count how many free variables remain.
Variables are $a_1..a_m$ and $b_1..b_m$.
Constraints: $\tau(r) = s \implies a_{\lceil r/2 \rceil} = b_{\lceil s/2 \rceil}$.
So we have 2m "variable nodes" ($a_1..a_m, b_1..b_m$).
For each $r \in 1..2m$:
   u = $\lceil r/2 \rceil$ (index of a)
   v = $\lceil \tau(r)/2 \rceil$ (index of b)
   add edge between $a_u$ and $b_v$.
Count connected components in this bipartite graph using standard DFS/BFS/UnionFind.
Exponent of d is the number of connected components.
"""
function integrate_indices_coe(
    indices::Vector{Tuple{Int,Int}},
    U_bar_indices::Vector{Tuple{Int,Int}},
    dim,
)
    # Re-structure arguments to match integrate_indices signature
    # process_term splits indices into U and U_bar.
    # Here, U_indices correspond to S terms, U_bar_indices correspond to S_bar terms.
    # Note: process_term treats S as U-type and S_bar as U_bar-type.

    n_s = length(indices)
    n_s_bar = length(U_bar_indices)

    if n_s != n_s_bar
        return 0
    end

    m = n_s # This is 'm' in comments above
    n = 2 * m # Total U matrices

    # Unpack indices
    # S_ij -> U_ia U_ja
    # S_pq_bar -> U_pb_bar U_qb_bar

    # Row indices for the 2n U matrices
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

    # We iterate over permutations valid for ROWS
    valid_sigmas = get_matching_permutations(U_rows, U_bar_rows)
    if isempty(valid_sigmas)
        return 0
    end

    n_fact = try
        factorial(n)
    catch
        # Overflow protection if n is huge, though n is small here
        BigInt(1)
    end
    
    # Check if we are summing over the full group
    is_full_group = length(valid_sigmas) == n_fact

    permutations_n = collect(permutations(1:n))

    total = 0 // 1

    if is_full_group
        # Optimization: Factorize sum
        # \sum_{\sigma, \tau} Wg(\sigma \tau^{-1}) d^{loops(\tau)}
        # = (\sum_z Wg(z)) * (\sum_\tau d^{loops(\tau)})
        
        # 1. Sum over Weingarten values (grouped by cycle type)
        sum_wg = 0 // 1
        for p_type in partitions(n)
            # Count permutations with this cycle type
            # formula: n! / (prod k^c_k * c_k!)
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

        # 2. Sum over loop weights (grouped by loop count)
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
        # Fallback: Full double sum
        # Group terms by cycle type to minimize Symbolics overhead
        wg_coeffs = Dict{Vector{Int}, Any}()
        
        for tau in permutations_n
            # Calculate loop weight for columns
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
    integrate_indices_cse(indices, U_bar_indices, dim)

Integration of CSE terms.
Mapping:
```math
S_{ij} = \sum_k U_{ik} (U^R)_{kj} = \sum_{k,x,y} U_{ik} J_{kx} U_{yx} (-J_{yj})
```
Mapping
```math
\bar{S}_{pq} = \sum_{z,w,v} \bar{U}_{pz} \bar{J}_{zw} \bar{U}_{vw} (-\bar{J}_{v_q})
```

Wait, this expansion is getting messy.
Let's use a simpler mapping.
```math
S_{ij} = (U U^R)_{ij}
```
We map monomials in S to monomials in U.
Pairs of indices for U:
S term $k$: $S_{i_k j_k}$.
Becomes:
$U_{i_k a_k}$ (standard U)
$U_{j_k b_k}$ ??? No.
$(U^R)_{a_k j_k}$. This is entry of $U^R$.
```math
(U^R)_{aj} = (J U^T J^T)_{aj} = \sum_{xy} J_{ax} U_{yx} (J^T)_{yj} = \sum_{xy} J_{ax} U_{yx} (-J_{jy})
```
Since $J$ has only one non-zero entry per row, the sums collapse.
$J_{ax}$ is non-zero when $x = pair(a)$ with sign.
```math
(U^R)_{aj} = sign(a) sign(j) U_{pair(j), pair(a)}
```
Where pair(k) maps $k \to k+N$ (or similar) and sign handles the symplectic part.

Let's define the symplectic stricture explicitly.
Assume dim = 2N.
$J = \begin{pmatrix} 0 & I_N \\ -I_N & 0 \end{pmatrix}$.
Indices $1..2N$.
pair(k): if $k \le N$, $k+N$. If $k > N$, $k-N$.
sign(k): if $k \le N$, +1. If $k > N$, -1. (Actually $J_{k, pair(k)}$)
```math
J_{xy} = \delta_{y, pair(x)} \cdot sign(x)
```

Then
```math
(U^R)_{aj} = \sum_{xy} (\delta_{x, pair(a)} sign(a)) U_{yx} (\delta_{y, pair(j)} sign(j) (-1))
```
Check $J^T = -J$, so $(J^T)_{yj} = -J_{yj} = -(\delta_{j, pair(y)} sign(y))$.
Wait, $J_{yj} = -J_{jy}$.
So $(J^T)_{yj} = J_{jy}$.
Let's stick to standard def. $(J^T)_{yj} = J_{jy}$.
$J_{jy} = \delta_{y, pair(j)} sign(j)$.
So
```math
(U^R)_{aj} = \sum_{xy} \delta_{x, pair(a)} sign(a) U_{yx} \delta_{y, pair(j)} sign(j)
```
Substitute x, y:
x = pair(a)
y = pair(j)
Term = $sign(a) sign(j) U_{pair(j), pair(a)}$.

So
```math
S_{ij} = \sum_a U_{ia} (U^R)_{aj} = \sum_a U_{ia} (sign(a) sign(j) U_{pair(j), pair(a)})
```
Variables: a sums from 1 to 2N.
We have product of $U_{i,a}$ and $U_{pair(j), pair(a)}$.
Coefficients: $sign(a) sign(j)$.

Similarly for $\bar{S}_{pq}$.
$\bar{S} = \overline{U U^R} = \bar{U} \overline{U^R} = \bar{U} (U^R)^*$.
$(U^R)^* = (U^R)^T\dagger = (J U^T J^T)^* = J U^\dagger J^T = J \bar{U}^T J^T = \bar{U^R}$.
So
```math
\bar{S}_{pq} = \sum_b \bar{U}_{pb} (\bar{U}^R)_{bq}
```
And assuming $\bar{U}^R$ follows same logic (just with bar):
```math
(\bar{U}^R)_{bq} = sign(b) sign(q) \bar{U}_{pair(q), pair(b)}
```

So mapping:
Each $S_{ij}$ becomes sum over a:
 $(sign(j) sign(a)) \cdot U_{i, a} \cdot U_{pair(j), pair(a)}$.
Each $\bar{S}_{pq}$ becomes sum over b:
 $(sign(q) sign(b)) \cdot \bar{U}_{p, b} \cdot \bar{U}_{pair(q), pair(b)}$.

We need to implement this transformation.
1. Generate the U indices and U_bar indices.
2. Track the coefficients (signs).
3. Sum over dummy variables a, b.

Summation logic is similar to COE:
Integration result is sum over valid sigmas (connecting external row indices) and valid taus (connecting column indices).
Weight: product of signs * sum over dummy a,b of delta-contractions.
Delta contractions link a, pair(a), b, pair(b).
Graph vertices: 2m dummy pairs $(a_k, pair(a_k))$ and $(b_k, pair(b_k))$.
Wait, $a_k$ determines $pair(a_k)$. The variable is just $a_k$.
Terms involving $a_k$:
1. $U$-col at position (2k-1) is $a_k$.
2. $U$-col at position (2k) is $pair(a_k)$.
3. Coefficient $sign(a_k)$.

Terms involving $b_k$:
1. $Ubar$-col at position (2k-1) is $b_k$.
2. $Ubar$-col at position (2k) is $pair(b_k)$.
3. Coefficient $sign(b_k)$.

We iterate over taus (permutations of 1..2n).
Tau connects U-cols to Ubar-cols.
Constraints on loop:
U-col(2k) = pair(U-col(2k-1)).
Ubar-col(2k) = pair(Ubar-col(2k-1)).
Tau connects U-col to Ubar-col.

Also we must sum
```math
\sum_{a_1...} \prod sign(a_k) sign(b_k) \dots
```
This "signed loop" count is tricky.
For COE, it was just $d^{loops}$.
Here, the value of $a_k$ affects the sign.
$sign(x)$ is +1 for $1..N$, -1 for $N+1..2N$.
$pair(x)$ flips between the two halves.
So $sign(pair(x)) = -sign(x)$.

Let's look at the constraints.
We have a graph where nodes are "variables" $a_1..a_m, b_1..b_m$.
Constraint from U-cols: $col(2k) = pair(col(2k-1))$.
Constraint from Ubar-cols: $col(2k) = pair(col(2k-1))$.
Constraint from Tau: $col(r) = col(\tau(r))$ (U to Ubar).

Let's analyze dependency chain.
Start with $a_1$.
Equation 1: $val(U, 2) = pair(val(U, 1)) = pair(a_1)$.
Equation 2: $val(U, 1) = a_1$.
Equation 3: $val(U, r) = val(Ubar, \tau(r))$.
This links a's and b's with 'pair' operations.
$x = pair(y)$ means they are coupled.
Any connected component of variables determines all variables in terms of one free variable (say $x$).
Relations are always of form $x = y$ or $x = pair(y)$.
If we traverse a loop and get $x = x$, it's consistent.
If $x = pair(x)$, it's impossible (since $k \ne k+N$). Component kills the term (0).
If consistent, we have 1 free variable $x \in 1..2N$.
We need to sum $\prod sign(...) $.
The signs depend on a's and b's.
Each $a_k$ contributes $sign(a_k)$. Each $b_k$ contributes $sign(b_k)$.
External signs $sign(j)$ and $sign(q)$ are constant factors.
We need to track how many "sign-flips" (pair operations) we traverse.
Actually, $sign(pair(x)) = -sign(x)$.
So we can express every variable in the component as $x$ or $pair(x)$, i.e., having sign $s_0$ or $-s_0$.
Then we check if the product of signs over the component sums to non-zero.
Sum over $x \in 1..2N$ of $(sign(x))^K$.
If K is even, sum is $2N$.
If K is odd, sum is $\sum sign(x) = (N * 1) + (N * -1) = 0$.
So we just need to count parity of references to 'sign'.

Implementation detail:
Use a DisjointSet with parity?
"Weighted Union Find" or just BFS.
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
        ;
        return 0;
    end

    m = n_s
    n = 2 * m
    # Use phys_dim for symplectic structure J
    # dim is the integration parameter (which might be symbolic)
    half_dim = div(phys_dim, 2) # N

    # 1. Expand Indices
    # S_ij -> sgn(j)*[ U_{i, a}, U_{pair(j), pair(a)}, sgn(a) ]
    # Bar S_pq -> sgn(q)*[ Ubar_{p, b}, Ubar_{pair(q), pair(b)}, sgn(b) ]

    U_rows = Vector{Int}(undef, n)
    # Signs from fixed indices
    fixed_sign_coeff = 1

    # Helper for sign/pair
    get_sign(idx) = (idx <= half_dim ? 1 : -1)
    get_pair(idx) = (idx <= half_dim ? idx + half_dim : idx - half_dim)

    for k = 1:m
        i, j = indices[k]
        # U term 1: U_{i, a} -> row i. col index 2k-1 (map to a_k)
        U_rows[2k-1] = i
        # U term 2: U_{pair(j), pair(a)} -> row pair(j). col index 2k (map to pair(a_k))
        U_rows[2k] = get_pair(j)

        fixed_sign_coeff *= get_sign(j)
    end

    U_bar_rows = Vector{Int}(undef, n)
    for k = 1:m
        p, q = U_bar_indices[k]
        # Ubar term 1: Ubar_{p, b} -> row p. col index 2k-1 (map to b_k)
        U_bar_rows[2k-1] = p
        # Ubar term 2: Ubar_{pair(q), pair(b)} -> row pair(q). col index 2k (map to pair(b_k))
        U_bar_rows[2k] = get_pair(q)

        fixed_sign_coeff *= get_sign(q)
    end

    # Valid Sigmas (Row permutations)
    valid_sigmas = get_matching_permutations(U_rows, U_bar_rows)

    if isempty(valid_sigmas)
        return 0
    end

    permutations_n = collect(permutations(1:n))
    total_val = 0 // 1

    # Group terms by cycle type to minimize Symbolics overhead
    wg_coeffs = Dict{Vector{Int}, Any}()

    for tau in permutations_n
        # Compute Loop Weight for Columns once per tau
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
        # Count components and handle signs
        processed = zeros(Bool, 2*m)
        for i = 1:(2*m)
            if !processed[i]
                root, _ = find_local(i)
                # All members of this component
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

        # Inner loop over sigma: record cycle types
        inv_tau = invperm(tau)
        for sigma in valid_sigmas
            P = [sigma[inv_tau[i]] for i = 1:n]
            cycles = get_cycle_type(P)
            wg_coeffs[cycles] = get(wg_coeffs, cycles, 0) + term_weight
        end
    end

    # Final assembly of symbolic result
    for (cycles, coeff) in wg_coeffs
        total_val += coeff * weingarten(cycles, dim)
    end
    return total_val

    return fixed_sign_coeff * total_val
end

# End of circular ensemble measures
