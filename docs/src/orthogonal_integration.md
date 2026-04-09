# Orthogonal and Symplectic Integration

This section details the integration of polynomial functions over the real Orthogonal group $O(d)$ and
the Symplectic group $Sp(d)$.

## Orthogonal Group $O(d)$

Integration over $O(d)$ differs from the unitary case because the entries $O_{ij}$ are real.
Consequently, there is no distinction between $U$ and $\bar{U}$, and we integrate monomials
of total even degree.

### Theory: Pair Partitions

The integration formula relies on summing over **Pair Partitions** (also known as matchings) rather than permutations.
Let $M_{2n}$ denote the set of all partitions of the set $\{1, \dots, 2n\}$ into $n$ pairs.
The size of this set is $(2n-1)!!$.

The general formula is [Collins & Śniady, 2006]:

```math
\int_{O(d)} O_{i_1 j_1} \dots O_{i_{2n} j_{2n}} dO = \sum_{\pi_1, \pi_2 \in M_{2n}} \Delta_{\pi_1}(i) \Delta_{\pi_2}(j) \text{Wg}^{O}(\pi_1, \pi_2, d)
```

where:
*   The integrand must have an even total degree ($2n$). If the degree is odd, the integral vanishes.
*   $\Delta_\pi(i)$ is a delta function that equals 1 if the indices $i$ match the pairing $\pi$ (i.e., $i_a = i_b$ for all pairs $\{a, b\} \in \pi$), and 0 otherwise.
*   $\text{Wg}^O(\pi_1, \pi_2, d)$ is the Orthogonal Weingarten function. It depends on the loop structure of the graph formed by superimposing the two matchings $\pi_1$ and $\pi_2$.

### 1. Basic Integration using `@integrate`

The `@integrate` macro automatically identifies `O` as the random orthogonal matrix.

```julia
using IntU, Symbolics
@variables d
# E[O_11^2]
@integrate O[1, 1]^2 dO(d)
# Output: 1/d
```

### 2. Manual Integration

For more control, you can declare symbols explicitly.

```julia
using IntU, Symbolics
@variables d
O = SymbolicMatrix(:O, :O)
integrate(O[1, 1]^2, dO(d))
```

### 3. Higher Moments

The `@integrate` macro can be used for higher moments too.

```julia
using IntU, Symbolics
@variables d
@integrate O[1, 1]^4 dO(d)
# Output: 3 / (d*(d + 2))

# Example O1: high powers, single row
@integrate O[1, 1]^2 * O[1, 2]^4 * O[1, 3]^6 dO(d)

# Example O2: mixed rows/cols
@integrate O[1, 1]^2 * O[2, 2]^4 * O[1, 3]^6 dO(d)
```

## Symplectic Group $Sp(d)$

The Symplectic group $Sp(d)$ (Compact Symplectic Group $USp(d)$) consists of unitary matrices that preserve an antisymmetric bilinear form $J$.
Crucially, the dimension $d$ must be even ($d=2n$).

### The Symplectic Form J

IntU.jl assumes the standard symplectic form $J$:
```math
J = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}
```
Or in tensor block structure. This matrix is used to relate the conjugate entries to the entries themselves:
```math
\bar{S}_{ij} = (J S J^T)_{ij}
```
IntU.jl handles `conj(S[i,j])` by automatically substituting it with the appropriate linear combination of $S$ entries and $J$ factors before integration.

### Theory

The integration formula is analogous to the orthogonal case but involves the symplectic metric.

```math
\int_{Sp(d)} S_{i_1 j_1} \dots S_{i_{2n} j_{2n}} dS = \sum_{\pi_1, \pi_2 \in M_{2n}} \Delta^J_{\pi_1}(i) \Delta^J_{\pi_2}(j) \text{Wg}^{Sp}(\pi_1, \pi_2, d)
```

where $\Delta^J_\pi(i)$ contracts the indices $i$ using the symplectic metric $J$ (introducing signs $\pm 1$).
The Weingarten function $\text{Wg}^{Sp}$ is related to $\text{Wg}^O$ by the dimensional shift $d \to -d$ (and sign adjustments).

### 1. Basic Integration using `@integrate`

The `@integrate` macro is also available for the symplectic group. It identifies `Sp` as the random symplectic matrix.

```julia
using IntU, Symbolics
@variables d
# E[|Sp_11|^2]
@integrate abs(Sp[1, 1])^2 dSp(2)
# Output: 1/2
```

### 2. Manual Integration

```julia
using IntU, Symbolics
@variables d
S = SymbolicMatrix(:S, :Sp)
integrate(abs(S[1, 1])^2, dSp(d))
```

### 3. Higher Moments

```julia
using IntU, Symbolics
@variables d
@integrate abs(Sp[1, 1])^4 dSp(d)
# Output: 2 / (d + d^2)

# Example Sp2: Symplectic mixed moments
@integrate abs(Sp[1, 1])^2 * abs(Sp[1, 2])^2 dSp(d)

# Example Sp1: high powers, single row
@integrate abs(Sp[1, 1])^2 * abs(Sp[1, 2])^4 * abs(Sp[1, 3])^6 dSp(d)
```

## Potential Pitfalls

> [!IMPORTANT]
> ### Symbolic (d) Pitfalls
> - **Small Dimensions**: For Haar-related measures (Unitary, Orthogonal, Circular), element-wise results are rational functions with poles at small $d$ (typically $d < n$ for degree $n$ moments). Pure trace moments $|\mathrm{tr}(U)|^{2k}$ are an exception: they depend on $d$ as a step function and require a concrete integer dimension.
> - **Removable Singularities**: Substituting numeric values can yield $0/0$ forms (e.g., at $d=1$ or $d=2$).
> - **Automatic Handling**: `IntU.jl`'s `evaluate` function automatically simplifies expressions to resolve removable singularities when a denominator evaluates to zero.

- **High-Performance Symbolic Solver**: For symbolic dimensions $d$, `IntU.jl` uses a specialized univariate polynomial solver based on the **Bareiss algorithm**.
- **Exact Rational Summation**: Orthogonal integration uses an internal exact rational arithmetic engine to compute Weingarten sums.
- **Cycle Type Grouping**: The integration engine automatically groups pair partitions by their loop cycle types, reducing the number of symbolic operations.

## References

1.  **Collins, B., & Śniady, P. (2006).** Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups. *Communications in Mathematical Physics*, 264(3), 773-795.
2.  **Matsuki, T.** (1990). The orbits of affine symmetric spaces under the action of parabolic subgroups. *Hiroshima Mathematical Journal*. (Relevant for O/Sp symmetry groups).

## See Also

- [Unitary Integration](unitary_integration.md) — Weingarten calculus for $U(d)$; shares the same index-contraction framework
- [Circular Ensembles](circular_ensembles.md) — COE arises from $O(d)$, CSE from $Sp(d)$
- [Asymptotic Expansions](asymptotic.md) — large-$d$ expansions of orthogonal/symplectic results
- [Integral Library](integral_library.md) — targeted pre-computed moments for $O(d)$ and $Sp(d)$
