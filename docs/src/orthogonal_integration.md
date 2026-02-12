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

### Usage

```julia
using IntU, Symbolics

# Define dimension
@variables d
# Define Orthogonal matrix
O = SymbolicMatrix(:O, :O)
# Integrate squared element
@integrate O[1,1]^2 dO(d)
println(res)
# Output: 1/d

# 2. Fourth power
@integrate O[1,1]^4 dO(d)
println(res4)
# Output: 3 / (d*(d + 2))
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

### Usage

```julia
using IntU, Symbolics

@variables d
# Define Symplectic matrix (Complex entries, unlike O(d))
S = SymbolicMatrix(:S, :Sp)

# |S_{1,1}|^2 integration
@integrate abs(S[1,1])^2 dSp(d)
println(res)
# Output: 1/d

# |S_{1,1}|^4 integration
@integrate abs(S[1,1])^4 dSp(d)
println(res4)
# Output: 2 / ((d + 1)*(d - 2))
```

## Implementation Details & Pitfalls

- **Automatic Conjugation**: For $Sp(d)$, the code treats `conj(S)` non-trivially. It uses the relation $\bar{S} = -J S J$ to rewrite conjugate entries in terms of $S$ entries (and J factors). This allows using the efficient Weingarten formula for products of $S$ only.
- **Dimension Parity**: For $Sp(d)$, $d$ must be even. The symbolic result is valid for even $d$.
- **Computational Complexity**: The sum over pair partitions $M_{2n}$ grows as $(2n)!!$. This is significantly faster than $(n!)^2$ for $U(d)$ but still grows combinatorially.

## References

1.  **Collins, B., & Śniady, P. (2006).** Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups. *Communications in Mathematical Physics*, 264(3), 773-795.
2.  **Matsuki, T.** (1990). The orbits of affine symmetric spaces under the action of parabolic subgroups. *Hiroshima Mathematical Journal*. (Relevant for O/Sp symmetry groups).
