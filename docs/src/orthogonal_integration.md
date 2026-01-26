# Orthogonal and Symplectic Integration

This section details the integration over the real Orthogonal group $O(d)$ and
the Symplectic group $Sp(d)$.

## Orthogonal Group $O(d)$

Integration over $O(d)$ is similar to $U(d)$ but involves only real indices
(since $O_{ij}$ are real entries).

### Formula

```math
\int_{O(d)} O_{i_1 j_1} \dots O_{i_{2n} j_{2n}} dO = \sum_{\pi \in M_{2n}} \delta_{i, \pi} \delta_{j, \pi} \text{Wg}^{O}(\pi, d)
```

where the integrand must have an even total degree ($2n$) (if odd, the integral is 0), $M_{2n}$ is the set of pair partitions (pairings) of $2n$ elements, and the Weingarten function $\text{Wg}^O$ differs from the Unitary one.

### Usage

```julia
using IntU, Symbolics
@variables d
@variables O_mat[1:d, 1:d]::Real
measure_O = dO(O_mat, d)

expr = O_mat[1,1]^2
integrate(expr, measure_O)
# Output: 1/d
```

## Symplectic Group $Sp(d)$

The Symplectic group preserves the symplectic form defined by the matrix $J$:

```math
J = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}
```

where $d=2n$.

**Important**: $d$ must be even.

### Formula

```math
\int_{Sp(d)} S_{i_1 j_1} \dots S_{i_{2n} j_{2n}} dS = \sum_{\pi \in M_{2n}} \mathcal{J}(i, \pi) \mathcal{J}(j, \pi) \text{Wg}^{Sp}(\pi, d)
```

where $\mathcal{J}$ denotes the contraction of indices with the symplectic
metric $J$.

### Usage

```julia
using IntU, Symbolics
@variables d
@variables S_mat[1:d, 1:d]::Complex
measure_S = dSp(S_mat, d)

# |S_{1,1}|^2 integration
expr = abs(S_mat[1,1])^2
integrate(expr, measure_S)
# Output: 1/d
```

## Implementation Details & Pitfalls

- **Pair Partitions**: Instead of permutations $\sigma, \tau$, we sum over pair
  partitions. This means the combinatorics grow differently (faster) than
  $U(d)$.
- **Dimension Parity**: For $Sp(d)$, $d$ is implicitly assumed to be even.
- **Metric $J$**: The definitions of indices in $Sp(d)$ integration heavily
  rely on the antisymmetric matrix $J$. Ensure your manual index checks align
  with the standard block structure (tensor identity):
  
```math
J = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}
```

## References

- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar
  measure on unitary, orthogonal and symplectic groups. *Communications in
  Mathematical Physics*.
