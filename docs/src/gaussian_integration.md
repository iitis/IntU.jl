# Gaussian Unitary Ensemble (GUE) Integration

This section details the integration of polynomial functions over the Gaussian Unitary Ensemble (GUE) of random matrices.

## Overview

IntU.jl allows for evaluating integrals over Hermitian random matrices $H$ drawn from the GUE. The measure is defined such that the entries are independent complex Gaussian variables (subject to Hermiticity $H = H^\dagger$).

```julia
measure = dGUE(H, d)
```

## Theory

The integration relies on **Wick's Theorem** (or Isserlis' theorem) for complex Gaussian variables.

The fundamental contraction rule is:

```math
\langle H_{ij} H_{kl} \rangle = \delta_{il} \delta_{jk}
```

This normalization implies:
- The variance of off-diagonal entries ($i \neq j$) is $\langle |H_{ij}|^2 \rangle = 1$.
- The variance of diagonal entries is $\langle H_{ii}^2 \rangle = 1$.
- $\langle \text{Tr}(H^2) \rangle = d^2$.

(Note: Some conventions scale the variance by $1/d$. In `IntU.jl`, we use the unscaled variance for symbolic simplicity. You can rescale the result by replacing $H \to H/\sqrt{d}$ if needed.)

For a general polynomial, the integral (expectation value) is computed by summing over all **pair partitions** (perfect matchings) of the factors.

```math
\langle H_{i_1 j_1} \dots H_{i_{2k} j_{2k}} \rangle = \sum_{\pi \in \mathcal{P}_{2k}} \prod_{(u, v) \in \pi} \langle H_{i_u j_u} H_{i_v j_v} \rangle
```

where the sum runs over all pairings $\pi$ of the indices $\{1, \dots, 2k\}$.

## Implementation Details

IntU.jl automates the following steps:
1.  **Index Collection**: Parses the expression to find all occurrences of $H$. Note that `conj(H[i,j])` is treated as $H[j,i]$ due to Hermiticity.
2.  **Pair Partitioning**: Generates all ways to pair up the $H$ factors. If the number of factors is odd, the integral is 0.
3.  **Contraction**: For each pair, checks if the contraction is non-zero (i.e., indices match according to the delta functions).
4.  **Summation**: Sums the contributions from all valid pairings.

## Examples

```julia
using IntU, Symbolics

@variables d
H = SymbolicMatrix(d, d, :H)
measure = dGUE(H, d)

# 1. Second Moment
# < Tr(H^2) > = sum_{ij} < H_{ij} H_{ji} > = sum_{ij} 1 = d^2
expr = tr(H^2)
integrate(expr, measure)
# Output: d^2

# 2. Fourth Moment
# < Tr(H^4) > = 2d^3 + d
# (Dominated by planar diagrams 2d^3, plus non-planar crossing term d)
expr = tr(H^4)
integrate(expr, measure)
# Output: 2d^3 + d
```
