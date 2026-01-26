# Gaussian Ensembles Integration

This section details the integration of polynomial functions over Gaussian random matrix ensembles: the Gaussian Unitary Ensemble (GUE) and the Gaussian Orthogonal Ensemble (GOE).

## Gaussian Unitary Ensemble (GUE)

### Overview

IntU.jl allows for evaluating integrals over Hermitian random matrices $H$ drawn from the GUE. The measure is defined such that the entries are independent complex Gaussian variables (subject to Hermiticity $H = H^\dagger$).

```julia
measure = dGUE(H, d)
```

### Theory

The integration relies on **Wick's Theorem** for complex Gaussian variables. The fundamental contraction rule is:

```math
\langle H_{ij} H_{kl} \rangle = \delta_{il} \delta_{jk}
```

This normalization implies:
- The variance of off-diagonal entries ($i \neq j$) is $\langle |H_{ij}|^2 \rangle = 1$.
- The variance of diagonal entries is $\langle H_{ii}^2 \rangle = 1$.
- $\langle \text{Tr}(H^2) \rangle = d^2$.

## Gaussian Orthogonal Ensemble (GOE)

### Overview

The GOE consists of real symmetric random matrices $H = H^T$. IntU.jl supports integration over these matrices via:

```julia
measure = dGOE(H, d)
```

### Theory

For the GOE, the contraction rule reflects the real symmetry:

```math
\langle H_{ij} H_{kl} \rangle = \delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}
```

This normalization implies:
- The variance of off-diagonal entries ($i \neq j$) is $\langle H_{ij}^2 \rangle = 1$.
- The variance of diagonal entries is $\langle H_{ii}^2 \rangle = 2$.
- $\langle \text{Tr}(H^2) \rangle = d^2 + d$.

## Implementation Details

IntU.jl automates the following steps:
1.  **Index Collection**: Parses the expression to find all occurrences of $H$.
    - For GUE, `conj(H[i,j])` is treated as $H[j,i]$.
    - For GOE, `conj(H[i,j])` is treated as $H[i,j]$.
2.  **Pair Partitioning**: Generates all ways to pair up the $H$ factors. If the number of factors is odd, the integral is 0.
3.  **Contraction**: For each pair, checks if the contraction is non-zero according to the ensemble-specific rules.
4.  **Summation**: Sums the contributions from all valid pairings.

## Examples

### GUE Integration

```julia
using IntU, Symbolics

@variables d
H = SymbolicMatrix(d, d, :H)
measure = dGUE(H, d)

# < Tr(H^4) > = 2d^3 + d
expr = tr(H^4)
integrate(expr, measure)
# Output: 2d^3 + d
```

### GOE Integration

```julia
using IntU, Symbolics

@variables d
H = SymbolicMatrix(d, d, :H)
measure = dGOE(H, d)

# < Tr(H^2) > = d^2 + d
expr = tr(H^2)
integrate(expr, measure)
# Output: d^2 + d

# < Tr(H^4) > = 2d^3 + 5d^2 + 5d
expr = tr(H^4)
integrate(expr, measure)
# Output: 2d^3 + 5d^2 + 5d
```
