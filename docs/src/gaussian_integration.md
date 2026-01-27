# Gaussian Ensembles Integration

This section details the integration of polynomial functions over Gaussian random matrix ensembles: the Gaussian Unitary Ensemble (GUE), the Gaussian Orthogonal Ensemble (GOE), and the Gaussian Symplectic Ensemble (GSE).

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

## Gaussian Symplectic Ensemble (GSE)

### Overview

The GSE consists of quaternionic Hermitian random matrices $H$. In complex notation, these are $d \times d$ matrices where $d = 2n$, satisfying $H = H^\dagger$ and $H = -J H^T J$. IntU.jl enables integration over these via:

```julia
measure = dGSE(H, d)
```

### Theory

The contraction rule for GSE involves the symplectic form $J$:

```math
\langle H_{ij} H_{kl} \rangle = \delta_{il} \delta_{jk} + (J)_{ik} (J)_{jl}
```

This normalization implies:
- $\langle \text{Tr}(H^2) \rangle = d^2 - d$.
- The ensemble follows the **$d \to -d$ duality** with GOE. Specifically:
  $$\langle \text{Tr}(H^k) \rangle_{GSE}(d) = (-1)^{k/2 + 1} \langle \text{Tr}(H^k) \rangle_{GOE}(-d)$$

### Implementation Details

IntU.jl automates the following steps:
1.  **Index Collection**: Parses the expression to find all occurrences of $H$.
    - For GUE, `conj(H[i,j])` is treated as $H[j,i]$.
    - For GOE, `conj(H[i,j])` is treated as $H[i,j]$.
    - For GSE, `conj(H[i,j])` is treated as $H[j,i]$ (since it is Hermitian).
2.  **Pair Partitioning**: Generates all ways to pair up the $H$ factors. If the number of factors is odd, the integral is 0.
3.  **Contraction**: For each pair, checks if the contraction is non-zero according to the ensemble-specific rules.
4.  **Summation**: Sums the contributions from all valid pairings.

## Examples

### GUE Integration

```julia
using IntU, Symbolics

@variables d
H = SymbolicMatrix(:H)
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
H = SymbolicMatrix(:H)
measure = dGOE(H, d)

# < Tr(H^2) > = d^2 + d
expr = tr(H^2)
integrate(expr, measure)
# Output: d^2 + d

# <Tr(H^4) > = 2d^3 + 5d^2 + 5d
expr = tr(H^4)
integrate(expr, measure)
# Output: 2d^3 + 5d^2 + 5d
```

### GSE Integration

```julia
using IntU, Symbolics

@variables d
H = SymbolicMatrix(:H)
measure = dGSE(H, d)

# < Tr(H^2) > = d^2 - d
expr = tr(H^2)
integrate(expr, measure)
# Output: d^2 - d

# < Tr(H^4) > = 2d^3 - 5d^2 + 5d
expr = tr(H^4)
integrate(expr, measure)
# Output: 2d^3 - 5d^2 + 5d
```

## Pre-computed Moments

For common moments like $\langle \text{Tr}(H^2) \rangle$, $\langle \text{Tr}(H^4) \rangle$, and $\langle \text{Tr}(H^6) \rangle$, `IntU.jl` uses a [Pre-computed Integral Library](integral_library.md) to provide results instantly.

## References

- Mehta, M. L. (2004). *Random Matrices* (Vol. 142). Elsevier.
- Anderson, G. W., Guionnet, A., & Zeitouni, O. (2010). *An Introduction to Random Matrices* (No. 118). Cambridge University Press.
- Forrester, P. J. (2010). *Log-Gases and Random Matrices* (LMS-34). Princeton University Press.
- Wick, G. C. (1950). The evaluation of the collision matrix. *Physical Review*, 80(2), 268.
