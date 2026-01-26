# Unitary Group Integration

This section details the integration of polynomial functions over the unitary
group $U(d)$ with respect to the Haar measure.

## Overview

IntU.jl allows for evaluating integrals of the form:

```math
\int_{U(d)} U_{i_1 j_1} \dots U_{i_n j_n} \bar{U}_{k_1 l_1} \dots \bar{U}_{k_n l_n} dU
```

The result is expressed in terms of the dimension $d$ and Kronecker deltas
matching the indices.

## Theory

The integration relies on the **Weingarten Calculus**. The general formula for
the integral is given by:

```math
\int_{U(d)} U_{i_1 j_1} \dots U_{i_n j_n} \bar{U}_{k_1 l_1} \dots \bar{U}_{k_n l_n} dU = \sum_{\sigma, \tau \in S_n} \delta_{i, k_\sigma} \delta_{j, l_\tau} \text{Wg}(\sigma \tau^{-1}, d)
```

where $S_n$ is the symmetric group of degree $n$ (permutations of $n$ elements),
$\delta_{i, k_\sigma} = \prod_{m=1}^n \delta_{i_m, k_{\sigma(m)}}$ is the
contraction of row indices, $\delta_{j, l_\tau} = \prod_{m=1}^n \delta_{j_m,
l_{\tau(m)}}$ is the contraction of column indices, and $\text{Wg}(\pi, d)$ is
the **Weingarten function**, which depends only on the cycle structure of the
permutation $\pi$ and the dimension $d$.

## Implementation Details

IntU.jl automates the following steps:
1.  **Index Identification**: It parses the symbolic expression to identify
    which variables correspond to elements of $U$ and $\bar{U}$, extracting
    their indices ($i, j, k, l$).
2.  **Degree Matching**: It checks if the number of $U$ factors matches the
    number of $\bar{U}$ factors. If they differ ($n \neq m$), the integral
    vanishes (returns 0) due to phase invariance.
3.  **Symbolic Summation**: It generates the sum over permutations
    $\sigma, \tau \in S_n$, computing the Kronecker delta products symbolically.
4.  **Weingarten Evaluation**: It computes the values of $\text{Wg}(\pi, d)$
    using character theory (via `Combinatorics.jl` and `murnaghan_nakayama`
    rule).

## Potential Pitfalls

-   **Symbolic vs Numeric Dimension**: The dimension $d$ can be symbolic
    (`@variables d`). However, for the Weingarten function to be well-defined,
    $d$ must typically be "large enough" ($d \ge n$). IntU.jl generally provides
    the rational function form valid for large $d$. Singularities may occur if
    you substitute small integer values for $d$ into the final symbolic result
    (poles of the Weingarten function).
-   **Computational Complexity**: The sum involves $(n!)^2$ terms. While optimized
    to group cycles, integrals with high degrees ($n > 6$) can become
    computationally expensive.
-   **Ambiguity with Conjugates**: Ensure you use `conj()` correctly. In Julia,
    `A'` is the conjugate transpose. IntU.jl handles `conj(U[i,j])` as
    $\bar{U}_{ij}$.

## Examples

```julia
using IntU, Symbolics

@variables d
@variables U[1:d, 1:d]::Complex
measure = dU(U, d)

# 1. Norm of a matrix element
expr = abs(U[1,1])^2 # U[1,1] * conj(U[1,1])
integrate(expr, measure)
# Output: 1/d

# 2. Trace moments (large d)
tr_U = IntU.tr(U)
expr = abs(tr_U)^4
integrate(expr, measure)
# Output: 2 (as d -> Infinity, converges to Gaussian moment)
```

## References
- Collins, B. (2003). Moments and Cumulants of Polynomial random variables on
  unitary groups, the Itzykson-Zuber integral and free probability.
  *International Mathematics Research Notices*.
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar
  measure on unitary, orthogonal and symplectic groups. *Communications in
  Mathematical Physics*.
