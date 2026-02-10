# Unitary Group Integration

This section details the integration of polynomial functions over the unitary
group $U(d)$ with respect to the Haar measure. It covers the theoretical foundations
based on Weingarten Calculus and provides practical examples using `IntU.jl`.

## Overview

IntU.jl allows for evaluating integrals of the form:

```math
\int_{U(d)} U_{i_1 j_1} \dots U_{i_n j_n} \bar{U}_{k_1 l_1} \dots \bar{U}_{k_n l_n} dU
```
where $U$ is a $d \times d$ unitary matrix ($U U^\dagger = I$). The integral is taken with respect
to the Haar measure, which is the unique translation-invariant probability measure on the compact
Lie group $U(d)$.

The result is expressed in terms of the dimension $d$ and Kronecker deltas matching the indices.

## Theory: Weingarten Calculus

The integration relies on the **Weingarten Calculus**, a combinatorial method for evaluating
integrals over compact groups.

### The Integration Formula

The general formula for the integral is given by [Collins & Śniady, 2006]:

```math
\int_{U(d)} U_{i_1 j_1} \dots U_{i_n j_n} \bar{U}_{k_1 l_1} \dots \bar{U}_{k_n l_n} dU = \sum_{\sigma, \tau \in S_n} \delta_{i, k_\sigma} \delta_{j, l_\tau} \text{Wg}(\sigma \tau^{-1}, d)
```

where:
*   $S_n$ is the symmetric group of permutations of $\{1, \dots, n\}$.
*   $\delta_{i, k_\sigma} = \prod_{m=1}^n \delta_{i_m, k_{\sigma(m)}}$ contracts row indices according to $\sigma$.
*   $\delta_{j, l_\tau} = \prod_{m=1}^n \delta_{j_m, l_{\tau(m)}}$ contracts column indices according to $\tau$.
*   $\text{Wg}(\pi, d)$ is the **Weingarten function**, which depends only on the cycle structure of the permutation $\pi$ and the dimension $d$.

### The Weingarten Function

The Weingarten function $\text{Wg}(\pi, d)$ is a rational function of $d$. It is defined via
the inverse of the Gram matrix of Schur functions or via character theory.

For small $n$, the values are:
*   **n=1**: $\text{Wg}([1], d) = \frac{1}{d}$
*   **n=2**:
    *   Identity $\text{Wg}([1,1], d) = \frac{1}{d^2-1}$
    *   Transposition $\text{Wg}([2], d) = -\frac{1}{d(d^2-1)$
*   **n=3**:
    *   $\text{Wg}([1,1,1], d) = \frac{d^2-2}{d(d^2-1)(d^2-4)}$
    *   $\text{Wg}([2,1], d) = -\frac{1}{(d^2-1)(d^2-4)}$
    *   $\text{Wg}([3], d) = \frac{2}{d(d^2-1)(d^2-4)}$

IntU.jl computes these values symbolically for any $n$, using the Murnaghan-Nakayama rule
to evaluate characters of the symmetric group.

## Implementation Details

IntU.jl automates the following steps:
1.  **Index Identification**: It parses the symbolic expression to identify
    which variables correspond to elements of $U$ and $\bar{U}$, extracting
    their indices ($i, j, k, l$).
2.  **Degree Matching**: It checks if the number of $U$ factors matches the
    number of $\bar{U}$ factors. If they differ ($n \neq m$), the integral
    vanishes (returns 0) due to phase invariance ($U \to e^{i\theta}U$).
3.  **Symbolic Summation**: It generates the sum over permutations
    $\sigma, \tau \in S_n$, computing the Kronecker delta products symbolically.
4.  **Weingarten Evaluation**: It computes the values of $\text{Wg}(\pi, d)$
    using `Combinatorics.jl` characters.

### Symbolic Dimension

A key feature of IntU.jl is the ability to leave the dimension $d$ as a symbolic variable.
This is achieved through the `SymbolicUnitary` type, which represents a Unitary matrix
of arbitrary (symbolic) size.

The macro `@symbolic_dimension` facilitates the creation of such matrices.

## Examples

### 1. Basic Integration

```julia
using IntU, Symbolics

# Define a unitary matrix U of symbolic dimension d
@variables d
@symbolic_dimension U[1:d, 1:d]
measure = dU(U)

# 1. Norm of a matrix element
# Integral of |U_{11}|^2
expr = abs(U[1,1])^2 
res = integrate(expr, measure)
println(res)
# Output: 1/d
```

### 2. Higher Unitary Moments

```julia
# 2. Fourth moment
# Integral of |U_{11}|^4
expr2 = abs(U[1,1])^4
res2 = integrate(expr2, measure)
println(res2)
# Output: 2 / (d*(d + 1))
```

### 3. Trace Moments

```julia
# 3. Trace moments
tr_U = IntU.tr(U)
# Integral of |Tr(U)|^2
res3 = integrate(abs(tr_U)^2, measure)
println(res3)
# Output: 1
```

### 4. Matrix Integration

New in v0.2: You can integrate matrix-valued expressions directly. The function `integrate` will element-wise integrate any `AbstractArray` passed to it.

```julia
using LinearAlgebra

# Integrate U * U' (should be identity)
# We collect the symbolic expression to a Matrix{Num} to ensure it's treated as an array of expressions
expr_mat = collect(U * U')
res_mat = integrate(expr_mat, measure)
# Result is the Identity matrix
```

### 5. HCIZ Integrals

New in v0.2: Direct support for **Harish-Chandra-Itzykson-Zuber (HCIZ)** integrals.

```julia
using IntU, LinearAlgebra

# Define matrices
A = diagm([1.0, 2.0])
B = diagm([0.5, 1.5])

# Compute ∫ dU exp(Tr(A U B U'))
res = hciz(A, B)
println(res)
# Output: 20.9329...
```

See the [API Reference](api.md) for more details.

## Potential Pitfalls

-   **Symbolic vs Numeric Dimension**: The dimension $d$ can be symbolic.
    However, the Weingarten function has poles at small integers ($d < n$).
    The symbolic result assumes $d$ is generic/large. Substituting discrete values $d < n$
    into the rational function may result in division by zero, although the integral itself
    is well-defined.
-   **Computational Complexity**: The sum involves $(n!)^2$ terms. While optimized
    to group cycles, integrals with high degrees ($n > 6$) can become
    computationally expensive. The number of terms grows factorially.

## References

1.  **Collins, B. (2003).** Moments and Cumulants of Polynomial random variables on unitary groups, the Itzykson-Zuber integral and free probability. *International Mathematics Research Notices*, 2003(17), 953-982.
2.  **Collins, B., & Śniady, P. (2006).** Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups. *Communications in Mathematical Physics*, 264(3), 773-795. [arXiv:math-ph/0402073](https://arxiv.org/abs/math-ph/0402073)
3.  **Puchala, Z., & Miszczak, J. A. (2017).** Symbolic integration with respect to the Haar measure on the unitary groups. *Bulletin of the Polish Academy of Sciences. Technical Sciences*, 65(1), 21-27.
