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

> [!NOTE]
> **Special Unitary Group $SU(d)$**: For all currently supported "balanced" polynomial expressions (where the number of $U$ and $\bar{U}$ factors are equal), the integration over $SU(d)$ is equivalent to $U(d)$. Non-stable-range effects involving $\epsilon$-tensors for specific small $d$ are not currently covered. In the current backend, non-balanced queries are evaluated with the same phase-invariance rule as $U(d)$ and therefore return `0`.

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
    *   Transposition $\text{Wg}([2], d) = -\frac{1}{d(d^2-1)}$
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
This is achieved through the `SymbolicMatrix` type, which represents a matrix
of arbitrary (symbolic) size.

## Unitary Designs
 
In many applications, full Haar integration is not required. Instead, one uses **Unitary $t$-designs**, which are ensembles that mimic the first $t$ moments of the Haar measure. A unitary $t$-design is a probability distribution $\nu$ such that for any polynomial $P$ of degree at most $t$ in $U$ and $\bar{U}$:
```math
\mathbb{E}_{U \sim \nu} [P(U, \bar{U})] = \mathbb{E}_{U \sim \text{Haar}} [P(U, \bar{U})]
```
 
### Usage: dDesign
 
Use the `dDesign(d, t)` measure to define a unitary $t$-design of dimension $d$ and order $t$.
 
```julia
using IntU, Symbolics
@variables d
# E[|U_11|^2] using a 2-design (degree 1 in U, 1 in U*)
@integrate abs(U[1,1])^2 dDesign(d, 2)
# Output: 1/d
```
 
### Integration Behavior and Guards
 
`IntU.jl` enforces validity of integration over $t$-designs with the
following behaviour:

1. **Balanced integrands** (equal degree $q$ in $U$ and $U^\dagger$):
   - $q \le t$: returns the exact Haar-averaged result via Weingarten calculus.
   - $q > t$: throws an `ArgumentError` (e.g., `Integrand degree (3) exceeds design order t=2`), preventing accidental reliance on non-guaranteed values.
2. **Unbalanced integrands** (different degree in $U$ and $U^\dagger$):
   returns zero, which is the Haar-correct value. This is guaranteed correct
   for total degree $\le t$; for higher-degree unbalanced monomials, a
   $t$-design does not guarantee this value, but no error is raised.

This guard mechanism helps ensure that physical simulations and protocol
verifications (such as randomized benchmarking) remain mathematically
rigorous for balanced polynomial integrands.

## Examples

### 1. Basic Integration using `@integrate`

The `@integrate` macro provides a convenient way to integrate expressions without manually declaring variables. It uses heuristics to identify random unitaries (usually `U`) and dimensions.

```julia
using IntU, Symbolics
@variables d
# E[|U_11|^2]
@integrate abs(U[1, 1])^2 dU(d)
# Output: 1/d

# High-degree moments (Example A: different rows/cols)
res_a = @integrate abs(U[1, 1])^2 * abs(U[2, 2])^4 * abs(U[1, 3])^6 dU(d)

# High-degree moments (Example B: same row, Dirichlet behavior)
res_b = @integrate abs(U[1, 1])^2 * abs(U[1, 2])^4 * abs(U[1, 3])^6 dU(d)

# Example C: "2-design style" correlator identity
res_c = @integrate U[1, 1] * conj(U[1, 2]) * U[2, 2] * conj(U[2, 1]) dU(d)
```

### 2. Manual Integration

For more control, or when dealing with multiple matrices, you can declare symbols explicitly.


```julia
using IntU, Symbolics
@variables d
U = SymbolicMatrix(:U, :U)

# 1. Norm of a matrix element
# Integral of |U_{11}|^2
integrate(abs(U[1, 1])^2, dU(d))
# Output: 1/d
```

### 3. Higher Unitary Moments

```julia
using IntU, Symbolics
@variables d
U = SymbolicMatrix(:U, :U)

# 2. Fourth moment
# Integral of |U_{11}|^4
integrate(abs(U[1, 1])^4, dU(d))
# Output: 2 / (d*(d + 1))
```

### 4. Trace Moments

```julia
using IntU, Symbolics
@variables d
U = SymbolicMatrix(:U, :U)

# 3. Trace moments
# Integral of |Tr(U)|^2
integrate(abs(tr(U))^2, dU(d))
# Output: 1
```

### 5. Matrix Integration

You can integrate matrix-valued expressions directly. The function `integrate` will element-wise integrate any `AbstractArray` (including `SymbolicMatrix` and `SymbolicMatrixProduct`) passed to it.

```julia
using IntU, Symbolics
@variables d
U = SymbolicMatrix(:U, :U)
# E[tr(U A U' B)] = tr(A) * tr(B) / d
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)
integrate(tr(U * A * U' * B), dU(d))
```

### 6. Complex Trace Powers

`IntU.jl` supports integrating complex powers and absolute values of traces, such as $|tr(U)|^n$. These are processed efficiently using a lazy evaluation engine.

```julia
using IntU, Symbolics
@variables d
U = SymbolicMatrix(:U, :U)

# Integral of |Tr(U)|^4
integrate(abs(tr(U))^4, dU(d))
# Output: 2
```

### 7. HCIZ Integrals

IntU.jl provides direct support for **Harish-Chandra-Itzykson-Zuber (HCIZ)** integrals.

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

> [!IMPORTANT]
> ### Symbolic (d) Pitfalls
> - **Small Dimensions**: For Haar-related measures (Unitary, Orthogonal, Circular), results are rational functions with poles at small $d$ (typically $d < n$ for degree $n$ moments).
> - **Removable Singularities**: Substituting numeric values can yield $0/0$ forms (e.g., at $d=1$ or $d=2$).
> - **Automatic Handling**: `IntU.jl`'s `evaluate` function automatically simplifies expressions to resolve removable singularities when a denominator evaluates to zero.

- **Computational Complexity**: The sum involves $(n!)^2$ terms. While optimized
    to group cycles, integrals with high degrees ($n > 6$) can become
    computationally expensive. The number of terms grows factorially.

## References

1.  **Collins, B. (2003).** Moments and Cumulants of Polynomial random variables on unitary groups, the Itzykson-Zuber integral and free probability. *International Mathematics Research Notices*, 2003(17), 953-982.
2.  **Collins, B., & Śniady, P. (2006).** Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups. *Communications in Mathematical Physics*, 264(3), 773-795. [arXiv:math-ph/0402073](https://arxiv.org/abs/math-ph/0402073)
3.  **Puchala, Z., & Miszczak, J. A. (2017).** Symbolic integration with respect to the Haar measure on the unitary groups. *Bulletin of the Polish Academy of Sciences. Technical Sciences*, 65(1), 21-27.

## See Also

- [Symbolic Trace Logic](symbolic_trace.md) — index-free trace integration built on top of this engine
- [Orthogonal & Symplectic](orthogonal_integration.md) — Weingarten calculus for $O(d)$ and $Sp(d)$
- [Circular Ensembles](circular_ensembles.md) — COE, CUE, CSE using the same Haar machinery
- [Integral Library](integral_library.md) — pre-computed results that bypass this engine
- [Asymptotic Expansions](asymptotic.md) — large-$d$ expansions of Weingarten results
- [Pure States](pure_states.md) — first-column Weingarten specialisation
- [Stiefel Manifolds](stiefel_manifold.md) — first-$k$-column Weingarten specialisation
