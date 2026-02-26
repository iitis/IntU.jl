# Gaussian Ensembles Integration

# Gaussian Random Matrix Integration

This section details the integration of polynomial functions over Gaussian Random Matrix Ensembles:
the Gaussian Unitary Ensemble (GUE), the Gaussian Orthogonal Ensemble (GOE), and the Gaussian Symplectic Ensemble (GSE).

## Theory: Wick's Theorem

Unlike Haar measure integration which requires Weingarten calculus, Gaussian integration relies on **Wick's Theorem** (or Isserlis' Theorem).
The integral of a product of centered Gaussian random variables is given by the sum over all possible pair contractions (matchings).

### GUE (Gaussian Unitary Ensemble)

$H$ is a complex Hermitian matrix ($H = H^\dagger$). The entries are independent complex Gaussian variables (subject to Hermiticity).
The contraction rule is:
```math
\langle H_{ij} \bar{H}_{kl} \rangle_{GUE} = \delta_{il} \delta_{jk}
```
This effectively "connects" the indices in a specific way corresponding to the unitary symmetry.

### GOE (Gaussian Orthogonal Ensemble)

$H$ is a real symmetric matrix ($H = H^T$). The entries are real Gaussian variables.
The contraction rule includes an extra term due to symmetry ($H_{kl} = H_{lk}$):
```math
\langle H_{ij} H_{kl} \rangle_{GOE} = \delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}
```

### GSE (Gaussian Symplectic Ensemble)

$H$ is a Hermitian quaternionic self-dual matrix ($d=2n$).
The integrals relate to GOE via specific duality relations (often mapping $d \to -d$ or $d \to 2d$).
The contraction rule involves the symplectic form $J$:
```math
\langle H_{ij} H_{kl} \rangle_{GSE} = \delta_{il} \delta_{jk} + (J)_{ik} (J)_{jl}
```
IntU.jl implements GSE integration by mapping it to contractions involving the definition of the symplectic metric.

## Ginibre Ensembles

Ginibre ensembles consist of non-Hermitian matrices where each entry is an independent Gaussian random variable.

### GinUE (Complex Ginibre Ensemble)

Matrices $G$ with i.i.d. complex Gaussian entries. The contraction rule is:
```math
\langle G_{ij} \bar{G}_{kl} \rangle_{GinUE} = \delta_{ik} \delta_{jl}
```
Note that only contractions between $G$ and its complex conjugate $\bar{G}$ are non-zero.

### GinOE (Real Ginibre Ensemble)

Matrices $G$ with i.i.d. real Gaussian entries. The contraction rule is:
```math
\langle G_{ij} G_{kl} \rangle_{GinOE} = \delta_{ik} \delta_{jl}
```

### GinSE (Symplectic Ginibre Ensemble)

Matrices $G$ with i.i.d. quaternionic Gaussian entries. Integrals are computed using duality relations.

## Usage

You can define the Gaussian measures using `dGUE`, `dGOE`, `dGSE`, `dGinUE`, `dGinOE`, and `dGinSE`.

### GUE Example

```julia
using IntU, Symbolics

# GUE Measure with symbolic dimension
# Average Trace of H^2
# < Tr(H^2) > = d^2
res = @integrate tr(H^2) dGUE(d)
println(res)
# Output: d^2
```

### GinUE Example

```julia
# GinUE Measure with symbolic dimension
# Average Trace of G G'
# < Tr(G G') > = d^2
res_ginue = @integrate tr(G * G') dGinUE(d)
println(res_ginue)
# Output: d^2
```

### GOE Example

```julia
# GOE Measure
# Average Trace of H^2
# < Tr(H^2) > = d^2 + d
res_goe = @integrate tr(H^2) dGOE(d)
println(res_goe)
# Output: d^2 + d
```

### GSE Example

```julia
# GSE Measure
# Average Trace of H^2
# < Tr(H^2) > = d^2 - d
res_gse = @integrate tr(H^2) dGSE(d)
println(res_gse)
# Output: d^2 - d
```

## Scaling Conventions

IntU.jl computes the raw Gaussian moments (combinatorial counts).
This corresponds to the normalization where the variance of off-diagonal entries is 1.

*   **Standard Physics Normalization** (Wigner's semicircle law radius 2): Requires scaling $H \to H/\sqrt{d}$.
    $\langle \text{Tr}(H^2) \rangle \to d$.
*   **IntU.jl Normalization**:
    $\langle \text{Tr}(H^2) \rangle_{GUE} = d^2$.

## Implementation Details

IntU.jl automates the following steps:
1.  **Index Collection**: Parses the expression to find all occurrences of $H$ (or $G$).
2.  **Pair Partitioning**: Generates all ways to pair up the matrix factors.
3.  **Contraction**: For each pair, applies the specific ensemble contraction rule.
4.  **Summation**: Sums the contributions.

## References

1.  **Mehta, M. L.** (2004). *Random Matrices*. Elsevier.
2.  **Livan, G., Novaes, M., & Vivo, P.** (2018). *Introduction to Random Matrices: Theory and Practice*. Springer.
3.  **Wick, G. C.** (1950). The evaluation of the collision matrix. *Physical Review*, 80(2), 268.
4.  **Ginibre, J.** (1965). Statistical ensembles of complex, real, and quaternionic matrices. *Journal of Mathematical Physics*, 6(3), 440-449.

## Pre-computed Moments

For common moments like $\langle \text{Tr}(H^2) \rangle$, $\langle \text{Tr}(H^4) \rangle$, and $\langle \text{Tr}(H^6) \rangle$, `IntU.jl` uses a [Pre-computed Integral Library](integral_library.md) to provide results instantly.
