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

## Usage

You can define the Gaussian measures using `dGUE`, `dGOE`, and `dGSE`.

### GUE Example

```julia
using IntU, Symbolics

@variables d
# GUE Measure with symbolic dimension
H = SymbolicMatrix(:H)
measure_GUE = dGUE(H, d)

# Average Trace of H^2
# < Tr(H^2) > = d^2
expr = IntU.tr(H^2)
res = integrate(expr, measure_GUE)
println(res)
# Output: d^2
```

### GOE Example

```julia
# GOE Measure
measure_GOE = dGOE(H, d)

# Average Trace of H^2
# < Tr(H^2) > = d^2 + d
res_goe = integrate(IntU.tr(H^2), measure_GOE)
println(res_goe)
# Output: d^2 + d
```

### GSE Example

```julia
# GSE Measure
measure_GSE = dGSE(H, d)

# Average Trace of H^2
# < Tr(H^2) > = d^2 - d
res_gse = integrate(IntU.tr(H^2), measure_GSE)
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
1.  **Index Collection**: Parses the expression to find all occurrences of $H$.
2.  **Pair Partitioning**: Generates all ways to pair up the $H$ factors ($\sim (2k-1)!!$ terms).
3.  **Contraction**: For each pair, applies the specific ensemble contraction rule (GUE, GOE, or GSE).
4.  **Summation**: Sums the contributions.

## References

1.  **Mehta, M. L.** (2004). *Random Matrices*. Elsevier.
2.  **Livan, G., Novaes, M., & Vivo, P.** (2018). *Introduction to Random Matrices: Theory and Practice*. Springer.
3.  **Wick, G. C.** (1950). The evaluation of the collision matrix. *Physical Review*, 80(2), 268.

## Pre-computed Moments

For common moments like $\langle \text{Tr}(H^2) \rangle$, $\langle \text{Tr}(H^4) \rangle$, and $\langle \text{Tr}(H^6) \rangle$, `IntU.jl` uses a [Pre-computed Integral Library](integral_library.md) to provide results instantly.
