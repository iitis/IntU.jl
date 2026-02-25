# Diagonal Unitary Matrices (Torus Group)

This section describes integration over the group of diagonal unitary matrices,
also known as the torus group $T^d$.

## Overview

A diagonal unitary matrix $D$ of dimension $d$ has the form:

```math
D = \text{diag}(e^{i\theta_1}, e^{i\theta_2}, \dots, e^{i\theta_d})
```

where each $\theta_k \in [0, 2\pi]$ is an independent phase. Integration over
this group corresponds to independent averaging over each phase:

```math
\int_{T^d} f(D) dD = \prod_{k=1}^d \left( \frac{1}{2\pi} \int_0^{2\pi} f(e^{i\theta_k}) d\theta_k \right)
```

## Integration Formula

For monomials in the matrix entries $D_{ij}$, the integral is non-zero only if:
1. All indices are diagonal ($i=j$ for all factors).
2. For each diagonal entry $D_{kk}$, the number of $D_{kk}$ factors matches the number of $\bar{D}_{kk}$ factors.

```math
\int_{T^d} D_{i_1 i_1} \dots D_{i_n i_n} \bar{D}_{j_1 j_1} \dots \bar{D}_{j_n j_n} dD = \begin{cases} 1 & \text{if } \{i_1, \dots, i_n\} = \{j_1, \dots, j_n\} \text{ as multisets} \\ 0 & \text{otherwise} \end{cases}
```

## Usage in IntU.jl

Use the `dDiagUnitary` measure to perform these integrations.

1. **Basic Integration using `@integrate`**

The `@integrate` macro automatically identifies `D` as the random diagonal unitary.

```julia
using IntU, Symbolics
@variables d
# E[|D_11|^2] = 1
@integrate abs(D[1, 1])^2 dDiagUnitary(d)
# Output: 1
```

2. **Manual Integration**

```julia
using IntU, Symbolics

@variables d
D = SymbolicMatrix(:D, :DiagUnitary)

# E[|D_11|^2] = 1
integrate(abs(D[1, 1])^2, dDiagUnitary(d))

# E[D_11 * D_22^*] = 0 (independent phases)
integrate(D[1, 1] * conj(D[2, 2]), dDiagUnitary(d))

# Non-diagonal entries are zero by definition
integrate(abs(D[1,2])^2, dDiagUnitary(d)) # Output: 0
```

## Performance Note

Integration over the diagonal group is significantly faster than over the full
unitary group $U(d)$, as it avoids the combinatorial complexity of Weingarten
functions and reduces to simple multiset comparisons.
