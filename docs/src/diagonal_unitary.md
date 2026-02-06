# Diagonal Unitary Matrices (Torus Group)

This section describes integration over the group of diagonal unitary matrices,
also known as the torus group $T^d$.

## Overview

A diagonal unitary matrix $V$ of dimension $d$ has the form:

```math
V = \text{diag}(e^{i\theta_1}, e^{i\theta_2}, \dots, e^{i\theta_d})
```

where each $\theta_k \in [0, 2\pi]$ is an independent phase. Integration over
this group corresponds to independent averaging over each phase:

```math
\int_{T^d} f(V) dV = \prod_{k=1}^d \left( \frac{1}{2\pi} \int_0^{2\pi} f(e^{i\theta_k}) d\theta_k \right)
```

## Integration Formula

For monomials in the matrix entries $V_{ij}$, the integral is non-zero only if:
1. All indices are diagonal ($i=j$ for all factors).
2. For each diagonal entry $V_{kk}$, the number of $V_{kk}$ factors matches the number of $\bar{V}_{kk}$ factors.

```math
\int_{T^d} V_{i_1 i_1} \dots V_{i_n i_n} \bar{V}_{j_1 j_1} \dots \bar{V}_{j_n j_n} dV = \begin{cases} 1 & \text{if } \{i_1, \dots, i_n\} = \{j_1, \dots, j_n\} \text{ as multisets} \\ 0 & \text{otherwise} \end{cases}
```

## Usage in IntU.jl

Use the `dDiagUnitary` measure to perform these integrations.

```julia
using IntU, Symbolics

@variables d
@variables V[1:d, 1:d]::Complex
measure = dDiagUnitary(V, d)

# E[|V_11|^2] = 1
integrate(abs(V[1,1])^2, measure)

# E[V_11 * V_22^*] = 0 (independent phases)
integrate(V[1,1] * conj(V[2,2]), measure)

# Non-diagonal entries are zero by definition
integrate(abs(V[1,2])^2, measure) # Output: 0
```

## Performance Note

Integration over the diagonal group is significantly faster than over the full
unitary group $U(d)$, as it avoids the combinatorial complexity of Weingarten
functions and reduces to simple multiset comparisons.
