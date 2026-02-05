# Permutation Group Integration

IntU.jl supports symbolic integration over the Symmetric Group $S_d$ (the group of $d \times d$ permutation matrices) and the ensemble of centered permutation matrices.

## Permutation Matrices ($S_d$)

Integration over the Haar measure of the permutation group computes the average of products of matrix entries $P_{i,j}$. Since each entry $P_{ij}$ is either 0 or 1, and there is exactly one '1' in each row and column, the integral is non-zero only if the indices are consistent with a permutation.

### Usage

Use `dPerm(d)` or `dPerm(P, d)` to define the measure.

```julia
using IntU, Symbolics

@variables d
@variables P[1:2, 1:2]
measure = dPerm(P, d)

# Average of a single entry
integrate(P[1,1], measure)
# Output: 1 / d

# Average of a product of entries
integrate(P[1,1] * P[2,2], measure)
# Output: 1 / (d * (d - 1))
```

### Integration Rule

The integral of a monomial $P_{i_1, j_1} P_{i_2, j_2} \dots P_{i_k, j_k}$ is:
- $0$ if any two indices $i_m, i_n$ are equal while $j_m \neq j_n$ (or vice versa).
- $\frac{(d-k)!}{d!}$ if all row indices and all column indices are distinct among the $k$ unique pairs $(i_m, j_m)$.

## Centered Permutation Matrices

Centered permutation matrices $Y$ are defined as $Y = P - J/d$, where $J$ is the all-ones matrix. These matrices satisfy $\sum_i Y_{ij} = \sum_j Y_{ij} = 0$.

### Usage

Use `dCPerm(d)` or `dCPerm(Y, d)` to define the measure.

```julia
@variables Y[1:2, 1:2]
m_centered = dCPerm(Y, d)

# First moment is zero by centering
integrate(Y[1,1], m_centered)
# Output: 0

# Second moment (variance-like)
integrate(Y[1,1]^2, m_centered)
# Output: (d - 1) / d^2
```

### Implementation Detail

Integration for centered permutations is handled by substituting $Y_{ij} = P_{ij} - 1/d$ and expanding the resulting polynomial, which is then integrated using the permutation group rules.
