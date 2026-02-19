# Permutation Group Integration

IntU.jl provides support for symbolic integration over two related ensembles: the **Symmetric Group** $S_d$ and the **Centered Permutation Ensemble**.

## Ordinary Permutation Matrices ($S_d$)

The **Symmetric Group** $S_d$ consists of all $d!$ possible $d \times d$
permutation matrices. A permutation matrix $P$ has exactly one entry equal to
$1$ in each row and column, with all other entries being $0$.

Integration over $S_d$ in IntU.jl computes the uniform average over this
discrete set:
```math
\mathbb{E}_{P \in S_d}[f(P)] = \frac{1}{d!} \sum_{P \in S_d} f(P)
```

### Usage

Use `dPerm(d)` or `dPerm(P, d)` to define the measure.

```julia
using IntU, Symbolics

@variables d
P = SymbolicMatrix(:P, :Perm)

# Expected value of a single entry: E[P_ij] = 1/d
integrate(P[1,1], measure)
# Output: 1 / d

# Expected value of a product: E[P_11 * P_22] = 1 / (d(d-1))
integrate(P[1, 1] * P[2, 2], dPerm(d))
# Output: 1 / (d * (d - 1))
```

### Integration Rule

The integral of a monomial $P_{i_1, j_1} P_{i_2, j_2} \dots P_{i_k, j_k}$
follows simple combinatorial rules:
1. **Consistency**: If the set of pairs $\{(i_m, j_m)\}$ overlaps in rows but
   not columns (e.g., $P_{11}P_{12}$) or vice versa, the integral is **$0$**
   because no permutation matrix can have two $1$s in the same row or column.
   $i_m$ are distinct and all $j_m$ are distinct), the result is:
```math
\mathbb{E}[P_{i_1, j_1} \dots P_{i_k, j_k}] = \frac{(d-k)!}{d!}
```

---

## Centered Permutation Matrices

The **Centered Permutation Ensemble** consists of matrices $Y$ obtained by
subtracting the "flat" matrix $J/d$ (where $J$ is the all-ones matrix) from a
permutation matrix $P \in S_d$:
```math
Y = P - \frac{1}{d} J, \quad Y_{ij} = P_{ij} - \frac{1}{d}
```

These matrices are "centered" because they satisfy:
```math
\sum_{i=1}^d Y_{ij} = 0, \quad \sum_{j=1}^d Y_{ij} = 0
```
This ensemble is particularly useful in studying
fluctuations and correlations in permutations.

### Usage

Use `dCPerm(d)` or `dCPerm(Y, d)` to define the measure.

```julia
Y = SymbolicMatrix(:Y, :CPerm)

# The first moment is zero by definition
integrate(Y[1,1], m_centered)
# Output: 0

# The second moment (variance) is E[(P_11 - 1/d)^2] = 1/d - 1/d^2 = (d-1)/d^2
integrate(Y[1, 1]^2, dCPerm(d))
# Output: (d - 1) / d^2
```

### Implementation Detail

Integration for centered permutations is handled by substituting $Y_{ij} = P_{ij} - 1/d$ and expanding the resulting polynomial, which is then integrated using the permutation group rules.

## Symbolic Traces

IntU.jl supports integration of traces involving permutation matrices. While these can be integrated using `Symbolics.jl` arrays and `scalarize`, the combinatorial nature of $S_d$ means that correlations are correctly handled even for large expressions.

```julia
@variables A[1:d, 1:d]
# E[tr(P * A)]
expr = Symbolics.scalarize(IntU.tr(P * A))
integrate(expr, measure)
# Output: Sum(A_ij) / d
```
