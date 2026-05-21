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

Use `dPerm(d)` to define the measure.

1. **Basic Integration using `@integrate`**

The `@integrate` macro automatically identifies `P` as the random permutation matrix.

```julia
using IntU, Symbolics
@variables d
# E[P_11] = 1/d
@integrate P[1,1] dPerm(d)
# Output: 1 / d
```

2. **Manual Integration**

```julia
using IntU, Symbolics

@variables d
P = SymbolicMatrix(:P, :Perm)

# Expected value of a single entry: E[P_ij] = 1/d
integrate(P[1,1], dPerm(d))
# Output: 1 / d

# Expected value of a product: E[P_11 * P_22] = 1 / (d(d-1))
integrate(P[1, 1] * P[2, 2], dPerm(d))
# Output: 1 / (d * (d - 1))
```

### Integration Rule

The integral of a monomial $P_{i_1, j_1} P_{i_2, j_2} \dots P_{i_k, j_k}$
follows simple combinatorial rules:
1. **Deduplication**: Repeated identical factors are collapsed first, i.e., the
   rule is applied to the set of distinct pairs $\{(i_m, j_m)\}$.
2. **Consistency**: If this set overlaps in rows but
   not columns (e.g., $P_{11}P_{12}$) or vice versa, the integral is **$0$**
   because no permutation matrix can have two $1$s in the same row or column.
3. **Monomial Integration**: If the distinct pairs are consistent (all rows
   distinct and all columns distinct), the result is:
```math
\mathbb{E}[P_{i_1, j_1} \dots P_{i_k, j_k}] = \frac{(d-k)!}{d!}
```
where $k$ is the number of distinct index pairs after deduplication.

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

Use `dCPerm(d)` to define the measure.

1. **Basic Integration using `@integrate`**

The `@integrate` macro automatically identifies `Y` as the centered permutation matrix.

```julia
using IntU, Symbolics
@variables d
# Variance E[Y_11^2] = (d-1)/d^2
@integrate Y[1, 1]^2 dCPerm(d)
# Output: (d - 1) / d^2
```

2. **Manual Integration**

```julia
using IntU, Symbolics
@variables d
Y = SymbolicMatrix(:Y, :CPerm)

# The first moment is zero by definition
integrate(Y[1,1], dCPerm(d))
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
using IntU, Symbolics
@variables d
# E[tr(P * A)]
@integrate tr(P * A) dPerm(d)
# Output: Sum(A_ij) / d
```

## See Also

- [Unitary Integration](unitary_integration.md) — Weingarten calculus for $U(d)$
- [Diagonal Unitaries](diagonal_unitary.md) — integration over the torus (phase matrices)
- [Asymptotic Expansions](asymptotic.md) — large-$d$ limits
- [Symbolic Trace Logic](symbolic_trace.md) — index-free trace expressions
