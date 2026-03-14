# Circular Ensembles

IntU.jl provides support for the three classical Circular Ensembles of Random Matrix Theory:
- **CUE (Circular Unitary Ensemble)**: Corresponds to the Haar measure on the Unitary group $U(d)$.
- **COE (Circular Orthogonal Ensemble)**: Ensemble of symmetric unitary matrices ($S^T = S$).
- **CSE (Circular Symplectic Ensemble)**: Ensemble of self-dual unitary matrices ($S^R = J S^T J^T = S$) of even dimension $d=2N$.

These ensembles are defined on the unitary group, unlike the Gaussian ensembles which are defined on the space of Hermitian matrices.

## Measures

```@docs
dCUE
dCOE
dCSE
```

## Examples

### COE (Circular Orthogonal Ensemble)

For the COE, the matrix $S$ is symmetric unitary. The diagonal entries have different statistical properties than off-diagonal entries.

1. **Basic Integration using `@integrate`**

The `@integrate` macro automatically identifies `S` as the random matrix when `dCOE` is used.

```julia
using IntU, Symbolics
@variables d
# COE diagonal 2nd moment E[|S_{1,1}|^2]
@integrate abs(S[1, 1])^2 dCOE(d)
# Output: 2 / (d + 1)

# COE off-diagonal 2nd moment E[|S_{1,2}|^2]
@integrate abs(S[1, 2])^2 dCOE(d)
# Output: 1 / (d + 1)

# COE diagonal 4th moment E[|S_{1,1}|^4]
@integrate abs(S[1, 1])^4 dCOE(d)
# Output: 6 / ((d + 1) * (d + 3))

# COE off-diagonal 4th moment E[|S_{1,2}|^4]
@integrate abs(S[1, 2])^4 dCOE(d)
# Output: 3 / ((d + 1) * (d + 3))

# COE correlation moment E[|S_{1,1}|^2 |S_{1,2}|^2]
@integrate abs(S[1, 1])^2 * abs(S[1, 2])^2 dCOE(d)
# Output: 2 / ((d + 1) * (d + 3))
```

2. **Manual Integration**

```julia
using IntU, Symbolics
@variables d
S = SymbolicMatrix(:S, :COE)
# COE moment E[|S_{1,1}|^2]
integrate(abs(S[1, 1])^2, dCOE(d))
# Output: 2 / (d + 1)
```

### CSE (Circular Symplectic Ensemble)

For the CSE, the matrix $S$ is defined on a space of dimension $2N$ and satisfies $S = J S^T J^T$.

1. **Basic Integration using `@integrate`**

```julia
using IntU, Symbolics
@variables d
# CSE diagonal 2nd moment E[|S_{1,1}|^2]
@integrate abs(S[1, 1])^2 dCSE(d)
# Output: 1 / (d - 1)

# CSE diagonal 4th moment E[|S_{1,1}|^4]
@integrate abs(S[1, 1])^4 dCSE(d)
# Output: 2 / (d * (d - 1))
```

2. **Manual Integration**

```julia
using IntU, Symbolics
@variables d
S = SymbolicMatrix(:S, :CSE)
# CSE moment E[|S_{1,1}|^2]
integrate(abs(S[1, 1])^2, dCSE(d))
# Output: 1 / (d - 1)
```

### CUE (Circular Unitary Ensemble)

The CUE is statistical identical to the standard Unitary Haar measure.

1. **Basic Integration using `@integrate`**

```julia
using IntU, Symbolics
@variables d
# CUE moment E[|U_{1,1}|^2]
@integrate abs(U[1, 1])^2 dCUE(d)
# Output: 1 / d
```

2. **Manual Integration**

```julia
using IntU, Symbolics
@variables d
U = SymbolicMatrix(:U, :U)
# CUE moment E[|U_{1,1}|^2]
res = integrate(abs(U[1, 1])^2, dCUE(d))
# Output: 1 / d
```

## Potential Pitfalls

> [!IMPORTANT]
> ### Symbolic (d) Pitfalls
> - **Small Dimensions**: For Haar-related measures (Unitary, Orthogonal, Circular), results are rational functions with poles at small $d$ (typically $d < n$ for degree $n$ moments).
> - **Removable Singularities**: Substituting numeric values can yield $0/0$ forms (e.g., at $d=1$ or $d=2$).
> - **Automatic Handling**: `IntU.jl`'s `evaluate` function automatically simplifies expressions to resolve removable singularities when a denominator evaluates to zero.

## See Also

- [Unitary Integration](unitary_integration.md) — CUE is the Haar measure on $U(d)$; Weingarten calculus background
- [Orthogonal & Symplectic](orthogonal_integration.md) — COE arises from $O(d)$, CSE from $Sp(d)$
- [Asymptotic Expansions](asymptotic.md) — large-$d$ limit of circular ensemble moments
- [Integral Library](integral_library.md) — pre-computed trace moments
