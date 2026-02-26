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
# COE moment E[|S_{1,1}|^2]
@integrate abs(S[1, 1])^2 dCOE(d)
# Output: 2 / (d + 1)
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
# CSE moment E[|S_{1,1}|^2]
@integrate abs(S[1, 1])^2 dCSE(d)
# Output: 1 / (d - 1)
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
