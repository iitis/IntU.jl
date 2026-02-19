# Special Unitary Group Integration

IntU.jl supports integration over the Special Unitary Group $SU(d)$.

## Overview

The Special Unitary group $SU(d)$ consists of $d \times d$ unitary matrices with determinant 1.
In the regime where $d$ is large (stable range) or symbolic, the integration over $SU(d)$ for polynomial functions matches the integration over $U(d)$ for "balanced" polynomials (where the number of $U$ and $\bar{U}$ indices is equal). Unbalanced polynomials vanish in this regime (though for specific finite $d$, terms involving $\epsilon$-tensors might survive, which are not currently covered involving explicit non-stable range logic).

## Usage

Use the `dSU` measure constructor.

```julia
using IntU, Symbolics
@variables d
U = SymbolicMatrix(:U, :U)
# E[|U_{1,1}|^2]
res = integrate(abs(U[1, 1])^2, dSU(d))
# Output: 1/d

# Unbalanced moment
integrate(U[1, 1], dSU(d))
# Output: 0
```

## Functions

```@docs
dSU
SpecialUnitary
```
