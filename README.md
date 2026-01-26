# IntU.jl

**IntU.jl** is a Julia package for performing **symbolic integration over the
Haar measure** of Unitary ($U(d)$), Orthogonal ($O(d)$), and Symplectic ($Sp(d)$) groups, as well as random pure quantum states. It leverages **Weingarten Calculus** to compute integrals of polynomial functions of matrix elements exactly, supporting both concrete and **symbolic dimensions** ($d$).

## Features

- **Symbolic Integration**: Integrate polynomials of matrix elements $U_{ij}$, $O_{ij}$, $S_{ij}$.
- **Multiple Groups**: Support for Unitary $U(d)$, Orthogonal $O(d)$, and Symplectic $Sp(d)$.
- **Symbolic Dimension**: Results can depend on a symbolic variable $d$, allowing for large-$d$ analysis.
- **Pure States**: Integration over Haar-random pure states $|\psi\rangle$.
- **Asymptotic Expansions**: Compute Taylor series expansions of integrals in powers of $1/d$.
- **Quantum Information Helpers**: Built-in functions for calculating average purity, fidelity, and partial traces.
- **Symbolic Trace Logic**: Index-free integration of traces like $\operatorname{tr}(U A U^\dagger B)$.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/iitis/IntU.jl")
```

## Usage

### 1. Basic Integration over $U(d)$

Calculate moments of unitary entries, e.g., $\int |U_{11}|^2 dU$.

```julia
using IntU, Symbolics

@variables d
@variables U[1:d, 1:d]::Complex

# Define the Haar measure for U(d)
measure = dU(U, d)

# Define an integrand: |U_{1,1}|^2
expr = abs(U[1,1])^2

# Integrate
result = integrate(expr, measure)
println(result)
# Output: 1/d
```

### 2. Symbolic Traces

Perform symbolic integration of traces without explicit indices.

```julia
# Define symbolic matrices
U_sym = SymbolicMatrix(:U, false, :U)
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)

# Compute Integral of Tr(U A U' B)
# Note: Use IntU.tr to avoid conflict with LinearAlgebra
expr = IntU.tr(U_sym * A * U_sym' * B)
res = integrate(expr, measure)
println(res)
# Output simplified: tr(A) * tr(B) / d
```

### 3. Orthogonal and Symplectic

```julia
# Orthogonal O(d)
@variables O_mat[1:d, 1:d]::Real
m_O = dO(O_mat, d)
integrate(O_mat[1,1]^2, m_O) # -> 1/d

# Symplectic Sp(d) (d must be even)
@variables S_mat[1:d, 1:d]::Complex
m_Sp = dSp(S_mat, d)
integrate(abs(S_mat[1,1])^2, m_Sp) # -> 1/d
```

## Documentation

Full documentation is available in the `docs/` directory.

## License

Apache 2.0
