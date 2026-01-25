# IntU.jl

**IntU.jl** is a Julia package for performing **symbolic integration over the
Haar measure** of Unitary ($U(d)$) groups and pure quantum states. It leverages
**Weingarten Calculus** to compute integrals of polynomial functions of matrix
elements exactly, supporting both concrete and **symbolic dimensions** ($d$).

## Features

- **Symbolic Integration**: Integrate polynomials of unitary matrix elements
  $U_{ij}$ and $\bar{U}_{kl}$.
- **Symbolic Dimension**: Results can depend on a symbolic variable $d$,
  allowing for large-$d$ analysis.
- **Pure States**: Integration over Haar-random pure states $|\psi\rangle$
  (equivalent to the first column of a random unitary).
- **Asymptotic Expansions**: Compute Taylor series expansions of integrals in
  powers of $1/d$.
- **Quantum Information Helpers**: Built-in functions for calculating average
  purity, fidelity, and partial traces of symbolic densitiy matrices.
- **Automated Weingarten Calculus**: Handles the combinatorial complexity of
  Weingarten functions (`Wg`) automatically.
- **Symbolics.jl Integration**: Built on top of `Symbolics.jl` for powerful
  expression manipulation.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/iitis/IntU.jl")
```

## Quick Start

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
println(result) # Output: 1/d
```

Check out the [Examples](@ref) section for more advanced usage.
