# IntU.jl

**IntU.jl** is a Julia package for performing **symbolic integration over the
Haar measure** of Unitary ($U(d)$) groups and pure quantum states. It leverages
**Weingarten Calculus** to compute integrals of polynomial functions of matrix
elements exacty, supporting both concrete and **symbolic dimensions** ($d$).

## Features

- **Symbolic Integration**: Integrate polynomials of unitary matrix elements
  $U_{ij}$ and $\bar{U}_{kl}$.
- **Symbolic Dimension**: Results can simply depend on a symbolic variable $d$,
  allowing for large-$d$ analysis.
- **Pure States**: Integration over Haar-random pure states $|\psi\rangle$
  (equivalent to the first column of a random unitary).
- **Automated Weingarten Calculus**: Handles the combinatorial complexity of
  Weingarten functions (`Wg`) automatically.
- **Symbolics.jl Integration**: Built on top of `Symbolics.jl` for powerful
  expression manipulation.

## Installation

```julia
using Pkg
Pkg.add("IntU") # Pending registration
# Or from source:
Pkg.add(url="https://github.com/iitis/IntU.jl")
```

## Usage

### 1. Basic Integration over $U(d)$

Calculate moments of unitary entries, e.g., $\int |U_{11}|^2 dU$ or traces.

```julia
using IntU
using Symbolics

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

Compute integrals involving traces, such as $\int |\text{Tr}(U)|^4 dU$.

```julia
expr = abs(tr(U))^4
result = integrate(expr, measure)
# Result: 2 (for d >= 2)
```

### 3. Integration over Pure States

Integrate over random pure states $|\psi\rangle \in \mathbb{C}^d$.

```julia
@variables psi[1:d]::Complex
measure_psi = dPsi(psi, d)

# Average purity of a subsystem?
# Expectation of |<psi|phi>|^2 (fidelity)
@variables phi[1:d]::Complex
expr = abs(psi' * phi)^2
result = integrate(expr, measure_psi)
# Output simplified: (sum(|phi_i|^2)) / d
```

### 4. Determinants and Minors

While `det(U)` expands symbolically, `IntU` can integrate minors.

```julia
# Explicitly collect into matrix for determinant expansion
U_mat = collect(U) 
expr = abs(det(U_mat))^2 
# ... integration logic
```

## How It Works

`IntU.jl` parses input expressions to identify indices of unitary variables. It
matches $U_{ij} \dots \bar{U}_{kl} \dots$ terms and applies the **Weingarten
formula**:

$$ \int_{U(d)} U_{i_1 j_1} \dots U_{i_n j_n} \bar{U}_{k_1 l_1} \dots
\bar{U}_{k_n l_n} dU = \sum_{\sigma, \tau \in S_n} \delta_{i, k_\sigma}
\delta_{j, l_\tau} \text{Wg}(\sigma \tau^{-1}, d) $$

Where $\text{Wg}$ is the Weingarten function for the unitary group.

## License

Apache 2.0
