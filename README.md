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
- **Asymptotic Expansions**: Compute Taylor series expansions of integrals in powers of $1/d$.
- **Quantum Information Helpers**: Built-in functions for calculating average purity,
  fidelity, and partial traces of symbolic densitiy matrices.
- **Automated Weingarten Calculus**: Handles the combinatorial complexity of
  Weingarten functions (`Wg`) automatically.
- **Symbolics.jl Integration**: Built on top of `Symbolics.jl` for powerful
  expression manipulation.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/iitis/IntU.jl")
```

## Usage

### 1. Basic Integration over $U(d)$

Calculate moments of unitary entries, e.g., $\int |U_{11}|^2 dU$.

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

### 2. Quantum Information Helpers

Compute average purity, fidelity, and symbolic partial traces.

```julia
# Average Purity of a Haar-randomly rotated state
rho_fixed = [1.0 0.0; 0.0 0.0]
rho_random = U * rho_fixed * U'
avg_pur = average_purity(rho_random, measure)

# Symbolic Partial Trace
M = [0.5 0 0 0.5; 0 0 0 0; 0 0 0 0; 0.5 0 0 0.5] # Bell state
rho_a = partial_trace(M, (2, 2), 2) # Tracing out subsystem 2
```

### 3. Integration over Pure States

Integrate over random pure states $|\psi\rangle \in \mathbb{C}^d$.

```julia
@variables psi[1:d]::Complex
measure_psi = dPsi(psi, d)

# Average fidelity with a fixed state phi
@variables phi[1:d]::Complex
expr = abs(psi' * phi)^2
result = integrate(expr, measure_psi)
# Output simplified: (sum(|phi_i|^2)) / d
```

### 4. Large-$d$ Asymptotic Expansions

Compute series expansions of integrals around $d \to \infty$.

```julia
# Expand |U_{1,1}|^4 up to order 1/d^4
expr = abs(U[1,1])^4
res = asymptotic(expr, measure, 4)
println(res)
# Output: 2/d^2 - 2/d^3 + 2/d^4
```

## Development and Verification

The package follows a modular architecture for maintainability. You can run all
tests, benchmarks, and examples using the provided scripts:

- **Tests**: `julia --project=. test/runtests.jl` (verbose output enabled).
- **Examples**: `./examples/runexamples.sh`
- **Benchmarks**: `./benchmarks/runbenchmarks.sh`

## How It Works

`IntU.jl` parses input expressions to identify indices of unitary variables,
applying the **Weingarten formula**:

$$ \int_{U(d)} U_{i_1 j_1} \dots U_{i_n j_n} \bar{U}_{k_1 l_1} \dots
\bar{U}_{k_n l_n} dU = \sum_{\sigma, \tau \in S_n} \delta_{i, k_\sigma}
\delta_{j, l_\tau} \text{Wg}(\sigma \tau^{-1}, d) $$

## License

Apache 2.0
