# Examples

## Basic Integration over $U(d)$

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
# Output: 1/d
```

## Integration over Pure States

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

## Large-$d$ Asymptotic Expansions

Compute series expansions of integrals around $d \to \infty$.

```julia
# Expand |U_{1,1}|^4 up to order 1/d^4
expr = abs(U[1,1])^4
res = asymptotic(expr, measure, 4)
# Output: 2/d^2 - 2/d^3 + 2/d^4
```

## Symbolic Trace Logic

Perform symbolic integration of traces of matrix products involving unitaries and constant matrices, without explicit indices.

```julia
# Define symbolic matrices
U_sym = SymbolicMatrix(:U, false, :U)
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)

# Compute Integral of Tr(U A U' B)
expr = tr_lazy(U_sym * A * U_sym' * B)
res = integrate(expr, measure)
# Output simplified: tr(A) * tr(B) / d
```

## Pure State Asymptotics

Asymptotic expansions follow the same syntax as the Haar measure.

```julia
# Expand average fideltity up to order 1/d^2
expr = abs(psi' * phi)^2
res = asymptotic(expr, measure_psi, 2)
# Output: 1/d
```

## Quantum Information Tasks

### Average Purity

Compute the average purity of a state evolved under a random unitary.

```julia
rho_fixed = [1.0 0.0; 0.0 0.0]
rho_random = U * rho_fixed * U'
avg_pur = average_purity(rho_random, measure)
```

### Partial Trace

Symbolic partial trace of a multipartite system.

```julia
M = [0.5 0 0 0.5; 0 0 0 0; 0 0 0 0; 0.5 0 0 0.5] # Bell state
rho_a = partial_trace(M, (2, 2), 2) # Tracing out subsystem 2
```
