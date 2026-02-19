# Quantum Information Helpers

A collection of utility functions to assist with calculations common in Quantum
Information.

## Functions

### `purity(rho)`
Calculates the purity $\gamma = \mathrm{tr}(\rho^2)$.
- For pure states, $\gamma = 1$.
- For maximally mixed states, $\gamma = 1/d$.

### `fidelity(rho, sigma)`
Calculates the fidelity between two states.
- Note: This implements the "Uhlmann fidelity" squared form for states commonly
  used in some contexts, or simply overlap $\mathrm{tr}(\rho \sigma)$.
  **Check implementation**: Currently defined as `tr(rho * sigma)`. For pure
  states $|\psi\rangle, |\phi\rangle$, this equals
  $|\langle \psi | \phi \rangle|^2$.

### `partial_trace(M, dims, subsystem)`
Symbolically computes the partial trace of a composite system.
- `M`: Matrix to trace.
- `dims`: Tuple of dimensions of subsystems, e.g., `(2, 2)` for two qubits.
- `subsystem`: Index of the subsystem to trace out (1 or 2, etc.).

## Example: Average Purity

Calculating the average purity of a subsystem when the global system is in a random pure state.

```julia
using IntU, Symbolics

# System Dimensions: d_A = 2, d_B = 2
d_A = 2
d_B = 2
d = d_A * d_B

# Random Unitary on full system
U = SymbolicMatrix(:U, :U, d)
measure = dU(d)

# Pure state |psi> = U |00> (first column of U)
# We form the density matrix rho = |psi><psi|
# rho_{ij} = U_{i,1} * conj(U_{j,1})

# We can conceptually construct the partial trace.
# But IntU provides helpers if you work with explicit indices.
# Or we can compute the integral of Purity(rho_A).

# < Purity(rho_A) > = (d_A + d_B) / (d + 1)
# For d_A=d_B=2, d=4 -> (2+2)/5 = 0.8

val = average_purity((d_A, d_B), 1) 
# Note: average_purity is a helper that returns the analytical formula or value
println(val)
# Output: 4//5 (for d_A=2, d_B=2)
```

## References

- Nielsen, M. A., & Chuang, I. L. (2010). *Quantum computation and quantum
  information*. Cambridge university press.
- Watrous, J. (2018). *The theory of quantum information*. Cambridge University
  Press.
