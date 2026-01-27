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

```julia
using IntU, Symbolics
@variables U[1:2, 1:2]::Complex
measure = dU(U, 2)

rho_in = [1 0; 0 0] # Pure state |0><0|
rho_out = U * rho_in * U'

# Partial trace over subsystem B of a 2-qubit system?
# (Example requires 4x4 U for 2-qubit)
```

## References

- Nielsen, M. A., & Chuang, I. L. (2010). *Quantum computation and quantum
  information*. Cambridge university press.
- Watrous, J. (2018). *The theory of quantum information*. Cambridge University
  Press.
