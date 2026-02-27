# Quantum Information Helpers

A collection of utility functions to assist with calculations common in Quantum
Information.

## Functions

### `partial_trace(M, dims, subsystem)`
Symbolically computes the partial trace of a composite system.
- `M`: Matrix to trace.
- `dims`: Tuple of dimensions of subsystems, e.g., `(2, 2)` for two qubits.
- `subsystem`: Index of the subsystem to trace out (1 or 2, etc.).

## Example: Average Purity

Calculating the average purity of a subsystem when the global system is in a random pure state. For a bipartite system $d = d_A d_B$, the average purity is given by $\langle \gamma \rangle = \frac{d_A + d_B}{d + 1}$.

```julia
using IntU, Symbolics

# System Dimensions: 2 qubits
d_A, d_B = 2, 2
d = d_A * d_B

# Random Unitary on full system
U = SymbolicMatrix(:U, :U, d)

# Pure state |psi> = U |1> (first column of U)
psi = U[:, 1]
rho = psi * adjoint(psi)

# Reduced density matrix rho_A by tracing out subsystem 2
rho_A = partial_trace(rho, (d_A, d_B), 2)

# Calculate Haar-average purity manually
# purity(rho_A) = tr(rho_A^2)
avg_pury = integrate(tr(rho_A * rho_A), dU(d))

println(avg_pury)
# Output: 4//5 (matches (2+2)/(4+1))
```

## References

- Nielsen, M. A., & Chuang, I. L. (2010). *Quantum computation and quantum
  information*. Cambridge university press.
- Watrous, J. (2018). *The theory of quantum information*. Cambridge University
  Press.
