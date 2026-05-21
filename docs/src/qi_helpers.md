# Quantum Information Helpers

Utility functions for calculations common in quantum information theory.
The current helper is `partial_trace`; additional functions (fidelity,
concurrence, relative entropy) are planned for future releases.

## `partial_trace`

```julia
partial_trace(M, dims, subsystem)
```

Computes the partial trace of a bipartite (or multipartite) composite
system. Supports symbolic matrix entries with concrete integer dimensions.

- **`M`**: the density matrix (or arbitrary matrix) to trace over, of size
  $d_1 d_2 \times d_1 d_2$.
- **`dims`**: a tuple of subsystem dimensions, e.g. `(d_A, d_B)`.
- **`subsystem`**: the index of the subsystem to trace *out* (1-based).

The subsystem dimensions must be concrete integers. The matrix entries
may be symbolic (e.g., products of `SymbolicMatrix` elements), so the
result retains symbolic dependence for further integration.

## Example: Average Purity of a Random Bipartite State

For a random pure state $|\psi\rangle$ on $\mathbb{C}^{d_A} \otimes \mathbb{C}^{d_B}$,
Page's formula gives the average purity of subsystem $A$:

$$\langle \mathrm{tr}(\rho_A^2) \rangle = \frac{d_A + d_B}{d_A d_B + 1}.$$

The following reproduces this result from first principles:

```julia
using IntegrateUnitary, Symbolics

# System dimensions: two qubits
d_A, d_B = 2, 2
d = d_A * d_B    # total dimension = 4

# Random unitary on the full system
U = SymbolicMatrix(:U, :U, d)

# Pure state |psi> = U|0> (first column of U)
psi = U[:, 1]
rho = psi * adjoint(psi)

# Reduced density matrix: trace out subsystem B
rho_A = partial_trace(rho, (d_A, d_B), 2)

# Average purity tr(rho_A^2) under the Haar measure
avg_purity = integrate(tr(rho_A * rho_A), dU(d))
println(avg_purity)
# Output: 4//5  (matches (2+2)/(4+1))
```

The result $4/5$ agrees with the Page formula $(d_A + d_B)/(d_A d_B + 1) = 4/5$.

## Example: Asymptotic Purity

The Page formula can be combined with `asymptotic` to recover the
large-system behaviour directly:

```julia
using IntegrateUnitary, Symbolics
@variables n

# Page formula for equal-size subsystems (d_A = d_B = n)
page_purity = 2n / (n^2 + 1)
asymptotic(page_purity, n, 5)
# Output: 2/n - 2/n^3 + 2/n^5
# Leading term: 2/n (same O(1/n) scaling as maximally mixed purity 1/n)
```

> [!NOTE]
> For large subsystems, the purity is order $1/n$:
> $\mathbb{E}[\mathrm{tr}(\rho_A^2)] = 2/n + O(1/n^3)$ for equal bipartitions.
> A maximally mixed $n$-dimensional state has purity exactly $1/n$, so this
> result indicates highly mixed (and nearly maximally entangled) subsystems,
> not equality with the maximally mixed value.

## See Also

- [Pure States](pure_states.md) — integration over Haar-random pure states
- [Stiefel Manifolds](stiefel_manifold.md) — generalisation to $k$-frames
- [Asymptotic Expansions](asymptotic.md) — large-$d$ behaviour
- [Integral Library](integral_library.md) — pre-computed moments

## References

- Page, D. N. (1993). Average entropy of a subsystem. *Physical Review Letters*,
  71(9), 1291–1294.
- Nielsen, M. A., & Chuang, I. L. (2010). *Quantum computation and quantum
  information*. Cambridge University Press.
- Watrous, J. (2018). *The theory of quantum information*. Cambridge University
  Press.

See [`partial_trace`](@ref) in the [API Reference](api.md).
