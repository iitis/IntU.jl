# IntU.jl Documentation

Welcome to the documentation for **IntU.jl**, a generic symbolic integration engine for random matrices and quantum states.

## Getting Started

IntU.jl is designed to make "back-of-the-envelope" calculations in Random Matrix Theory and Quantum Information exact and automated.

### Installation

```julia
using Pkg
Pkg.add(url="https://github.com/iitis/IntU.jl")
```

### Reproducible Setup From a Checkout

For paper-grade reproducibility, use a tagged checkout (for this manuscript:
`v1.0.0`) and instantiate the pinned project environments used by the
examples and benchmarks:

```julia
using Pkg
Pkg.activate("examples"); Pkg.instantiate()
Pkg.activate("benchmarks"); Pkg.instantiate()
```

This uses the versioned `Manifest.toml` files in `examples/` and `benchmarks/`.

## The @integrate Macro

For more intuitive symbolic integration, **IntU.jl** provides the `@integrate` macro. It automatically identifies random matrices based on the measure and manages variable declarations.

```julia
@integrate expr measure
```

> [!TIP]
> **Symbol Scope and Redefinition**: The macro manages a persistent symbolic state. If a symbol is used in one context (e.g., as a random matrix for `dU`) and then in another (e.g., as a constant for `dO`), the macro automatically re-binds it to the correct type. This "Safety Rebind" prevents silent mathematical errors during sequential execution.

## Common Interface

The primary way to interact with **IntU.jl** is through the `integrate` function. It provides a unified interface for all supported groups and ensembles.

```julia
integrate(expr, measure)
```

- **`expr`**: A symbolic expression (e.g., product of `SymbolicMatrix` elements) or an array of such expressions.
- **`measure`**: A measure object defining the group/ensemble and its dimension (e.g., `dU(d)`, `dO(d)`, `dGUE(d)`).

## Choosing a Measure

| If you want to average over… | Measure | Page |
|---|---|---|
| Generic random unitary matrices | `dU(d)` | [Unitary Integration](unitary_integration.md) |
| Special unitary (balanced polynomials) | `dSU(d)` | [Unitary Integration](unitary_integration.md) |
| Real orthogonal matrices | `dO(d)` | [Orthogonal & Symplectic](orthogonal_integration.md) |
| Symplectic matrices ($d$ even) | `dSp(d)` | [Orthogonal & Symplectic](orthogonal_integration.md) |
| Haar-random unitary matrices (CUE) | `dCUE(d)` | [Circular Ensembles](circular_ensembles.md) |
| Symmetric unitary matrices (COE) | `dCOE(d)` | [Circular Ensembles](circular_ensembles.md) |
| Self-dual unitary matrices (CSE) | `dCSE(d)` | [Circular Ensembles](circular_ensembles.md) |
| Hermitian random matrices | `dGUE(d)` / `dGOE(d)` / `dGSE(d)` | [Gaussian Ensembles](gaussian_integration.md) |
| Non-Hermitian random matrices | `dGinUE(d)` / `dGinOE(d)` / `dGinSE(d)` | [Gaussian Ensembles](gaussian_integration.md) |
| Random permutation matrices | `dPerm(d)` | [Permutation Groups](permutation_integration.md) |
| Centered permutation matrices | `dCPerm(d)` | [Permutation Groups](permutation_integration.md) |
| Random pure states | `dPsi(d)` | [Pure States](pure_states.md) |
| Orthonormal $k$-frames (Stiefel) | `dStiefel(d, k)` | [Stiefel Manifolds](stiefel_manifold.md) |
| Independent diagonal phases | `dDiagUnitary(d)` | [Diagonal Unitaries](diagonal_unitary.md) |
| Moments up to order $t$ only | `dDesign(d, t)` | [Unitary Integration](unitary_integration.md) |

## Manual

- [Unitary Integration](unitary_integration.md): Core functionality for $U(d)$ and $SU(d)$ (balanced polynomials).
- [Diagonal Unitaries](diagonal_unitary.md): Integration over the Torus group.
- [Orthogonal & Symplectic](orthogonal_integration.md): Integration over real groups.
- [Gaussian Ensembles](gaussian_integration.md): GUE, GOE, and GSE.
- [Circular Ensembles](circular_ensembles.md): COE, CUE, and CSE.
- [Permutation Groups](permutation_integration.md): Symmetric and centered permutations.
- [Pure States](pure_states.md): Integration over random vectors.
- [Stiefel Manifolds](stiefel_manifold.md): Integration over Stiefel manifolds $V_k(\mathbb{C}^d)$.
- [Symbolic Trace Logic](symbolic_trace.md): Index-free matrix integration.
- [Asymptotic Expansions](asymptotic.md): Large-$d$ limit analysis.
- [Integral Library](integral_library.md): Pre-computed standard results.
- [QI Helpers](qi_helpers.md): Partial trace and quantum information helpers.
- [ITensors Integration](itensors.md): Symmetric integration of tensor networks.

## API Reference

See the [API Reference](api.md) for detailed function signatures.
