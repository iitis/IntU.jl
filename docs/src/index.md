# IntU.jl Documentation

Welcome to the documentation for **IntU.jl**, a generic symbolic integration engine for random matrices and quantum states.

## Getting Started

INTU.jl is designed to make "back-of-the-envelope" calculations in Random Matrix Theory and Quantum Information exact and automated.

### Installation

```julia
using Pkg
Pkg.add(url="https://github.com/iitis/IntU.jl")
```

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

## Manual

- [Unitary Integration](unitary_integration.md): Core functionality for $U(d)$ and $SU(d)$.
- [Diagonal Unitaries](diagonal_unitary.md): Integration over the Torus group.
- [Orthogonal & Symplectic](orthogonal_integration.md): Integration over real groups.
- [Gaussian Ensembles](gaussian_integration.md): GUE, GOE, and GSE.
- [Circular Ensembles](circular_ensembles.md): COE, CUE, and CSE.
- [Permutation Groups](permutation_integration.md): Symmetric and centered permutations.
- [Pure States](pure_states.md): Integration over random vectors.
- [Symbolic Trace Logic](symbolic_trace.md): Index-free matrix integration.
- [Asymptotic Expansions](asymptotic.md): Large-$d$ limit analysis.
- [Integral Library](integral_library.md): Pre-computed standard results.
- [QI Helpers](qi_helpers.md): Tools for Purity, Fidelity, etc.
- [ITensors Integration](itensors.md): Symmetric integration of tensor networks.

## API Reference

See the [API Reference](api.md) for detailed function signatures.
