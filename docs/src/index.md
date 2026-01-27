# IntU.jl Documentation

Welcome to the documentation for **IntU.jl**, a generic symbolic integration engine for random matrices and quantum states.

## Getting Started

INTU.jl is designed to make "back-of-the-envelope" calculations in Random Matrix Theory and Quantum Information exact and automated.

### Installation

```julia
using Pkg
Pkg.add(url="https://github.com/iitis/IntU.jl")
```

## Manual

- [Unitary Integration](unitary_integration.md): Core functionality for $U(d)$.
- [Orthogonal & Symplectic](orthogonal_integration.md): Integration over real groups.
- [Gaussian Ensembles](gaussian_integration.md): GUE, GOE, and GSE.
- [Pure States](pure_states.md): Integration over random vectors.
- [Symbolic Trace Logic](symbolic_trace.md): Index-free matrix integration.
- [Asymptotic Expansions](asymptotic.md): Large-$d$ limit analysis.
- [Integral Library](integral_library.md): Pre-computed standard results.
- [QI Helpers](qi_helpers.md): Tools for Purity, Fidelity, etc.

## API Reference

See the [API Reference](api.md) for detailed function signatures.
