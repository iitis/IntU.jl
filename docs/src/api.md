# API Reference

## Core Integration

```@docs
integrate
@integrate
evaluate
asymptotic
hciz
vandermonde_det
```

> [!NOTE]
> `integrate(expr, measure)` is the universal entry point for all calculations
> in IntegrateUnitary.jl. It handles matrix-valued expressions and library lookups, and
> supports symbolic dimensions for a broad class of entry-wise and
> trace-polynomial workflows. Some paths require concrete integer dimensions
> (including `|tr(U)|^(2k)` for `k > 1`, `hciz` on `SymbolicMatrix` inputs, and
> direct matrix-valued integration of `SymbolicMatrix` /
> `SymbolicMatrixProduct` expressions).

## Measures

```@docs
IntegrateUnitary.AbstractMeasure
IntegrateUnitary.measure_info
```

### Unitary Group
```@docs
dU
dSU
```

### Orthogonal & Symplectic
```@docs
dO
dSp
```

### Circular Ensembles

See [Circular Ensembles](circular_ensembles.md) for `dCOE`, `dCSE`, `dCUE`.

### Gaussian Ensembles
```@docs
dGUE
dGOE
dGSE
dGinUE
dGinOE
dGinSE
```

### Pure States
```@docs
dPsi
```

See [Stiefel Manifolds](stiefel_manifold.md) for `dStiefel`.

### Permutation Groups
```@docs
dPerm
dCPerm
```

### Diagonal Unitary Matrices
```@docs
dDiagUnitary
IntegrateUnitary.DiagonalUnitaryMeasure
```

## Symbolic Helpers

```@docs
IntegrateUnitary.SymbolicMatrix
tr
symbolic_unitary
symbolic_orthogonal
symbolic_symplectic
symbolic_pure_state
symbolic_permutation
```

## Quantum Information Utilities

```@docs
partial_trace
```

## Internal / Advanced

These functions are part of the internal machinery but documented for development reference.

### Integration Engine & Helpers

```@docs
integrate_indices
IntegrateUnitary.integrate_indices_symplectic
IntegrateUnitary.integrate_indices_permutation
IntegrateUnitary.integrate_indices_diagonal
IntegrateUnitary.integrate_indices_coe
IntegrateUnitary.integrate_indices_gue
IntegrateUnitary.integrate_indices_goe
IntegrateUnitary.integrate_indices_ginue
IntegrateUnitary.integrate_indices_ginoe
IntegrateUnitary.integrate_indices_ginse
IntegrateUnitary.tr_lazy
LazyTrace
LazySum
IntegrateUnitary.check_library
IntegrateUnitary.tr_val
IntegrateUnitary._expand_asymptotic
IntegrateUnitary._poly_degree
IntegrateUnitary._ensure_symbolic_dim
IntegrateUnitary._try_numeric
IntegrateUnitary._try_extract_int
IntegrateUnitary.robust_substitute
IntegrateUnitary.get_full_cycle_type
IntegrateUnitary.get_weingarten_reduced_data
```

### Matcher and Logic

```@docs
AbstractIndexMatcher
MetadataMatcher
_integrate_core
process_term
weingarten
IntegrateUnitary.ParityUnionFind
```

### Weingarten & Combinatorics

```@docs
IntegrateUnitary.weingarten_orthogonal_val
IntegrateUnitary.weingarten_symplectic_val
IntegrateUnitary.get_pair_partitions
IntegrateUnitary.get_matching_pair_partitions_filtered
IntegrateUnitary.canonicalize_pair_partition
IntegrateUnitary.conjugate_partition
IntegrateUnitary.murnaghan_nakayama
IntegrateUnitary.character_at_id
IntegrateUnitary.irrep_dimension
IntegrateUnitary.compute_symplectic_contraction
IntegrateUnitary.weingarten_orthogonal_val_canonical
IntegrateUnitary.INTEGRATION_RULES
```
