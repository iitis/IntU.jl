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
> `integrate(expr, measure)` is the universal entry point for all calculations in IntU.jl. It automatically handles symbolic dimensions, matrix-valued expressions, and library lookups.

## Measures

```@docs
IntU.AbstractMeasure
IntU.measure_info
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
IntU.DiagonalUnitaryMeasure
```

## Symbolic Helpers

```@docs
IntU.SymbolicMatrix
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
IntU.integrate_indices_symplectic
IntU.integrate_indices_permutation
IntU.integrate_indices_diagonal
IntU.integrate_indices_coe
IntU.integrate_indices_gue
IntU.integrate_indices_goe
IntU.integrate_indices_ginue
IntU.integrate_indices_ginoe
IntU.integrate_indices_ginse
IntU.tr_lazy
LazyTrace
LazySum
IntU.check_library
IntU.tr_val
IntU._expand_asymptotic
IntU._poly_degree
IntU._ensure_symbolic_dim
IntU._try_numeric
IntU._try_extract_int
IntU.robust_substitute
IntU.get_full_cycle_type
IntU.get_weingarten_reduced_data
```

### Matcher and Logic

```@docs
AbstractIndexMatcher
MetadataMatcher
_integrate_core
process_term
weingarten
IntU.ParityUnionFind
```

### Weingarten & Combinatorics

```@docs
IntU.weingarten_orthogonal_val
IntU.weingarten_symplectic_val
IntU.get_pair_partitions
IntU.get_matching_pair_partitions_filtered
IntU.canonicalize_pair_partition
IntU.conjugate_partition
IntU.murnaghan_nakayama
IntU.character_at_id
IntU.irrep_dimension
IntU.compute_symplectic_contraction
IntU.weingarten_orthogonal_val_canonical
IntU.INTEGRATION_RULES
```
