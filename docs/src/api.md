# API Reference

## Core Integration

```@docs
integrate
asymptotic
```

## Measures

### Unitary Group
```@docs
dU
symbolic_dimension_unitary
@symbolic_dimension
```

### Orthogonal & Symplectic
```@docs
dO
dSp
```

### Gaussian Ensembles
```@docs
dGUE
dGOE
dGSE
```

### Pure States
```@docs
dPsi
```

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
SymbolicMatrix
tr
```

## Quantum Information Utilities

```@docs
purity
average_purity
fidelity
average_fidelity
partial_trace
```

## Internal / Advanced

These functions are part of the internal machinery but documented for development reference.

```@docs
integrate_indices
tr_lazy
LazyTrace
LazySum
IntU.integrate_indices_gue
IntU.integrate_indices_goe
IntU.integrate_indices_orthogonal
IntU.integrate_indices_symplectic
IntU.integrate_indices_coe
IntU.integrate_indices_cse
IntU.integrate_indices_permutation
IntU.integrate_indices_diagonal
IntU.fallback_integrate
IntU.check_library
IntU.tr_val
IntU._expand_asymptotic
IntU._poly_degree
```

### Matcher and Logic

```@docs
AbstractIndexMatcher
LookupMatcher
SymbolicUnitary
_integrate_core
process_term
weingarten
```

### Weingarten & Combinatorics

```@docs
IntU.weingarten_orthogonal_val
IntU.weingarten_symplectic_val
IntU.get_pair_partitions
IntU.get_matching_pair_partitions_filtered
IntU.canonicalize_pair_partition
IntU.conjugate_partition
IntU.count_loops
IntU.murnaghan_nakayama
IntU.character_at_id
IntU.irrep_dimension
IntU.get_weingarten_orthogonal_data
IntU.compute_symplectic_contraction
```
