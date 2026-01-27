# API Reference

## Measures

```@docs
IntU.dU
IntU.dPsi
IntU.dGUE
IntU.dGOE
IntU.dGSE
```

## Integration Functions

```@docs
integrate
asymptotic
```

## Quantum Information Helpers

```@docs
purity
average_purity
fidelity
average_fidelity
partial_trace
```

## Symbolic Trace Logic

```@docs
SymbolicMatrix
tr_lazy
LazyTrace
```

## Internal Utilities

```@docs
integrate_indices
tr
```

## Real Groups

```@docs
dO
dSp
```

## Internal & Weingarten Utilities

These functions are part of the internal machinery but documented for development reference.

```@docs
IntU.weingarten_orthogonal_matrix
IntU.character_at_id
IntU._poly_degree
IntU.weingarten_symplectic_val
IntU.conjugate_partition
IntU.get_matching_pair_partitions_filtered
IntU.irrep_dimension
IntU.compute_symplectic_contraction
IntU.integrate_indices_symplectic
IntU.weingarten_orthogonal_val
IntU.orthogonal_gram_matrix
IntU.murnaghan_nakayama
IntU.tr_val
IntU.integrate_indices_orthogonal
IntU.get_pair_partitions
IntU.count_loops
IntU._expand_asymptotic
IntU.check_library
IntU.fallback_integrate
IntU.integrate_indices_gue
IntU.integrate_indices_goe
```
