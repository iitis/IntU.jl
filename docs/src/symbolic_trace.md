# Symbolic Trace Integration

IntU.jl simplifies the integration of expressions involving traces of products
of random matricies, such as $\mathrm{tr}(U A U^\dagger B)$.

## Concept

Instead of writing out indices explicitly ($ \sum_{ijkl} U_{ij} A_{jk}
\bar{U}_{lk} B_{li} \dots $), you can work with symbolic matrix objects. The
package provides a `tr_lazy` function (exported as `tr` via `IntU.tr` or
explicitly `LinearAlgebra.tr` overload on symbolic types in older versions) that
delays evaluation until integration time.

## Usage

1.  **Define Symbolic Matrices**: Use `SymbolicMatrix(name)`.
2.  **Define Matrix Variables**: $U$ (random) and $A, B, C$ (constants).
3.  **Construct Trace**: Use standard `*` and adjoint `'` operations.

```julia
using IntU, Symbolics

@variables d
# Dummy usage for measure definition
@variables u_dummy[1:1, 1:1]
measure = dU(u_dummy, d)

U = SymbolicMatrix(:U, false, :U) # Variable name U, is_constant=false, ID=:U
A = SymbolicMatrix(:A)            # Constant matrix A

expr = IntU.tr(U * A * U' * A)
integrate(expr, measure)
# Output: (tr(A)^2)/d
```

## Implementation

The system converts the lazy trace expression into a tensor network of indices,
automatically assigning input/output indices to each matrix multiplication. It
then feeds these indices into the standard Weingarten integration core.

This feature is particularly useful for checking identities in Quantum
Information Theory.

## Products and Sums of Traces

`IntU.jl` supports integration of products and linear combinations of symbolic traces:

- **Multiplication**: `tr(A) * tr(B)` creates a multi-cycle trace object.
- **Addition**: `tr(A) + tr(B)` creates a `LazySum` object.

### Example

```julia
# Product of traces
expr = tr(U * A) * tr(U' * B)
integrate(expr, measure)
# Output: tr(A B) / d

# Sum of traces
expr_sum = tr(U * A * U') + tr(B)
integrate(expr_sum, measure)
# Output: tr(A)/d * tr(I) + tr(B) = tr(A) + tr(B)
```

The underlying engine handles the "wiring" of indices across multiple trace cycles automatically.
