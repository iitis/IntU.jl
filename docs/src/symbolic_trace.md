# Symbolic Trace Integration

IntU.jl simplifies the integration of expressions involving traces of products
of random matrices, such as $\mathrm{tr}(U A U^\dagger B)$.

## Concept

Instead of writing out indices explicitly ($ \sum_{ijkl} U_{ij} A_{jk} \bar{U}_{lk} B_{li} \dots $), you can work with symbolic matrix objects.
The package provides a `tr` function (representing a lazy trace) that delays evaluation until integration time.

## Usage

1.  **Define Symbolic Matrices**: Use `SymbolicMatrix(name, is_adj, type)`.
    *   Set type to `:U` for the random unitary.
    *   Set type to `:Constant` for fixed matrices.
2.  **Define Measure**: Create a `HaarMeasure` using `dU`. The `dim` is crucial.
3.  **Construct Trace**: Use `IntU.tr` with standard matrix multiplication `*` and adjoint `'`.

```julia
using IntU, Symbolics

@variables d

# 1. Define Matrices
# Random Unitary U
U = SymbolicMatrix(:U, false, :U)
# Constant Matrix A
A = SymbolicMatrix(:A) 

# 2. Define Measure
measure = dU(U, d)

# 3. Construct Expression
# tr(U * A * U' * A)
expr = IntU.tr(U * A * U' * A)

# 4. Integrate
res = integrate(expr, measure)
println(res)
# Output: (tr(A)^2) / d
```

## Implementation Details

The system converts the lazy trace expression into a tensor network of indices,
automatically assigning input/output indices to each matrix multiplication. It
then feeds these indices into the standard Weingarten integration core.

This feature is particularly useful for checking identities in Quantum
Information Theory without getting bogged down in index hell.

## Products and Sums of Traces

`IntU.jl` supports integration of products and linear combinations of symbolic traces:

- **Multiplication**: `tr(A) * tr(B)` creates a multi-cycle trace object.
- **Addition**: `tr(A) + tr(B)` creates a `LazySum` object.

### Example

```julia
# Product of traces
expr = tr(U * A) * tr(U' * B)
integrate(expr, measure)
# Output: tr(A*B) / d

# Sum of traces
expr_sum = tr(U * A * U') + tr(B)
integrate(expr_sum, measure)
# Output: (tr(A) / d) * tr(I_d) + tr(B) = tr(A) + tr(B)
```

The underlying engine handles the "wiring" of indices across multiple trace cycles automatically.
