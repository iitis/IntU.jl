# Symbolic Trace Integration

IntU.jl simplifies the integration of expressions involving traces of products
of random matrices, such as $\mathrm{tr}(U A U^\dagger B)$.

## Concept

Instead of writing out indices explicitly ($ \sum_{ijkl} U_{ij} A_{jk} \bar{U}_{lk} B_{li} \dots $), you can work with symbolic matrix objects.
The package provides a `tr` function (representing a lazy trace) that delays evaluation until integration time.

## Usage

1. **Basic Integration using `@integrate`**

The `@integrate` macro works seamlessly with traces.

```julia
using IntU, Symbolics
@variables d
# E[tr(U A U' B)] = tr(A)*tr(B) / d
@integrate tr(U * A * U' * B) dU(d)
```

2. **Manual Integration**

For more control, or when dealing with complex expressions, you can declare symbols explicitly.

```julia
using IntU, Symbolics

@variables d

# 1. Define Matrices
# Random Unitary U
U = SymbolicMatrix(:U, :U)
# Constant Matrix A
A = SymbolicMatrix(:A) 

# 2. Define Measure
measure = dU(d)

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
using IntU, Symbolics
@variables d
# Product of traces
@integrate tr(U * A) * tr(U' * B) dU(d)
# Output: tr(A*B) / d
```

The underlying engine handles the "wiring" of indices across multiple trace cycles automatically.

## Advanced: 2-Design Identities

For deeper insights into quantum chaos or state designs, you can integrate more complex expressions. A well-known identity involving two $U$ and two $U^\dagger$ factors (a 2-design property) is:

```julia
using IntU, Symbolics
@variables d
# E[tr(U A U' B U C U' D)]
res = @integrate tr(U*A*U'*B*U*C*U'*D) dU(d)
# Simplified result:
# (tr(AC)tr(BD) - d*tr(A)tr(BD)tr(C) - d*tr(AC)tr(B)tr(D) + tr(A)tr(B)tr(C)tr(D)) / (d - d^3)
```

This showcases how the symbolic engine resolves complex Weingarten sums into readable rational functions of $d$ and traces of the constant matrices.
