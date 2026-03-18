# Symbolic Trace Integration

IntU.jl allows integration of expressions involving traces of products of random
matrices — such as $\mathrm{tr}(U A U^\dagger B)$ — without writing out tensor
indices explicitly.

## Concept

Index-free notation avoids manual bookkeeping like

$$\sum_{i,j,k,l} U_{ij}\, A_{jk}\, \bar{U}_{lk}\, B_{li}.$$

Instead, you work with symbolic matrix objects. The `tr()` function returns a
`LazyTrace` that defers index assignment to integration time, when the package
automatically reconstructs the Weingarten graph from the cycle structure.

## `tr` and `tr_lazy`

Two entry points build lazy trace expressions:

| Function | Use when |
|---|---|
| `tr(expr)` | Building expressions with `*` — the natural choice |
| `tr_lazy(factors)` | You have a pre-assembled `Vector` of matrix factors |

Both return a [`LazyTrace`](@ref) object. Products of traces (`tr(A) * tr(B)`)
produce a multi-cycle `LazyTrace`; sums (`tr(A) + tr(B)`) produce a
[`LazySum`](@ref).

## Basic usage

### With `@integrate`

```julia
using IntU, Symbolics
@variables d

# E[tr(U A U' B)] = tr(A)*tr(B) / d
@integrate tr(U * A * U' * B) dU(d)
```

### Manual construction

```julia
using IntU, Symbolics
@variables d

U = SymbolicMatrix(:U, :U)
A = SymbolicMatrix(:A)

expr = IntU.tr(U * A * U' * A)
integrate(expr, dU(d))
# Output: tr(A)^2 / d
```

### Using `tr_lazy` with a vector of factors

```julia
using IntU, Symbolics
@variables d
U = symbolic_unitary(:U, d)
A = SymbolicMatrix(:A)

# Equivalent to tr(U * A * U')
expr = tr_lazy([U, A, U'])
integrate(expr, dU(d))
# Output: tr(A) / d
```

## Multi-cycle traces (products of traces)

Multiplying two `LazyTrace` objects creates a two-cycle expression:

```julia
# E[tr(U*A) * tr(U'*B)] — two separate trace cycles
@integrate tr(U * A) * tr(U' * B) dU(d)
# Output: tr(A*B) / d
```

The engine wires the indices across cycles automatically using the Weingarten
matching rules.

## `LazySum` — linear combinations of traces

Addition creates a [`LazySum`](@ref), which is integrated term by term:

```julia
@variables d
# E[tr(U*A*U') + tr(U*B*U')]
@integrate tr(U * A * U') + tr(U * B * U') dU(d)
# Output: tr(A)/d + tr(B)/d
```

## Advanced: 2-design identity

Expressions with two $U$ and two $U^\dagger$ factors exercise the degree-2
Weingarten sum. The following identity is the defining property of a unitary
2-design:

```julia
@variables d
res = @integrate tr(U*A*U'*B*U*C*U'*D) dU(d)
# Output: (tr(A*C)*tr(B*D) + ...) / (d^2 - 1)
# (full expression is a rational function of d and traces of A,B,C,D)
```

## Library fast-path

Common patterns — particularly single-channel conjugations and trace moments —
are matched against the [Pre-computed Integral Library](integral_library.md)
before the full Weingarten engine is invoked. Matched expressions return in
$\mathcal{O}(1)$.

```julia
# tr(U A U' B) hits the library and returns immediately
@integrate tr(U * A * U' * B) dU(d)
# Output: tr(A)*tr(B) / d  (O(1) retrieval)
```

## Implementation note

When `integrate` receives a `LazyTrace`, it:

1. Identifies the positions of $U$ and $U^\dagger$ factors within each cycle.
2. Assigns effective indices to the intervening constant matrices.
3. Constructs the Weingarten contraction graph from the resulting index tuples.
4. Dispatches to the standard integration engine.

This avoids the $\mathcal{O}((k!)^2)$ cost of eager index expansion for
high-degree expressions.

## See Also

- [Integral Library](integral_library.md) — O(1) fast-paths for common patterns
- [Unitary Integration](unitary_integration.md) — Weingarten calculus background
- [ITensors Integration](itensors.md) — graphical Weingarten for tensor networks

See [`IntU.LazyTrace`](@ref), [`IntU.LazySum`](@ref), and [`tr_lazy`](@ref) in the [API Reference](api.md).
