# Pre-computed Integral Library

`IntU.jl` includes a library of pre-computed integrals for common random matrix measures. These results are retrieved instantly, avoiding the computational overhead of the Weingarten engine.

## Unitary Group $U(d)$

The library provides the following general result for the trace of a unitary conjugation:

```math
\int tr(U A U^\dagger B) dU = \frac{tr(A) tr(B)}{d}
```

This is particularly useful when working with [Symbolic Trace Logic](symbolic_trace.md).

### Example
```julia
using IntU
@variables d
U = SymbolicMatrix(:U, false, :U)
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)

# Retrieved instantly from library
res = integrate(tr(U * A * U' * B), dU(d))
# Result: (tr(A)*tr(B))/d
```

## Gaussian Ensembles

Low-order moments of Gaussian ensembles are available for both matrix and trace-level expressions.

| Moment | GUE ($d^2$) | GOE ($d^2+d$) | GSE ($d^2-d$) |
| :--- | :--- | :--- | :--- |
| $\langle Tr(H^2) \rangle$ | $d^2$ | $d^2 + d$ | $d^2 - d$ |
| $\langle Tr(H^4) \rangle$ | $2d^3 + d$ | $2d^3 + 5d^2 + 5d$ | $2d^3 - 5d^2 + 5d$ |
| $\langle Tr(H^6) \rangle$ | $5d^4 + 10d^2$ | - | - |

### Example
```julia
using IntU
@variables d
H = SymbolicMatrix(:H, :GUE)

# <tr(H^4)>_GUE
res = integrate(tr(H^4), dGUE(d))
# Result: 2d^3 + d
```

## Fallback Mechanism
If an expression is not found in the library, `IntU.jl` automatically falls back to full symbolic integration using Weingarten calculus (for Haar) or Wick's theorem (for Gaussian), ensuring that all valid integrals can still be calculated.
