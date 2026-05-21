# Asymptotic Expansions

For large Hilbert space dimension $d$, exact Weingarten results are often
complicated rational functions of $d$. IntU.jl provides utilities to expand
these into a Laurent series in $1/d$, extracting the physically relevant
scaling behaviour.

## API

```julia
asymptotic(expr, measure, order=1)   # integrate then expand
asymptotic(expr, d, order=1)         # expand a rational function directly
```

**Arguments**

| Argument | Description |
|---|---|
| `expr` | Symbolic expression to integrate, or a rational function of `d` |
| `measure` | Integration measure (`dU`, `dGUE`, `dO`, etc.) |
| `d` | Symbolic variable to expand in (for the rational-function form) |
| `order` | Maximum power of $1/d$ to retain (default: 1) |

The first form computes the exact symbolic integral first, then performs the
Laurent expansion. The second form expands any rational expression directly,
skipping integration.

> [!NOTE]
> `order` counts the *highest* negative power of $d$ retained. For example,
> `order=3` keeps terms up to and including $1/d^3$. For most physics
> applications `order=2` or `order=3` is sufficient.

## Unitary group examples

### Leading-order trace-polynomial moment

For trace-polynomial moments that are rational in $d$, `asymptotic(expr, measure, order)`
returns the large-$d$ scaling directly:

```julia
using IntU, Symbolics
@variables d
U = symbolic_unitary(:U, d)
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)

asymptotic(tr(U * A * U' * B), dU(d), 2)
# Output: tr(A)*tr(B)/d
```

Pure-trace moments $|\mathrm{tr}(U)|^{2k}$ are a special case: their exact value
is step-function in $d$, so they should be evaluated at concrete dimension with
`integrate`.

### Entry-wise high-degree moment

The 10th moment of a single entry has exact form
$120 / (d(d+1)(d+2)(d+3)(d+4))$. The asymptotic expansion reveals the
leading scaling immediately:

```julia
using IntU, Symbolics
@variables d
Ud = SymbolicMatrix(:Ud, :U, d)

asymptotic(abs(Ud[1,1])^10, dU(d), 6)
# Output: 120/d^5 - 1200/d^6 + ...
```

The leading $1/d^5$ matches the expectation from the stable-rank estimate.

## Gaussian ensemble examples

For GUE the second trace moment $\langle \mathrm{tr}(H^2) \rangle = d^2$
grows quadratically, while $\langle \mathrm{tr}(H)^4 \rangle = 3d^2$
matches the fourth moment of a Gaussian (three Wick pairings):

```julia
using IntU, Symbolics
@variables d
H = SymbolicMatrix(:H, :GUE, d)

# Leading large-d behaviour of the GUE 4th trace moment
asymptotic(abs(tr(H))^4, dGUE(d), 2)
# Output: 3d^2
```

## Orthogonal group example

```julia
using IntU, Symbolics
@variables d
O = SymbolicMatrix(:O, :O, d)

# E[O_{11}^4] = 3/(d*(d+2)); leading large-d behaviour:
asymptotic(O[1,1]^4, dO(d), 3)
# Output: 3/d^2 - 6/d^3 + ...
```

## Expanding a rational function directly

The second form of `asymptotic` accepts any rational expression and expands it
without performing any integration. This is useful when you already have an
exact result (e.g., from the integral library or a manual calculation):

```julia
using Symbolics
@variables d

# Page's purity formula: 2n/(n^2+1) for equal subsystems of size n
page = 2d / (d^2 + 1)
asymptotic(page, d, 5)
# Output: 2/d - 2/d^3 + 2/d^5
# Leading 2/d: purity is O(1/d); maximally mixed would be exactly 1/d
```

```julia
# Exact Weingarten result 2/(d*(d+1)) expanded to four terms:
asymptotic(2 / (d * (d + 1)), d, 4)
# Output: 2/d^2 - 2/d^3 + 2/d^4 - 2/d^5
```

> [!IMPORTANT]
> ### Symbolic (d) Pitfalls
> - **Poles at small dimensions**: For Haar-related measures (Unitary,
>   Orthogonal, Circular), element-wise results are rational functions with
>   poles at small integer $d$ (typically $d < n$ for degree-$n$ moments).
> - **Non-rational results**: Pure trace moments $|\mathrm{tr}(U)|^{2k}$
>   depend on $d$ as a step function and require a concrete integer dimension.
>   Prefer `integrate(abs(tr(U))^(2k), dU(n))` over `asymptotic(expr, measure)`
>   for these cases.
> - **Removable singularities**: substituting numeric values can yield $0/0$
>   forms (e.g. at $d = 1$ or $d = 2$). Use `evaluate` to resolve these
>   automatically.

## See Also

- [Integral Library](integral_library.md) — O(1) retrieval for common results
- [Unitary Integration](unitary_integration.md) — exact Weingarten engine
- [Gaussian Ensembles](gaussian_integration.md) — GUE/GOE/GSE moments
- [QI Helpers](qi_helpers.md) — combining `asymptotic` with `partial_trace`

## References

- Collins, B. (2003). Moments and cumulants of polynomial random variables on
  unitary groups, the Itzykson–Zuber integral and free probability.
  *International Mathematics Research Notices*.
- Puchała, Z., & Miszczak, J. A. (2017). Symbolic integration with respect to
  the Haar measure on the unitary group. *Bulletin of the Polish Academy of
  Sciences: Technical Sciences*.
- Page, D. N. (1993). Average entropy of a subsystem. *Physical Review Letters*,
  71(9), 1291–1294.

See the [API Reference](api.md) for the full `asymptotic` signature.
