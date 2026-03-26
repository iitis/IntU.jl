# Pre-computed Integral Library

IntU.jl maintains a library of pre-computed results for frequently encountered
integrals. Matching expressions are returned in $\mathcal{O}(1)$ time,
bypassing the Weingarten or Wick-contraction engines entirely.

> [!TIP]
> When an integral matches a library entry, retrieval is essentially
> instantaneous. For example, `|tr(U)|^{14}` with $d=10$ returns in 0.02 ms
> (see the benchmark table in the paper), compared to seconds for cold
> symbolic integration at that degree.

## Unitary Group $U(d)$

### Trace of a conjugation channel

The most commonly needed library entry is the first-moment identity:

$$\int_{U(d)} \mathrm{tr}(U A U^\dagger B)\, dU = \frac{\mathrm{tr}(A)\,\mathrm{tr}(B)}{d}.$$

```julia
using IntU, Symbolics
@variables d
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)

# Retrieved instantly from the library
@integrate tr(U * A * U' * B) dU(d)
# Output: tr(A)*tr(B) / d
```

### Trace moments

Pure trace moments $|\mathrm{tr}(U)|^{2k} = \mathrm{tr}(U)^k \cdot \mathrm{tr}(U^\dagger)^k$
are pre-computed for small $k$:

| Integral | Result |
|---|---|
| $\mathbb{E}[\,\|\mathrm{tr}(U)\|^2\,]$ | $1$ |
| $\mathbb{E}[\,\|\mathrm{tr}(U)\|^4\,]$ | $2$ |
| $\mathbb{E}[\,\|\mathrm{tr}(U)\|^6\,]$ | $6$ |
| $\mathbb{E}[\,\|\mathrm{tr}(U)\|^8\,]$ | $24$ |

The pattern $k!$ is the large-$d$ (stable-range) limit, reflecting the
Gaussian universality of $\mathrm{tr}(U)$ as $d \to \infty$. For **integer**
$d$, the library returns the exact value
$\sum_{\lambda \vdash k,\, \ell(\lambda) \le d} (f^\lambda)^2$, which equals
$k!$ when $d \ge k$ but is smaller for $d < k$. Because the dependence on $d$
is a step function (not a polynomial), trace moments require a concrete integer
dimension and will raise an error if called with symbolic $d$.

## Gaussian Ensembles

Low-order trace moments are pre-computed for GUE, GOE, and GSE.

| Moment | GUE | GOE | GSE |
|:---|:---|:---|:---|
| $\langle \mathrm{tr}(H^2) \rangle$ | $d^2$ | $d^2 + d$ | $d^2 - d$ |
| $\langle \mathrm{tr}(H^4) \rangle$ | $2d^3 + d$ | $2d^3 + 5d^2 + 5d$ | $2d^3 - 5d^2 + 5d$ |
| $\langle \mathrm{tr}(H^6) \rangle$ | $5d^4 + 10d^2$ | — | — |

```julia
using IntU, Symbolics
@variables d

# GUE 4th moment — retrieved from library
@integrate tr(H^4) dGUE(d)
# Output: 2d^3 + d

# GOE 2nd moment
@integrate tr(H^2) dGOE(d)
# Output: d^2 + d
```

The library also covers element-wise second moments:

```julia
# GUE: E[H_{ij} H_{kl}] = delta_{il}*delta_{jk}
# This feeds into higher-level moment matching automatically.
```

## Ginibre Ensembles

Trace moments for GinUE, GinOE, and GinSE are cached for low degrees:

| Moment | GinUE |
|:---|:---|
| $\langle \mathrm{tr}(G G^\dagger) \rangle$ | $d^2$ |
| $\langle \mathrm{tr}(G G^\dagger)^2 \rangle$ | $d^4 + d^2$ |
| $\langle \mathrm{tr}((G G^\dagger)^2) \rangle$ | $2d^3$ |

```julia
# GinUE second moment — O(1) retrieval
@integrate tr(G * G') dGinUE(d)
# Output: d^2
```

## Fallback Mechanism

If an expression is not found in the library, IntU.jl automatically falls back
to full symbolic integration (Weingarten calculus for compact groups, Wick
contraction for Gaussian/Ginibre ensembles). No user action is required — the
library is consulted transparently on every `integrate` call.

## See Also

- [Symbolic Trace Logic](symbolic_trace.md) — index-free trace expressions that
  benefit most from library matching
- [Gaussian Ensembles](gaussian_integration.md) — Wick-contraction details
- [Asymptotic Expansions](asymptotic.md) — expanding library results to large $d$
