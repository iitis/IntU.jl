# Stiefel Manifold Integration

IntU.jl supports integration over the **Stiefel manifold** $V_k(\mathbb{C}^d)$,
the set of all $d \times k$ matrices with orthonormal columns:

$$V_k(\mathbb{C}^d) = \{ V \in \mathbb{C}^{d \times k} \mid V^\dagger V = I_k \}.$$

## Mathematical Background

### Connection to the Unitary Group

The uniform (Haar) measure on $V_k(\mathbb{C}^d)$ is induced by the Haar
measure on $U(d)$. A Haar-random Stiefel frame is obtained by taking the first
$k$ columns of a Haar-random unitary:

$$V = U_{1:d,\; 1:k}, \quad U \sim \mathrm{Haar}(U(d)).$$

Consequently, every polynomial integral over the Stiefel manifold reduces to
unitary Weingarten calculus:

$$\int_{V_k(\mathbb{C}^d)} f(V)\, dV = \int_{U(d)} f\!\left(U_{1:d,\, 1:k}\right) dU.$$

IntU.jl exploits this reduction internally — the full symbolic Weingarten
engine is automatically available for Stiefel integrals.

### Special case: $k = 1$ and pure states

When $k = 1$, the Stiefel manifold reduces to the complex unit sphere
$S^{2d-1}$, and the Haar measure coincides with the Fubini–Study measure on
Haar-random pure quantum states. In code, `dStiefel(d, 1)` and `dPsi(d)` are
equivalent. See [Pure States](pure_states.md).

> [!IMPORTANT]
> The rank $k$ must satisfy $k \le d$. Passing `k > d` raises an
> `ArgumentError` at measure construction time.

## Usage

Use the `dStiefel(d, k)` measure. The `@integrate` macro automatically
identifies `V` as the random Stiefel matrix.

### Basic integration

```julia
using IntU, Symbolics
@variables d

# E[|V_{11}|^2] = 1/d
@integrate abs(V[1, 1])^2 dStiefel(d, 2)
# Output: 1/d
```

The result $1/d$ matches the unitary case: the first column of a Stiefel frame
is a Haar-random unit vector, so each component has average squared magnitude $1/d$.

### Off-diagonal element (same column, different rows)

```julia
# E[|V_{11}|^2 * |V_{21}|^2]: two entries from the same column
@integrate abs(V[1, 1])^2 * abs(V[2, 1])^2 dStiefel(d, 2)
# Output: 1 / (d*(d+1))
```

This is strictly less than $1/d^2$, reflecting the negative correlation
(anti-correlation) between entries of the same column imposed by the
unit-norm constraint.

### Cross-column correlation

```julia
# E[|V_{11}|^2 * |V_{12}|^2]: entries from different columns
@integrate abs(V[1, 1])^2 * abs(V[1, 2])^2 dStiefel(d, 2)
# Output: 1 / (d*(d+1))
```

The same value as the within-column result above follows from the unitary
invariance of the measure.

### Second-moment matrix

Matrix integration confirms the isotropic structure $\mathbb{E}[V V^\dagger] = \frac{k}{d} I_d$:

```julia
# E[V * V'] for d=3, k=2 (concrete dimensions required for matrix integration)
@integrate V * V' dStiefel(3, 2)
# Output: (2/3) * I_3
```

### Asymptotic expansion

Large-$d$ expansions are fully supported. Pass the exact result as a rational
function of `d` to `asymptotic`, or integrate and expand in one step:

```julia
using IntU, Symbolics
@variables d

# Exact: 1/(d*(d+1)); leading large-d behaviour:
asymptotic(1 / (d * (d + 1)), d, 3)
# Output: 1/d^2 - 1/d^3 + 1/d^4
```

### Manual integration

For more control, create the `SymbolicMatrix` explicitly:

```julia
using IntU, Symbolics
@variables d
V = SymbolicMatrix(:V, :V)
integrate(abs(V[1, 1])^2, dStiefel(d, 2))
# Output: 1/d
```

## Potential Pitfalls

> [!IMPORTANT]
> ### Symbolic (d) Pitfalls
> - Results are rational functions of $d$ with poles at small integer values
>   (typically $d < k$ for degree-$k$ moments).
> - Use `evaluate` to resolve removable singularities when substituting
>   specific numeric dimensions.

## See Also

- [Pure States](pure_states.md) — the $k = 1$ special case (`dPsi`)
- [Unitary Integration](unitary_integration.md) — underlying Weingarten engine
- [Asymptotic Expansions](asymptotic.md) — large-$d$ limit analysis

## References

- Edelman, A., Arias, T. A., & Smith, S. T. (1998). The geometry of algorithms
  with orthogonality constraints. *SIAM Journal on Matrix Analysis and
  Applications*, 20(2), 303–353.
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar
  measure on unitary, orthogonal and symplectic groups. *Communications in
  Mathematical Physics*, 264(3), 773–795.

## API Reference

```@docs
IntU.StiefelMeasure
dStiefel
```
