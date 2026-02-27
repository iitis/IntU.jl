# Asymptotic Expansions

For large Hilbert space dimension $d$, exact Weingarten results can be
complicated rational functions. IntU.jl provides utilities to expand these
results as a Taylor series in $1/d$.

## Usage

```julia
asymptotic(expr, measure, order=1)
```

- **expr**: The symbolic expression to integrate.
- **measure**: The integration measure (Haar, PureState, GinUE, etc.).
- **order**: The maximum power of $1/d$ to retain (default 1).
- **Mechanism**: The function computes the exact symbolic integral and then performs a Taylor expansion of the result in powers of $1/d$.

## Example: The "Painful" Showcase

Standard Weingarten integration involves a sum over permutations that grows as $(n!)^2$ for general expressions. For high-degree moments, the exact result can be a very complicated rational function, while the asymptotic expansion remains clean and easy to interpret.

### High-Degree Trace Moment

Consider the 12th moment of the trace $|\text{Tr}(U)|^{12}$ (degree $n=6$):

```julia
using IntU, Symbolics
@variables d
U = SymbolicMatrix(:U, :U, d)

# Integrating |tr(U)|^12 (degree 6 in U and 6 in U*)
# Exact integration is "painful" due to the combinatorics (n=6)
# but asymptotic provides a quick results:
as = asymptotic(abs(tr(U))^12, dU(d), 1)
# Output: 720.0 (matches n! for large d)
```

### Complex Rational Functions

For expressions that are not "pure traces", the advantage is even clearer. For example, the 10th moment of a single entry $|U_{11}|^{10}$:

```julia
res = asymptotic(abs(U[1,1])^10, dU(d), 6)
# Output: 120/d^5 - 1200/d^6 + 7800/d^7
```

While the exact result is the rational function $\frac{120}{d(d+1)(d+2)(d+3)(d+4)}$, the asymptotic form highlights the leading $1/d^5$ behavior and subsequent corrections which are often sufficient for physics applications.

> [!IMPORTANT]
> ### Symbolic (d) Pitfalls
> - **Small Dimensions**: For Haar-related measures (Unitary, Orthogonal, Circular), results are rational functions with poles at small $d$ (typically $d < n$ for degree $n$ moments).
> - **Removable Singularities**: Substituting numeric values can yield $0/0$ forms (e.g., at $d=1$ or $d=2$).
> - **Automatic Handling**: `IntU.jl`'s `evaluate` function automatically simplifies expressions to resolve removable singularities when a denominator evaluates to zero.

This approximation is useful for checking convergence properties or
leading-order behavior in high-dimensional quantum systems.

## References

- Collins, B. (2003). Moments and Cumulants of Polynomial random variables on
  unitary groups, the Itzykson-Zuber integral and free probability.
  *International Mathematics Research Notices*.
- Puchała, Z., & Miszczak, J. A. (2017). Symbolic integration with respect to
  the Haar measure on the unitary group. *Bulletin of the Polish Academy of
  Sciences: Technical Sciences*.
