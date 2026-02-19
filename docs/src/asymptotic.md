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

## Example

Evaluating the fourth moment of a matrix entry $|U_{11}|^4$:

```julia
@variables d
U = SymbolicMatrix(:U, :U, d)
res = asymptotic(abs(U[1,1])^4, dU(d), 4)
# Output: 2/d^2 - 2/d^3 + 2/d^4
```

This approximation is useful for checking convergence properties or
leading-order behavior in high-dimensional quantum systems.

## References

- Collins, B. (2003). Moments and Cumulants of Polynomial random variables on
  unitary groups, the Itzykson-Zuber integral and free probability.
  *International Mathematics Research Notices*.
- Puchała, Z., & Miszczak, J. A. (2017). Symbolic integration with respect to
  the Haar measure on the unitary group. *Bulletin of the Polish Academy of
  Sciences: Technical Sciences*.
