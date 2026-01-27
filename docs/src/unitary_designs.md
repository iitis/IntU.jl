# Unitary Designs

IntU support integration over Unitary $t$-designs. A unitary $t$-design is an ensemble of unitary matrices $\{U_k\}$ that reproduces the Haar measure averages for all polynomials of degree at most $t$ in the elements of $U$ and $U^\dagger$.

## Usage

Use the `dDesign(U, dim, t)` function to define a measure representing a unitary $t$-design.

```@docs
dDesign
```

```julia
using IntU, Symbolics


d = 3
@variables U[1:d, 1:d]::Complex
design = dDesign(U, d, 2) # A 2-design
```

## Integration Behavior

When integrating with `integrate(expr, design)`:

1.  **Degree Check**: The integrator explicitly calculates the degree of the integrand in $U$ and $U^\dagger$.
2.  **Valid Degrees**: If both degrees are $\le t$, the integral is computed using the standard unitary Weingarten calculus (same result as Haar measure).
3.  **Invalid Degrees**: If the degree exceeds $t$, an error is thrown to indicate that the design does not support this moment.

## Example

```julia
# 2-design supports 2nd moments
expr = abs(U[1,1])^2  # degree 1
integrate(expr, design) 
# Output: 1/3 (Matches Haar)

# 2-design supports 4th moments (degree 2 in U, 2 in U*)
expr = abs(U[1,1] * U[2,2])^2 # degree 2
integrate(expr, design)
# Output: 1/8 (Matches Haar)

# 2-design DOES NOT support 6th moments
expr = abs(U[1,1])^6 # degree 3
integrate(expr, design)
# Error: Integrand degree (3, 3) exceeds design order t=2
```
