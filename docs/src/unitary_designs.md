# Unitary Designs

IntU support integration over Unitary $t$-designs.

## Theory

A unitary $t$-design is a probability distribution (ensemble) $\nu$ over unitary matrices $U \in U(d)$ such that for any polynomial $P(U, U^\dagger)$ of degree at most $t$ in the entries of $U$ and at most $t$ in the entries of $U^\dagger$:

```math
\mathbb{E}_{U \sim \nu} [P(U, U^\dagger)] = \mathbb{E}_{U \sim \text{Haar}} [P(U, U^\dagger)]
```

In other words, a $t$-design perfectly simulates the Haar measure for the first $t$ moments.
This concept is crucial in Quantum Information for efficient randomization, benchmarking, and decoupling.

## Usage

Use the `dDesign(dim, t)` function to define a measure representing a unitary $t$-design.

```@docs
dDesign
```

1. **Basic Integration using `@integrate`**

The `@integrate` macro automatically identifies `U` as the random matrix.

```julia
using IntU, Symbolics
@variables d
# E[|U_11|^2] (degree 1 in U, 1 in U*)
@integrate abs(U[1,1])^2 dDesign(d, 2)
# Output: 1/d
```

2. **Manual Integration**

```julia
using IntU, Symbolics

# Define dimension and matrix
@variables d
U = SymbolicMatrix(:U, :U, d)

# Create a 2-design measure
design = dDesign(d, 2)
integrate(abs(U[1,1])^2, design)
```

## Integration Behavior

When integrating with `integrate(expr, design)`:

1.  **Degree Check**: The integrator explicitly calculates the degree of the integrand in $U$ and $U^\dagger$.
2.  **Valid Degrees**: If both degrees are $\le t$, the integral is computed using the standard unitary Weingarten calculus (same result as Haar measure).
3.  **Invalid Degrees**: If the degree exceeds $t$, an error is thrown to indicate that the design does not support this moment. This prevents accidental reliance on non-guaranteed values.

## Example

```julia
# 1. 2-design supports 4th moments (degree 2 in U, 2 in U*)
expr2 = abs(U[1,1] * U[2,2])^2
res2 = integrate(expr2, design)
println(res2)
# Output: 1 / (d^2 - 1) * ... (Matches Haar)

# 2. 2-design DOES NOT support 6th moments (degree 3)
try
    expr3 = abs(U[1,1])^6
    integrate(expr3, design)
catch e
    println(e)
    # Error: Integrand degree (3, 3) exceeds design order t=2
end
```

## References

1.  **Dankert, C., Cleve, R., Emerson, J., & Livine, E.** (2009). Exact and approximate unitary 2-designs and their application to fidelity estimation. *Physical Review A*, 80(1), 012304.
2.  **Gross, D., Audenaert, K., & Eisert, J.** (2007). Evenly distributed unitaries: On the structure of unitary designs. *Journal of Mathematical Physics*, 48(5), 052104.
3.  **Ambainis, A., & Emerson, J.** (2007). Quantum t-designs: t-wise independence in the quantum world. *Twenty-Second Annual IEEE Conference on Computational Complexity (CCC'07)*.
