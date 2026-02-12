# Special Unitary Measure Integration

"""
    SpecialUnitary{D}

Dummy type to represent the Special Unitary measure (SU(d)).
"""
struct SpecialUnitary{D}
    dim::D
end

@doc raw"""
    dSU(dim)
    dSU(U::SymbolicMatrix)

Defines the Haar measure for the Special Unitary group $SU(d)$.

If called with `dim`, it integrates entries tagged with `:U` via `SymbolicMatrix(:U, :U)`.

This integration engine currently assumes the **stable range** where the dimension $d$ is large
relative to the polynomial degree, or $d$ is symbolic. In this regime, integrals over $SU(d)$
coincide with integrals over $U(d)$ for "balanced" polynomials (equal number of $U$ and $\bar{U}$ entries).
Unbalanced polynomials integrate to 0.

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
dSU(dim) = SpecialUnitary(dim)
dSU(U::SymbolicMatrix) = SpecialUnitary(U.dim)

"""
    integrate(expr, measure::SpecialUnitary)
"""
function integrate(expr::AbstractArray, measure::SpecialUnitary)
    return map(e -> integrate(e, measure), expr)
end

# Resolve ambiguity with SymbolicMatrixProduct
function integrate(P::SymbolicMatrixProduct, measure::SpecialUnitary)
    return invoke(integrate, Tuple{SymbolicMatrixProduct, Any}, P, measure)
end

function IntU.measure_info(measure::SpecialUnitary)
    # Delegates to HaarMeasure info
    return IntU.measure_info(dU(measure.dim))
end
