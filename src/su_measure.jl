# Special Unitary Measure Integration

"""
    SpecialUnitary{M,D}

Dummy type to represent the Special Unitary measure (SU(d)).
"""
struct SpecialUnitary{M,D}
    U::M
    dim::D
end

"""
    dSU(U, dim)
    dSU(dim)

Defines the Haar measure for the Special Unitary group SU(d).

This integration engine currently assumes the **stable range** where the dimension \$d\$ is large
relative to the polynomial degree, or \$d\$ is symbolic. In this regime, integrals over \$SU(d)\$
coincide with integrals over \$U(d)\$ for "balanced" polynomials (equal number of \$U\$ and \$\\bar{U}\$ entries).
Unbalanced polynomials integrate to 0 (unless the imbalance is a multiple of \$d\$, which requires explicit
Levita-Civita tensor handling, not yet implemented).

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
dSU(U, dim) = SpecialUnitary(U, dim)
dSU(dim) = SpecialUnitary(nothing, dim)

"""
    integrate(expr, measure::SpecialUnitary)
"""
function integrate(expr::AbstractArray, measure::SpecialUnitary)
    return map(e -> integrate(e, measure), expr)
end

function IntU.measure_info(measure::SpecialUnitary)
    # Delegates to HaarMeasure info
    return IntU.measure_info(dU(measure.U, measure.dim))
end
