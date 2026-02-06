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

function fallback_integrate(expr, measure::SpecialUnitary)
    # Re-use the scalar integration logic from HaarMeasure
    # But we need to ensure we catch unbalanced terms first if we want to be strict.
    
    # However, the current Haar integration logic in `haar_measure.jl` (used by `fallback_integrate`)
    # ALREADY checks if n_U == n_U_bar and returns 0 if not.
    # For U(d), this is correct: \int U_{ij} dU = 0.
    # For SU(d), \int U_{ij} dSU is also 0.
    # The difference appears at order d. E.g. det(U) = 1 in SU(d), but average of det(U) in U(d) is 0.
    
    # Since our Haar logic (Weingarten) imposes balance, reusing it is safe for 
    # the subset of integrals that are non-zero in U(d).
    # For terms that are zero in U(d) but non-zero in SU(d) (like det(U)), 
    # our U(d) logic returns 0.
    # Implementing full SU(d) support would require detecting "det-like" structures.
    
    # For now, we simply map to a HaarMeasure and delegate.
    
    haar_measure = dU(measure.U, measure.dim)
    return fallback_integrate(expr, haar_measure)
end
