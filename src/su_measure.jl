# Special Unitary Measure Integration

"""
    SpecialUnitary{D}

Dummy type to represent the Special Unitary measure (SU(d)).
"""
struct SpecialUnitary{D} <: AbstractMeasure
    dim::D
end

@doc raw"""
    dSU(dim)

Defines the Haar measure for the Special Unitary group $SU(d)$.

Integration for $SU(d)$ is performed via the $U(d)$ measure. For balanced 
polynomials, the integrals over $SU(d)$ and $U(d)$ coincide.
"""
dSU(dim) = SpecialUnitary(dim)

function IntU.measure_info(measure::SpecialUnitary)
    # Delegates to HaarMeasure info
    return IntU.measure_info(dU(measure.dim))
end
