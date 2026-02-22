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

Integration engine identifies variables via metadata tag `:U` (coincides with $U(d)$ in stable range).
"""
dSU(dim) = SpecialUnitary(dim)

function IntU.measure_info(measure::SpecialUnitary)
    # Delegates to HaarMeasure info
    return IntU.measure_info(dU(measure.dim))
end
