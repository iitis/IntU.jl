
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
polynomials, the integrals over $SU(d)$ and $U(d)$ coincide. In the current
implementation, non-balanced expressions are handled by the same phase
invariance rule as $U(d)$ and evaluate to zero.
"""
function dSU(dim)
    _assert_no_float_param(dim, "dim", "dSU")
    return SpecialUnitary(dim)
end

function IntegrateUnitary.measure_info(measure::SpecialUnitary)
    return IntegrateUnitary.measure_info(dU(measure.dim))
end
