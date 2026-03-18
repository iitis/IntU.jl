
struct PureStateMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end

PureStateMeasure(dim) = PureStateMeasure(dim, nothing)
@doc raw"""
    dPsi(dim)

Defines the Fubini-Study measure for a **random pure state** |psi> 
distributed according to the Haar measure. 

Integration engine identifies variables via metadata tag `:psi`.
"""
dPsi(dim) = PureStateMeasure(dim)

IntU._measure_tag(::PureStateMeasure) = :psi

IntU._reconstruct_symbolic(::PureStateMeasure, d_asymp) = dPsi(d_asymp)
