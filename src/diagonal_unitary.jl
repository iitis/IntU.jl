
@doc raw"""
    DiagonalUnitaryMeasure{D}

Represents integration over the group of diagonal unitary matrices (the torus $T^d$).
For a diagonal unitary matrix $V$, the only non-zero entries are $V_{ii} = e^{i\theta_i}$.
Integration over $T^d$ is equivalent to independent phase integrations for each diagonal entry.
"""
struct DiagonalUnitaryMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end
DiagonalUnitaryMeasure(dim) = DiagonalUnitaryMeasure(dim, nothing)

@doc raw"""
    dDiagUnitary(dim)

Defines the measure for the group of diagonal unitary matrices (the torus $T^d$) of dimension `dim`.
Integration engine identifies variables via metadata tag `:DiagUnitary`.
"""
function dDiagUnitary(dim)
    _assert_no_float_param(dim, "dim", "dDiagUnitary")
    return DiagonalUnitaryMeasure(dim)
end

IntegrateUnitary._measure_tag(::DiagonalUnitaryMeasure) = :DiagUnitary
