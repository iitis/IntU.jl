# Diagonal Unitary Matrices (Torus group) Integration

@doc raw"""
    DiagonalUnitaryMeasure{D}

Represents integration over the group of diagonal unitary matrices (the torus $T^d$).
For a diagonal unitary matrix $V$, the only non-zero entries are $V_{ii} = e^{i\theta_i}$.
Integration over $T^d$ is equivalent to independent phase integrations for each diagonal entry.
"""
struct DiagonalUnitaryMeasure{D}
    dim::D
end

@doc raw"""
    dDiagUnitary(dim)
    dDiagUnitary(V::SymbolicMatrix)

Defines the measure for the group of diagonal unitary matrices (the torus $T^d$) of dimension `dim`.

If called with `dim`, it integrates entries tagged with `:DiagUnitary` via `SymbolicMatrix(:V, :DiagUnitary)`.
"""
dDiagUnitary(dim) = DiagonalUnitaryMeasure(dim)
dDiagUnitary(V::SymbolicMatrix) = DiagonalUnitaryMeasure(V.dim)

"""
    integrate(expr, measure::DiagonalUnitaryMeasure)
"""
function integrate(expr::AbstractArray, measure::DiagonalUnitaryMeasure)
    return map(e -> integrate(e, measure), expr)
end

function IntU.measure_info(measure::DiagonalUnitaryMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = MetadataMatcher(:DiagUnitary)
    return (subs_dict, matcher, measure.dim, :DiagUnitary)
end
