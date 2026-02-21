# Diagonal Unitary Matrices (Torus group) Integration

@doc raw"""
    DiagonalUnitaryMeasure{D}

Represents integration over the group of diagonal unitary matrices (the torus $T^d$).
For a diagonal unitary matrix $V$, the only non-zero entries are $V_{ii} = e^{i\theta_i}$.
Integration over $T^d$ is equivalent to independent phase integrations for each diagonal entry.
"""
struct DiagonalUnitaryMeasure{D,M}
    dim::D
    matcher::M
end
DiagonalUnitaryMeasure(dim) = DiagonalUnitaryMeasure(dim, nothing)

@doc raw"""
    dDiagUnitary(dim)

Defines the measure for the group of diagonal unitary matrices (the torus $T^d$) of dimension `dim`.
Integration engine identifies variables via metadata tag `:DiagUnitary`.
"""
dDiagUnitary(dim) = DiagonalUnitaryMeasure(dim)

"""
    integrate(expr, measure::DiagonalUnitaryMeasure)
"""
function integrate(expr::AbstractArray, measure::DiagonalUnitaryMeasure)
    return map(e -> integrate(e, measure), expr)
end

# Resolve ambiguities with SymbolicMatrix/SymbolicMatrixProduct
function integrate(expr::SymbolicMatrix, measure::DiagonalUnitaryMeasure)
    return invoke(integrate, Tuple{SymbolicMatrix,Any}, expr, measure)
end
function integrate(expr::SymbolicMatrixProduct, measure::DiagonalUnitaryMeasure)
    return invoke(integrate, Tuple{SymbolicMatrixProduct,Any}, expr, measure)
end

function IntU.measure_info(measure::DiagonalUnitaryMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = measure.matcher === nothing ? MetadataMatcher(:DiagUnitary) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, :DiagUnitary)
end
