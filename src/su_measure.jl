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

Defines the Haar measure for the Special Unitary group $SU(d)$.

Integration engine identifies variables via metadata tag `:U` (coincides with $U(d)$ in stable range).
"""
dSU(dim) = SpecialUnitary(dim)

"""
    integrate(expr, measure::SpecialUnitary)
"""
function integrate(expr::AbstractArray, measure::SpecialUnitary)
    return map(e -> integrate(e, measure), expr)
end

# Resolve ambiguity with SymbolicMatrix/SymbolicMatrixProduct
function integrate(expr::SymbolicMatrix, measure::SpecialUnitary)
    return invoke(integrate, Tuple{SymbolicMatrix,Any}, expr, measure)
end
function integrate(P::SymbolicMatrixProduct, measure::SpecialUnitary)
    return invoke(integrate, Tuple{SymbolicMatrixProduct,Any}, P, measure)
end

function IntU.measure_info(measure::SpecialUnitary)
    # Delegates to HaarMeasure info
    return IntU.measure_info(dU(measure.dim))
end
