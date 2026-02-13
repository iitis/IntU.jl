# Unitary t-designs

struct UnitaryDesign{D,T}
    dim::D
    t::T
end

@doc raw"""
    dDesign(dim, t)

Defines a measure representing a **unitary $t$-design**. 

Integration engine identifies variables via metadata tag `:U`.
"""
dDesign(dim, t) = UnitaryDesign(dim, t)

"""
    integrate(expr, measure::UnitaryDesign)
"""
function integrate(expr::AbstractArray, measure::UnitaryDesign)
    return map(e -> integrate(e, measure), expr)
end

# Resolve ambiguities with SymbolicMatrix/SymbolicMatrixProduct
function integrate(expr::SymbolicMatrix, measure::UnitaryDesign)
    return invoke(integrate, Tuple{SymbolicMatrix, Any}, expr, measure)
end
function integrate(expr::SymbolicMatrixProduct, measure::UnitaryDesign)
    return invoke(integrate, Tuple{SymbolicMatrixProduct, Any}, expr, measure)
end

function IntU.measure_info(measure::UnitaryDesign)
    subs_dict = Dict{Any,Any}()
    matcher = MetadataMatcher(:U)
    return (subs_dict, matcher, measure.dim, (:Design, measure.t))
end
