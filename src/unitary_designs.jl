# Unitary t-designs

struct UnitaryDesign{D,T,M}
    dim::D
    t::T
    matcher::M
end
UnitaryDesign(dim, t) = UnitaryDesign(dim, t, nothing)

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
    return invoke(integrate, Tuple{SymbolicMatrix,Any}, expr, measure)
end
function integrate(expr::SymbolicMatrixProduct, measure::UnitaryDesign)
    return invoke(integrate, Tuple{SymbolicMatrixProduct,Any}, expr, measure)
end

function IntU.measure_info(measure::UnitaryDesign)
    subs_dict = Dict{Any,Any}()
    matcher = measure.matcher === nothing ? MetadataMatcher(:U) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, (:Design, measure.t))
end
