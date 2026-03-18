
struct UnitaryDesign{D,T,M} <: AbstractMeasure
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

function IntU.measure_info(measure::UnitaryDesign)
    subs_dict = Dict{Any,Any}()
    matcher = measure.matcher === nothing ? MetadataMatcher(:U) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, (:Design, measure.t))
end
