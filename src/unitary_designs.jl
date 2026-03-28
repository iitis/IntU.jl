
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
function dDesign(dim, t)
    _assert_no_float_param(dim, "dim", "dDesign")
    _assert_no_float_param(t, "t", "dDesign")
    return UnitaryDesign(dim, t)
end

function IntU.measure_info(measure::UnitaryDesign)
    subs_dict = Dict{Any,Any}()
    matcher = measure.matcher === nothing ? MetadataMatcher(:U) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    dim = _assert_no_float_param(dim, "dim", "UnitaryDesign")
    t = _assert_no_float_param(measure.t, "t", "UnitaryDesign")
    return (subs_dict, matcher, dim, (:Design, t))
end
