# Pure state integration

# Type to represent integration over pure states |psi>
# Type to represent integration over pure states |psi>
struct PureStateMeasure{D,M} <: AbstractMeasure
    dim::D
    matcher::M
end

# Constructor for backward compatibility
PureStateMeasure(dim) = PureStateMeasure(dim, nothing)
@doc raw"""
    dPsi(dim)

Defines the Fubini-Study measure for a **random pure state** |psi> 
distributed according to the Haar measure. 

Integration engine identifies variables via metadata tag `:psi`.
"""
dPsi(dim) = PureStateMeasure(dim)

function IntU.measure_info(measure::PureStateMeasure)
    subs_dict = Dict{Any,Any}()
    matcher = measure.matcher === nothing ? MetadataMatcher(:psi) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, :psi)
end

"""
    asymptotic(expr, measure::PureStateMeasure, order=1)

Returns the series expansion of the integral in powers of `1/d`.
"""
function asymptotic(expr, measure::PureStateMeasure, order = 1)
    d = measure.dim
    if d isa SymbolicMatrix
        d = d.dim
    end

    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dPsi(d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end
