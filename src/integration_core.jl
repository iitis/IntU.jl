_to_Num(z::Complex) = Complex(Num(real(z)), Num(imag(z)))
_to_Num(z) = Num(z)

"""
    AbstractMeasure

Abstract base type for all integration measures (Haar, Gaussian, Circular, etc.).
Provides generic dispatch for element-wise integration of arrays and 
ambiguity resolution for `SymbolicMatrix` and `SymbolicMatrixProduct`.
"""
abstract type AbstractMeasure end


function integrate(expr::LazySum, measure::AbstractMeasure)
    n_terms = length(expr.terms)
    if n_terms > 1
        p = Progress(n_terms; dt=10.0, desc="Integrating lazy terms... ")
        res = zero(Num)
        for t in expr.terms
            res += integrate(t, measure)
            next!(p)
        end
        return res
    end
    return sum(t -> integrate(t, measure), expr.terms)
end

function integrate(lp::LazyPower, measure::AbstractMeasure)
    e = Symbolics.unwrap(lp.exponent)
    if e isa Integer
        return integrate(lp.base^Int(e), measure)
    elseif e isa Number && isinteger(e)
        return integrate(lp.base^Int(e), measure)
    end
    throw(IntegrationError("Non-integer power of trace integration not supported yet: ($(lp.base))^$(lp.exponent)"))
end

function _get_measure_dim(measure)
    if hasproperty(measure, :dim)
        d = measure.dim
        if d isa SymbolicMatrix
            return d.dim
        end
        return d
    end
    return nothing
end

"""
    integrate(expr::SymbolicMatrix, measure)

Integrate a SymbolicMatrix as a whole. Returns a matrix of results if the 
dimension in `measure` is a concrete integer.
"""
function integrate(A::SymbolicMatrix, measure::AbstractMeasure)
    dim = _get_measure_dim(measure)
    if dim isa Integer
        res = Matrix{Any}(undef, dim, dim)
        fill!(res, 0)
        p = Progress(dim*dim; dt=10.0, desc="Integrating matrix elements... ")
        for i = 1:dim
            for j = 1:dim
                res[i, j] = integrate(A[i, j], measure)
                next!(p)
            end
        end
        return res
    end
    error(
        "Direct integration of SymbolicMatrix requires a numeric dimension in the measure.",
    )
end

_get_integration_tag(m::MetadataMatcher) = m.type_tag

function _has_integration_variable(expr, tag::Symbol)
    if expr isa SymbolicMatrix
        return expr.special_type === tag
    elseif expr isa SymbolicKron
        return _has_integration_variable(expr.A, tag) || _has_integration_variable(expr.B, tag)
    elseif expr isa SymbolicMatrixProduct
        return any(f -> _has_integration_variable(f, tag), expr.factors)
    else
        return false
    end
end

"""
    integrate(expr::SymbolicMatrixProduct, measure)

Integrate a product of SymbolicMatrices. Returns a matrix of results if the 
dimension in `measure` is a concrete integer. It skips expansion and returns the product 
itself if it does not contain the integration variable.
"""
function integrate(P::SymbolicMatrixProduct, measure::AbstractMeasure)
    if isempty(P.factors)
        return Num(1)
    end

    # Fast path: if the product does not contain the integration variable at all,
    # it acts as a constant, so the integral is just the product itself.
    _, matcher, _, _ = measure_info(measure)
    
    # We define a helper to safely extract the tag from different matchers
    if matcher isa MetadataMatcher || (isdefined(Main, :SymbolicMatcher) && matcher isa Main.SymbolicMatcher) || hasproperty(matcher, :tag) || hasproperty(matcher, :type_tag)
        tag = hasproperty(matcher, :type_tag) ? matcher.type_tag : matcher.tag
        if !_has_integration_variable(P, tag)
            return P
        end
    end

    dim_measure = _get_measure_dim(measure)
    n_factors = length(P.factors)
    
    # inner_dims[i] is the shared dimension between factor i-1 and factor i
    # inner_dims[1] is number of rows, inner_dims[n+1] is number of columns
    inner_dims = Vector{Any}(fill(nothing, n_factors + 1))
    
    # Pass 1: Collect known dimensions
    for (i, f) in enumerate(P.factors)
        fr, fc = size(f)
        fr_un = Symbolics.unwrap(fr)
        fc_un = Symbolics.unwrap(fc)
        
        if fr_un isa Integer && fr_un != typemax(Int)
            inner_dims[i] = fr_un
        end
        if fc_un isa Integer && fc_un != typemax(Int)
            inner_dims[i+1] = fc_un
        end
    end
    
    # Pass 2: Fallback to measure dimension for any remaining unknowns
    for i in 1:(n_factors + 1)
        if inner_dims[i] === nothing
            inner_dims[i] = dim_measure
        end
    end

    nr = inner_dims[1]
    nc = inner_dims[end]
    
    nr_un = Symbolics.unwrap(nr)
    nc_un = Symbolics.unwrap(nc)

    if nr_un isa Integer && nc_un isa Integer
        nr = Int(nr_un)
        nc = Int(nc_un)

        # Pre-calculate factor matrices to avoid repeated getindex overhead
        mats = []
        for (i, f) in enumerate(P.factors)
            cur_r = Symbolics.unwrap(inner_dims[i])
            cur_c = Symbolics.unwrap(inner_dims[i+1])

            if !(cur_r isa Integer && cur_c isa Integer)
                error(
                    "Factor $f has non-numeric size ($cur_r, $cur_c). cannot expand product for matrix-valued integration. Use tr() for scalar results with symbolic dimensions.",
                )
            end

            push!(mats, [f[r, c] for r = 1:Int(cur_r), c = 1:Int(cur_c)])
        end
        # Calculate the full symbolic product once
        res_mat = reduce(*, mats)

        res = Matrix{Any}(undef, nr, nc)
        fill!(res, 0)
        p = Progress(nr*nc; dt=10.0, desc="Integrating matrix product elements... ")
        for i = 1:nr
            for j = 1:nc
                res[i, j] = integrate(res_mat[i, j], measure)
                next!(p)
            end
        end
        return res
    end
    error("Direct matrix-valued integration of SymbolicMatrixProduct requires numeric result dimensions (got $nr x $nc). Use tr() for scalar results.")
end

"""
    integrate(expr::AbstractArray, measure)

Performs element-wise integration of a matrix or array of expressions.
Returns an array of the same shape containing integrated values.
"""
function integrate(expr::AbstractArray, measure::AbstractMeasure)
    return map(t -> integrate(t, measure), expr)
end

"""
    integrate(expr, measure)

Top-level integration function. It first checks the [Pre-computed Integral Library](@ref) 
for instant results. If not found, it calls `fallback_integrate` for the specific measure.
"""
function integrate(expr, measure::AbstractMeasure)
    lib_res = check_library(expr, measure)
    if lib_res !== nothing
        return lib_res
    end

    return fallback_integrate(expr, measure)
end

function _integrate_core(
    expr::Complex{Num},
    dim,
    subs_dict,
    matcher::AbstractIndexMatcher,
    measure_type = :U,
)
    # Split into real and imaginary parts to avoid recursion
    re = Symbolics.real(expr)
    im_part = Symbolics.imag(expr)

    int_re = _integrate_core(re, dim, subs_dict, matcher, measure_type)
    int_im = _integrate_core(im_part, dim, subs_dict, matcher, measure_type)

    return int_re + 1im * int_im # standard complex number result
end

function fallback_integrate(expr, measure)
    info = measure_info(measure)
    if info !== nothing
        subs_dict, matcher, dim, measure_type = info
        return _integrate_core(expr, dim, subs_dict, matcher, measure_type)
    end
    _manual_fallback(expr, measure)
end

function _manual_fallback(expr, measure)
    error("Fallback integrate not implemented for this measure: $(typeof(measure))")
end

"""
    _measure_tag(measure)

Returns the integration tag symbol (e.g. `:U`, `:O`, `:GUE`) for a measure.
Concrete measure types should define a one-liner method:

    _measure_tag(::MyMeasure) = :MyTag

This is used by the generic `measure_info` default below.
"""
function _measure_tag end

"""
    _reconstruct_symbolic(measure, d_asymp)

Reconstruct a measure with a symbolic dimension variable `d_asymp`.
Used by the generic `asymptotic` method.
"""
function _reconstruct_symbolic end

"""
    measure_info(measure)

Returns `(subs_dict, matcher, dim, measure_type)` for a given measure.
Subtypes participate in the unified integration flow by defining `_measure_tag`.
Measures with custom logic can override `measure_info` directly.
"""
function measure_info(measure)
    return nothing
end

function measure_info(measure::AbstractMeasure)
    tag = _measure_tag(measure)
    subs_dict = Dict{Any,Any}()
    matcher = if hasproperty(measure, :matcher)
        measure.matcher === nothing ? MetadataMatcher(tag) : measure.matcher
    else
        MetadataMatcher(tag)
    end
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    return (subs_dict, matcher, dim, tag)
end

