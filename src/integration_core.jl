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
        p = Progress(n_terms; dt = 10.0, desc = "Integrating lazy terms... ")
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
    throw(
        IntegrationError(
            "Non-integer power of trace integration not supported yet: ($(lp.base))^$(lp.exponent)",
        ),
    )
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
    _validate_measure_discrete_params(measure)

Fail fast on float-like discrete measure parameters, even for code paths that
bypass `measure_info` (e.g. precomputed-library or specialized integration
methods). This keeps constructor and integration-time policy aligned.
"""
function _validate_measure_discrete_params(measure::AbstractMeasure)
    ctx = "measure $(nameof(typeof(measure)))"

    if hasproperty(measure, :dim)
        dim = getproperty(measure, :dim)
        if dim isa SymbolicMatrix
            dim = dim.dim
        end
        _assert_no_float_param(dim, "dim", ctx)
    end

    if hasproperty(measure, :k)
        _assert_no_float_param(getproperty(measure, :k), "k", ctx)
    end

    if hasproperty(measure, :t)
        _assert_no_float_param(getproperty(measure, :t), "t", ctx)
    end

    return nothing
end

"""
    integrate(expr::SymbolicMatrix, measure)

Integrate a SymbolicMatrix as a whole. Returns a matrix of results if the
matrix dimensions are concrete integers. Supports rectangular matrices
(e.g., pure states with size `(d, 1)` or Stiefel matrices with size `(d, k)`).
"""
function integrate(A::SymbolicMatrix, measure::AbstractMeasure)
    rows_raw, cols_raw = size(A)
    nr = _try_extract_int(rows_raw)
    nc = _try_extract_int(cols_raw)
    if nr !== nothing && nc !== nothing
        res = Matrix{Any}(undef, nr, nc)
        fill!(res, 0)
        p = Progress(nr * nc; dt = 10.0, desc = "Integrating matrix elements... ")
        for i = 1:nr
            for j = 1:nc
                res[i, j] = integrate(A[i, j], measure)
                next!(p)
            end
        end
        return res
    end
    throw(
        ArgumentError(
            "Direct integration of SymbolicMatrix requires concrete integer dimensions.",
        ),
    )
end

_get_integration_tag(m::MetadataMatcher) = m.type_tag

function _has_integration_variable(expr, tag::Symbol)
    if expr isa SymbolicMatrix
        return expr.special_type === tag
    elseif expr isa SymbolicKron
        return _has_integration_variable(expr.A, tag) ||
               _has_integration_variable(expr.B, tag)
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

    _, matcher, _, _ = measure_info(measure)

    if matcher isa MetadataMatcher ||
       (isdefined(Main, :SymbolicMatcher) && matcher isa Main.SymbolicMatcher) ||
       hasproperty(matcher, :tag) ||
       hasproperty(matcher, :type_tag)
        tag = hasproperty(matcher, :type_tag) ? matcher.type_tag : matcher.tag
        if !_has_integration_variable(P, tag)
            size(P)  # validate internal dimension consistency
            return P
        end
    end

    dim_measure = _get_measure_dim(measure)
    n_factors = length(P.factors)

    inner_dims = Vector{Any}(fill(nothing, n_factors + 1))

    for (idx, f) in enumerate(P.factors)
        fr, fc = size(f)
        fr_int = _try_extract_int(fr)
        if fr_int === nothing
            fr_v = Symbolics.unwrap(fr)
            fr_int = fr_v isa Integer && fr_v != typemax(Int) ? fr_v : nothing
        end
        fc_int = _try_extract_int(fc)
        if fc_int === nothing
            fc_v = Symbolics.unwrap(fc)
            fc_int = fc_v isa Integer && fc_v != typemax(Int) ? fc_v : nothing
        end

        if fr_int !== nothing
            prev = inner_dims[idx]
            if prev isa Integer && prev != fr_int
                throw(ArgumentError(
                    "Dimension mismatch in matrix product: factor $idx has $fr_int rows " *
                    "but preceding factor has $prev columns."
                ))
            end
            inner_dims[idx] = fr_int
        end
        if fc_int !== nothing
            nxt = inner_dims[idx+1]
            if nxt isa Integer && nxt != fc_int
                throw(ArgumentError(
                    "Dimension mismatch in matrix product: factor $idx has $fc_int columns " *
                    "but following factor has $nxt rows."
                ))
            end
            inner_dims[idx+1] = fc_int
        end
    end

    for i = 1:(n_factors+1)
        if inner_dims[i] === nothing
            dim_int = _try_extract_int(dim_measure)
            inner_dims[i] = dim_int !== nothing ? dim_int : dim_measure
        end
    end

    nr = inner_dims[1] isa Integer ? inner_dims[1] : _try_extract_int(inner_dims[1])
    nc = inner_dims[end] isa Integer ? inner_dims[end] : _try_extract_int(inner_dims[end])

    if nr isa Integer && nc isa Integer
        mats = []
        for (idx, f) in enumerate(P.factors)
            cur_r = inner_dims[idx]
            cur_c = inner_dims[idx+1]

            if !(cur_r isa Integer && cur_c isa Integer)
                throw(
                    ArgumentError(
                        "Factor $f has non-numeric size ($cur_r, $cur_c). Cannot expand product for matrix-valued integration. Use tr() for scalar results with symbolic dimensions.",
                    ),
                )
            end

            push!(mats, [f[r, c] for r = 1:Int(cur_r), c = 1:Int(cur_c)])
        end
        res_mat = reduce(*, mats)

        res = Matrix{Any}(undef, nr, nc)
        fill!(res, 0)
        p = Progress(nr*nc; dt = 10.0, desc = "Integrating matrix product elements... ")
        for i = 1:nr
            for j = 1:nc
                res[i, j] = integrate(res_mat[i, j], measure)
                next!(p)
            end
        end
        return res
    end
    throw(
        ArgumentError(
            "Direct matrix-valued integration of SymbolicMatrixProduct requires numeric result dimensions (got $nr x $nc). Use tr() for scalar results.",
        ),
    )
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
    _validate_measure_discrete_params(measure)

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
    re = Symbolics.real(expr)
    im_part = Symbolics.imag(expr)

    int_re = _integrate_core(re, dim, subs_dict, matcher, measure_type)
    int_im = _integrate_core(im_part, dim, subs_dict, matcher, measure_type)

    return int_re + 1im * int_im
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
    throw(
        ArgumentError(
            "Fallback integrate not implemented for this measure: $(typeof(measure))",
        ),
    )
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
    dim = _assert_no_float_param(dim, "dim", "measure $tag")
    return (subs_dict, matcher, dim, tag)
end
