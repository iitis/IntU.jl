_to_Num(z::Complex) = Complex(Num(real(z)), Num(imag(z)))
_to_Num(z) = Num(z)

"""
    AbstractMeasure

Abstract base type for all integration measures (Haar, Gaussian, Circular, etc.).
Provides generic dispatch for element-wise integration of arrays and 
ambiguity resolution for `SymbolicMatrix` and `SymbolicMatrixProduct`.
"""
abstract type AbstractMeasure end

function _symbolic_isequal(a, b)
    a_v = Symbolics.value(a)
    b_v = Symbolics.value(b)

    if a_v isa Number && b_v isa Number
        return a_v == b_v
    end

    a_un = Symbolics.unwrap(a_v)
    b_un = Symbolics.unwrap(b_v)

    if a_un isa Complex && b_un isa Complex
        return _symbolic_isequal(real(a_un), real(b_un)) &&
               _symbolic_isequal(imag(a_un), imag(b_un))
    end

    if a_un isa Complex
        return _symbolic_isequal(real(a_un), b_un) && _iszero(imag(a_un))
    end

    if b_un isa Complex
        return _symbolic_isequal(a_un, real(b_un)) && _iszero(imag(b_un))
    end

    if Symbolics.iscall(a_un) &&
       (Symbolics.operation(a_un) == complex || Symbolics.operation(a_un) == Base.complex)
        args = Symbolics.arguments(a_un)
        return _symbolic_isequal(args[1], b_un) && _iszero(args[2])
    end

    if Symbolics.iscall(b_un) &&
       (Symbolics.operation(b_un) == complex || Symbolics.operation(b_un) == Base.complex)
        args = Symbolics.arguments(b_un)
        return _symbolic_isequal(a_un, args[1]) && _iszero(args[2])
    end

    res = isequal(a_un, b_un)
    v = Symbolics.value(res)
    return v === true
end

function _iszero(x)
    v = Symbolics.value(x)
    if v isa Number
        return iszero(v)
    end
    s = Symbolics.simplify(x)
    sv = Symbolics.value(s)
    if sv isa Number
        return iszero(sv)
    end
    return _symbolic_isequal(v, 0)
end

"""
    _ensure_symbolic_dim(d)

Ensure dimension `d` is a proper Symbolics variable. If `d` unwraps to a plain
`Symbol`, wrap it via `Symbolics.variable`; otherwise return as-is.
"""
function _ensure_symbolic_dim(d)
    d_un = Symbolics.unwrap(d)
    return d_un isa Symbol ? Symbolics.variable(d_un) : d
end

"""
    _try_numeric(v)

Attempt to convert a value to a clean numeric form. Returns the converted value,
or `nothing` if conversion is not possible.
- `AbstractFloat` → rationalized
- `Real` → returned as-is
"""
function _try_numeric(v)
    if v isa AbstractFloat
        return rationalize(v, tol = 1e-13)
    end
    if v isa Real
        return v
    end
    return nothing
end

function _robust_real(x)
    if x isa AbstractArray
        return map(_robust_real, x)
    end

    x_un = Symbolics.unwrap(x)

    v = Symbolics.value(x_un)
    result = _try_numeric(v)
    result !== nothing && return result

    if v isa Complex
        rv = _robust_real(real(v))
        iv = _robust_real(imag(v))
        if iszero(iv)
            return rv
        end
        return Complex(rv, iv)
    end

    if Symbolics.iscall(x_un) &&
       (Symbolics.operation(x_un) == complex || Symbolics.operation(x_un) == Base.complex)
        args = Symbolics.arguments(x_un)
        if _iszero(args[2])
            return _robust_real(args[1])
        end

    end


    if x_un isa Real
        return x_un
    end

    if _is_manifestly_real(x_un)
        return x_un
    end


    nx = _safe_Num(x_un)

    if !(nx isa Num || nx isa Complex{Num})
        return x_un
    end

    v = Symbolics.value(nx)
    result = _try_numeric(v)
    result !== nothing && return result

    nx = Symbolics.simplify(nx)
    v = Symbolics.value(nx)
    result = _try_numeric(v)
    result !== nothing && return result

    if _iszero(Symbolics.simplify(imag(nx)))
        rx = Symbolics.simplify(real(nx))
        vx = Symbolics.value(rx)
        result = _try_numeric(vx)
        result !== nothing && return result
        return rx
    end
    return nx
end

function _is_manifestly_real(x)
    if x isa Real
        return true
    end
    if x isa AbstractArray

        for elem in x
            u = Symbolics.unwrap(elem)
            # Prevent infinite recursion for symbolic arrays with getindex(self, ...)
            if Symbolics.iscall(u) && Symbolics.operation(u) == getindex
                args = Symbolics.arguments(u)
                if args[1] === x || args[1] == x
                    continue
                end
            end

            if !_is_manifestly_real(u)
                return false
            end
        end
        return true
    end
    if x isa Complex
        return iszero(imag(x))
    end
    if SymbolicUtils.issym(x)
        return true
    end

    if SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        if op == complex || op == Base.complex || op == imag || op == Base.imag
            return false
        end

        args = SymbolicUtils.arguments(x)
        for arg in args
            if !_is_manifestly_real(arg)
                return false
            end
        end
        return true
    end

    v = Symbolics.value(x)
    if v !== x
        return _is_manifestly_real(v)
    end

    return false
end

"""
    AbstractIndexMatcher

Base type for index matching strategies. Subtypes must implement `match_index`.
"""
abstract type AbstractIndexMatcher end



"""
    MetadataMatcher(type_tag::Symbol)

A matcher that identifies random matrix entries based on metadata attached to the symbols.
The `type_tag` should match the `special_type` of a `SymbolicMatrix` (e.g., `:U`, `:O`, `:Sp`).
"""
struct MetadataMatcher <: AbstractIndexMatcher
    type_tag::Symbol # :U, :O, :Sp, etc.
end

function match_index(m::MetadataMatcher, t)
    s = Symbolics.unwrap(t)

    if s isa SymbolicMatrix
        if s.special_type === m.type_tag
            final_tag = s.is_adj ? Symbol(m.type_tag, :_bar) : m.type_tag
            return (final_tag, nothing, nothing)
        end
        return nothing
    end

    # Handle conj(U_i_j)
    is_conj = false
    if Symbolics.iscall(s) &&
       (Symbolics.operation(s) == conj || Symbolics.operation(s) == Base.conj)
        is_conj = true
        s = Symbolics.arguments(s)[1]
    end

    can_have_meta =
        (s isa SymbolicUtils.BasicSymbolic) && SymbolicUtils.hasmetadata(s, MatrixMetadata)

    if can_have_meta
        meta = SymbolicUtils.getmetadata(s, MatrixMetadata)
        if get(meta, :type, nothing) === m.type_tag
            indices = get(meta, :indices, nothing)
            if indices !== nothing
                i, j = indices
                # Combine is_conj (from call) and :is_adj (from metadata)
                final_is_conj = is_conj || get(meta, :is_adj, false)

                final_tag = final_is_conj ? Symbol(m.type_tag, :_bar) : m.type_tag
                return (final_tag, i, j)
            end
        end
    end
    return nothing
end

"""
    INTEGRATION_RULES

A dictionary mapping measure types (symbols) to their respective integration rule functions.
Each rule function should have the signature `(u_indices, u_bar_indices, dim, measure_type)`.
"""
const INTEGRATION_RULES = Dict{Symbol,Function}()

function _extract_coeff_core(term)
    if term isa Number
        return term, 1
    end

    if Symbolics.iscall(term) && Symbolics.operation(term) == (*)
        args = Symbolics.arguments(term)
        c = 1
        others = Any[]
        for a in args
            if a isa Number
                c *= a
            else
                push!(others, a)
            end
        end

        if isempty(others)
            core = 1
        elseif length(others) == 1
            core = others[1]
        else
            core = prod(others)
        end
        return c, core
    else
        return 1, term
    end
end

function _is_fn_sq(term, fn1, fn2)
    # Strip potential 1 * ... wrapper
    if Symbolics.iscall(term) && Symbolics.operation(term) == (*)
        args = Symbolics.arguments(term)
        if length(args) == 2 && isequal(args[1], 1)
            term = args[2]
        end
    end

    if Symbolics.iscall(term) && Symbolics.operation(term) == (^)
        args = Symbolics.arguments(term)
        base = args[1]
        expon = args[2]
        if isequal(expon, 2) &&
           Symbolics.iscall(base) &&
           (Symbolics.operation(base) == fn1 || Symbolics.operation(base) == fn2)
            return true, Symbolics.arguments(base)[1]
        end
    end
    return false, nothing
end

_is_real_sq(term) = _is_fn_sq(term, real, Base.real)
_is_imag_sq(term) = _is_fn_sq(term, imag, Base.imag)

"""
    robust_substitute(ex, dict)

Substitute symbolic variables in `ex` using `dict`. Falls back to manual
Postwalk traversal with unwrapped keys if `Symbolics.substitute` doesn't
produce a different expression.
"""
function robust_substitute(ex, dict)
    res = Symbolics.substitute(Symbolics.wrap(ex), dict)

    if !_symbolic_isequal(res, Symbolics.wrap(ex))
        return res
    end

    # Fallback: manual Postwalk traversal with unwrapped keys
    unwrapped_dict = Dict(Symbolics.unwrap(k) => Symbolics.unwrap(v) for (k, v) in dict)

    p_res = SymbolicUtils.Postwalk(x -> begin
        u = Symbolics.unwrap(x)
        if haskey(unwrapped_dict, u)
            return unwrapped_dict[u]
        end
        return x
    end)(
        Symbolics.unwrap(ex),
    )

    return Symbolics.wrap(p_res)
end

"""
    _integrate_core(expr, dim, subs_dict, matcher, measure_type=:U)

The internal integration engine. It performs several steps:
1.  **Normalization/Rewriting**: Expands `abs2(z)`, `real(z)`, `imag(z)` into explicit polynomials.
2.  **Substitution**: Replaces symbolic variables with internal atomic representatives via `subs_dict`.
3.  **Expansion**: Distributes products over sums to get a sum of monomials.
4.  **Monomial Integration**: For each monomial, invokes `process_term` to identify indices and apply the appropriate integration rule (Weingarten or Wick).
"""
function _integrate_core(
    expr,
    dim,
    subs_dict,
    matcher::AbstractIndexMatcher,
    measure_type = :U,
)
    if expr isa Complex
        val_re = _integrate_core(real(expr), dim, subs_dict, matcher, measure_type)
        val_im = _integrate_core(imag(expr), dim, subs_dict, matcher, measure_type)
        return _robust_real(val_re + 1im * val_im)
    end

    expr_un = Symbolics.unwrap(expr)
    if expr_un isa Number && (expr_un isa Real || expr_un isa Complex)
        return expr_un
    end


    if expr_un isa Complex
        return _integrate_core(expr_un, dim, subs_dict, matcher, measure_type)
    end


    function is_add(t)
        Symbolics.iscall(t) && Symbolics.operation(t) == (+)
    end

    r_abs2 = @rule abs2(~x) => (~x) * conj(~x)
    r_abs_pow = @rule abs(~x)^~n => begin
        n_un = Symbolics.unwrap(~n)
        if n_un isa Number && isinteger(n_un) && iseven(Int(n_un))
            k = Int(n_un) ÷ 2
            ((~x) * conj(~x))^k
        else
            nothing
        end
    end
    r_abs = @rule abs(~x) => hypot(real(~x), imag(~x))

    r_real = @rule real(~x) => (1//2) * (~x + conj(~x))
    r_real_base = @rule Base.real(~x) => (1//2) * (~x + conj(~x))
    r_imag = @rule imag(~x) => (1//(2im)) * (~x - conj(~x))
    r_imag_base = @rule Base.imag(~x) => (1//(2im)) * (~x - conj(~x))


    r_hypot_sum =
        @rule hypot(~x::is_add, ~y::is_add) => ((~x + im*~y) * (~x - im*~y))^(1//2)

    r_hypot_default = @rule hypot(~x, ~y) => ((~x)^2 + (~y)^2)^(1//2)


    r_complex = @rule complex(~x, ~y) => ~x + im*~y
    r_complex_base = @rule Base.complex(~x, ~y) => ~x + im*~y

    r_rev_abs2 = @rule real(~x)^2 + imag(~x)^2 => (~x) * conj(~x)
    r_rev_abs2_coeff = @rule ~c * real(~x)^2 + ~c * imag(~x)^2 => ~c * (~x) * conj(~x)

    r_rev_abs2_base = @rule Base.real(~x)^2 + Base.imag(~x)^2 => (~x) * conj(~x)
    r_rev_abs2_coeff_base =
        @rule ~c * Base.real(~x)^2 + ~c * Base.imag(~x)^2 => ~c * (~x) * conj(~x)

    r_conj_add = @rule conj(~x + ~y) => conj(~x) + conj(~y)
    r_conj_add_base = @rule Base.conj(~x + ~y) => conj(~x) + conj(~y)
    r_conj_mul = @rule conj(~x * ~y) => conj(~x) * conj(~y)
    r_conj_mul_base = @rule Base.conj(~x * ~y) => conj(~x) * conj(~y)
    r_conj_pow = @rule conj((~x)^~n) => (conj(~x))^~n
    r_conj_conj = @rule conj(conj(~x)) => ~x
    r_conj_conj_base = @rule Base.conj(Base.conj(~x)) => ~x
    r_conj_neg = @rule conj(-(~x)) => -(conj(~x))


    function power_simplifier(x, a, b)
        new_expon = a * b
        if isinteger(new_expon)
            return x^Int(new_expon)
        end
        return x^new_expon
    end
    r_pow_nested = @rule ((~x)^~a)^~b => power_simplifier(~x, ~a, ~b)


    r_hypot_pow = @rule hypot(~x, ~y)^~n => begin
        n_un = Symbolics.unwrap(~n)
        if n_un isa Number && isinteger(n_un) && iseven(Int(n_un))
            k = Int(n_un) ÷ 2
            ((~x)^2 + (~y)^2)^k
        else
            nothing
        end
    end


    r_float_to_int_pow = @rule (~x)^~a => begin
        a_un = Symbolics.unwrap(~a)
        if a_un isa AbstractFloat && isinteger(a_un)
            (~x)^Int(a_un)
        else
            nothing
        end
    end

    chain = SymbolicUtils.Chain([
        r_rev_abs2_coeff,
        r_rev_abs2,
        r_rev_abs2_coeff_base,
        r_rev_abs2_base,
        r_abs_pow,
        r_abs2,
        r_abs,
        r_real,
        r_real_base,
        r_imag,
        r_imag_base,
        r_hypot_sum,
        r_hypot_default,
        r_complex,
        r_complex_base,
        r_conj_add,
        r_conj_add_base,
        r_conj_mul,
        r_conj_mul_base,
        r_conj_pow,
        r_conj_conj,
        r_conj_conj_base,
        r_conj_neg,
        r_pow_nested,
        r_hypot_pow,
        r_float_to_int_pow,
    ])

    function pair_real_imag(expr)
        if !Symbolics.iscall(expr) || Symbolics.operation(expr) != (+)
            return expr
        end

        args = Symbolics.arguments(expr)
        n_args = length(args)
        new_args = []
        skip_indices = falses(n_args)

        for i = 1:n_args
            if skip_indices[i]
                continue
            end

            term_i = args[i]

            c_i, core_i = _extract_coeff_core(term_i)

            is_real_sq, x_real = _is_real_sq(core_i)
            is_imag_sq, x_imag = _is_imag_sq(core_i)

            matched = false

            if is_real_sq || is_imag_sq
                target_x = is_real_sq ? x_real : x_imag
                target_is_real = !is_real_sq

                for j = (i+1):n_args
                    if skip_indices[j]
                        continue
                    end

                    term_j = args[j]
                    c_j, core_j = _extract_coeff_core(term_j)

                    if !isequal(c_i, c_j)
                        continue
                    end

                    is_real_sq_j, x_real_j = _is_real_sq(core_j)
                    is_imag_sq_j, x_imag_j = _is_imag_sq(core_j)

                    if target_is_real
                        if is_real_sq_j && isequal(x_real_j, target_x)
                            # Found match!
                            push!(new_args, c_i * target_x * conj(target_x))
                            skip_indices[j] = true
                            matched = true
                            break
                        end
                    else
                        if is_imag_sq_j && isequal(x_imag_j, target_x)
                            # Found match!
                            push!(new_args, c_i * target_x * conj(target_x))
                            skip_indices[j] = true
                            matched = true
                            break
                        end
                    end
                end
            end

            if !matched
                push!(new_args, term_i)
            end
        end

        return sum(new_args)
    end

    function monomializer(x)
        ux = Symbolics.unwrap(x)
        if !Symbolics.iscall(ux)
            return x
        end
        op = Symbolics.operation(ux)
        args = Symbolics.arguments(ux)
        if op == conj || op == Base.conj
            inner = Symbolics.unwrap(args[1])
            if Symbolics.iscall(inner)
                inner_op = Symbolics.operation(inner)
                inner_args = Symbolics.arguments(inner)
                if inner_op == (+)
                    return sum(monomializer(conj(a)) for a in inner_args)
                elseif inner_op == (*)
                    return prod(monomializer(conj(a)) for a in inner_args)
                elseif inner_op == (/)
                    return monomializer(conj(inner_args[1])) /
                           monomializer(conj(inner_args[2]))
                elseif inner_op == (^)
                    return monomializer(conj(inner_args[1])) ^
                           monomializer(conj(inner_args[2]))
                elseif inner_op == conj || inner_op == Base.conj
                    return monomializer(Symbolics.arguments(inner)[1])
                end
            end
        elseif op == (/)
            # (A + B) / C -> A/C + B/C
            numerator = Symbolics.unwrap(args[1])
            if Symbolics.iscall(numerator) && Symbolics.operation(numerator) == (+)
                return sum(
                    monomializer(a / args[2]) for a in Symbolics.arguments(numerator)
                )
            end
        elseif op == complex || op == Base.complex
            return monomializer(args[1]) + im * monomializer(args[2])
        elseif op == real || op == Base.real
            return monomializer((args[1] + conj(args[1])) / 2)
        elseif op == imag || op == Base.imag
            return monomializer((args[1] - conj(args[1])) / (2im))
        end
        return Symbolics.wrap(
            Symbolics.iscall(ux) ?
            SymbolicUtils.maketerm(
                typeof(ux),
                op,
                [monomializer(a) for a in args],
                SymbolicUtils.metadata(ux),
            ) : ux,
        )
    end

    expr_rewritten = SymbolicUtils.Postwalk(
        x -> begin
            # Apply rules first
            res = SymbolicUtils.PassThrough(chain)(x)
            # Then distribute
            return monomializer(res)
        end;
        maketerm = (st, f, args, metadata; kwargs...) -> begin
            if f == complex || f == Base.complex
                return args[1] + im * args[2]
            end
            SymbolicUtils.maketerm(st, f, args, metadata)
        end,
    )(
        expr_un,
    )

    expr_rewritten = pair_real_imag(expr_rewritten)

    expr_num = _safe_Num(expr_rewritten)
    expr_num = Symbolics.expand(expr_num)

    expr_subbed = if expr_num isa Complex
        robust_substitute(Symbolics.unwrap(real(expr_num)), subs_dict) +
        im * robust_substitute(Symbolics.unwrap(imag(expr_num)), subs_dict)
    else
        robust_substitute(Symbolics.unwrap(expr_num), subs_dict)
    end



    expanded_expr = Symbolics.expand(_safe_Num(Symbolics.unwrap(expr_subbed)))

    # Flatten complex(...) calls introduced by substitution
    expanded_expr = _safe_Num(
        SymbolicUtils.Postwalk(
            x -> begin
                ux = Symbolics.unwrap(x)
                if Symbolics.iscall(ux) && (
                    Symbolics.operation(ux) == complex ||
                    Symbolics.operation(ux) == Base.complex
                )
                    args = Symbolics.arguments(ux)
                    return args[1] + im * args[2]
                end
                return x
            end,
        )(
            Symbolics.unwrap(expanded_expr),
        ),
    )


    expanded_expr = Symbolics.expand(expanded_expr)


    function process_term_wrapped(term)
        return process_term(term, matcher, dim, measure_type)
    end


    function integrate_num_expr(ex)
        ex_un = Symbolics.unwrap(ex)
        if ex_un isa Complex
            return integrate_num_expr(real(ex_un)) + im * integrate_num_expr(imag(ex_un))
        end
        if ex_un isa Number && (ex_un isa Real)
            return ex_un
        end
        if Symbolics.iscall(ex_un)
            op = Symbolics.operation(ex_un)
            if op == (+)
                terms = Symbolics.arguments(ex_un)
                n_terms = length(terms)
                if n_terms > 1
                    p = Progress(n_terms; dt=10.0, desc="Integrating terms... ")
                    res = zero(Num)
                    for t in terms
                        res += process_term_wrapped(t)
                        next!(p)
                    end
                    return res
                end
                return sum(process_term_wrapped, terms)
            end
            if op == complex || op == Base.complex
                args = Symbolics.arguments(ex_un)
                return integrate_num_expr(args[1]) + im * integrate_num_expr(args[2])
            end
        end
        return process_term_wrapped(ex)
    end


    final_res = integrate_num_expr(expanded_expr)
    res = _robust_real(final_res)
    return res
end

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
    measure_info(measure)

Returns `(subs_dict, matcher, dim, measure_type)` for a given measure.
Subtypes should implement this to participate in the unified integration flow.
"""
function measure_info(measure)
    return nothing
end

function _safe_Num(x)
    if x isa Num || x isa Complex{Num}
        return x
    end
    if x isa LazyTrace || x isa LazySum || x isa LazyPower
        return x
    end
    if x isa AbstractArray
        return map(_safe_Num, x)
    end
    if !(x isa Number)
        # Only wrap if it's a known symbolic object to avoid overhead/errors on obscure types
        if x isa SymbolicUtils.BasicSymbolic || x isa Symbolics.ComplexTerm
            return Num(x)
        end
        return x
    end
    x_un = Symbolics.unwrap(x)
    return _to_Num(x_un)
end

"""
    _try_extract_int(p_val)

Attempt to extract an integer from a symbolic or numeric value.
Handles Julia `Integer`, `AbstractFloat` with integer value, and
`SymbolicUtils.BasicSymbolic` types that print as integers.
Returns the `Int` value or `nothing`.
"""
function _try_extract_int(p_val)
    if p_val isa Integer
        return Int(p_val)
    end
    if p_val isa AbstractFloat && isinteger(p_val)
        return Int(p_val)
    end
    # Fallback for symbolic types (e.g. BasicSymbolic{SymReal}) that
    # represent integers but aren't Julia Integer subtypes
    return tryparse(Int, string(p_val))
end


function _robust_real_num(x)
    res = _robust_real(x)
    return _safe_Num(res)
end


"""
    process_term(term, matcher, dim, measure_type)

Integrates a single monomial term.
1.  **Index Collection**: Traverses the term to find all random matrix elements using the `matcher`.
2.  **Dispatch**: Calls specific index integration functions based on `measure_type`:
    - `:U`: `integrate_indices` (Weingarten)
    - `:O`: `integrate_indices_orthogonal`
    - `:Sp`: `integrate_indices_symplectic`
    - `:GUE`, `:GOE`, `:GSE`: Wick contraction rules.
"""
function process_term(term, matcher::AbstractIndexMatcher, dim, measure_type = :U)
    term = Symbolics.unwrap(term)
    if term isa LazyTrace
        # If it leaked here, it means it's not specialized for this measure.
        # Try a very basic expansion? No, let's error gracefully if not handled.
        error(
            "Graphical integration (LazyTrace) not implemented for measure type $measure_type. Try expanding traces element-wise.",
        )
    end

    if Symbolics.iscall(term)
        op = Symbolics.operation(term)
        args = Symbolics.arguments(term)
        if op == (+)
            return sum(t -> process_term(t, matcher, dim, measure_type), args)
        elseif op == complex || op == Base.complex
            return process_term(args[1], matcher, dim, measure_type) +
                   im * process_term(args[2], matcher, dim, measure_type)
        end
    end

    is_gaussian = (measure_type === :GUE || measure_type === :GOE || measure_type === :GSE)
    coeff = 1 // 1
    u_indices = Vector{Tuple{Any,Any}}()
    u_bar_indices = Vector{Tuple{Any,Any}}()

    function _push_matched_index!(match_res, conjugated)
        type, i, j = match_res
        is_bar = endswith(string(type), "_bar")
        final_is_bar = is_bar ⊻ conjugated
        if final_is_bar
            push!(u_bar_indices, (i, j))
        else
            push!(u_indices, (i, j))
        end
    end

    function traverse(t, conjugated = false)
        t_unwrapped = Symbolics.unwrap(t)

        # Try to match the variable directly (handles getindex, metadata, etc.)
        match_res = match_index(matcher, t_unwrapped)
        if match_res !== nothing
            _push_matched_index!(match_res, conjugated)
            return
        end

        if t isa Number && !(t isa Num) && !(t isa Complex{Num})
            coeff *= t
            return
        end

        if Symbolics.iscall(t_unwrapped)
            local op = Symbolics.operation(t_unwrapped)
            local args = Symbolics.arguments(t_unwrapped)

            # Multiplicative: traverse each factor
            if op == (*)
                for arg in args
                    traverse(arg, conjugated)
                end
                return
            # Power: repeat base p times
            elseif op == (^)
                base = args[1]
                p_val = Symbolics.unwrap(args[2])
                p = _try_extract_int(p_val)
                if p isa Integer
                    for _ = 1:p
                        traverse(base, conjugated)
                    end
                    return
                end
                traverse(base, conjugated)
                return
            # Division: traverse numerator, divide coefficient
            elseif op == (/)
                traverse(args[1], conjugated)
                coeff /= args[2]
                return
            # Conjugation: flip conjugated flag
            elseif op == conj || op == Base.conj
                traverse(args[1], !conjugated)
                return
            # real(x) = (x + conj(x)) / 2
            elseif op == real || op == Base.real
                coeff *= 1 // 2
                traverse(args[1], conjugated)
                traverse(args[1], !conjugated)
                return
            # imag(x) = (x - conj(x)) / (2im)
            elseif op == imag || op == Base.imag
                coeff *= 1 // (2im)
                traverse(args[1], conjugated)
                coeff *= -1
                traverse(args[1], !conjugated)
                return
            # Distributive: complex, +, -
            elseif op == complex || op == Base.complex || op == (+) || op == (-)
                for arg in args
                    traverse(arg, conjugated)
                end
                return
            end
        end

        coeff *= conjugated ? conj(t) : t
    end
    traverse(term, false)

    n_u = length(u_indices)
    n_bar = length(u_bar_indices)

    rule = get(INTEGRATION_RULES, measure_type, nothing)
    if rule === nothing && measure_type isa Tuple
        rule = get(INTEGRATION_RULES, first(measure_type), nothing)
    end

    if rule !== nothing
        val = rule(u_indices, u_bar_indices, dim, measure_type)
        if _symbolic_isequal(val, 0)
            return 0
        end
        return coeff * val
    else
        error("Unknown measure type: $measure_type")
    end
end

INTEGRATION_RULES[:U] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        length(u) != length(ub) ? 0 : (length(u) == 0 ? 1 : integrate_indices(u, ub, d))
    end

# Stiefel V_k(C^d) and pure state integration are same as Haar U(d) for entries
INTEGRATION_RULES[:V] = INTEGRATION_RULES[:U]
INTEGRATION_RULES[:psi] = INTEGRATION_RULES[:U]

INTEGRATION_RULES[:O] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_orthogonal(all_indices, d))
    end

INTEGRATION_RULES[:Sp] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        d_un = Symbolics.unwrap(d)
        if d_un isa Integer
            if isodd(d_un)
                # Sp(d) only exists for even d
                return 0
            end
            m = div(d_un, 2)
            sign_factor = 1
            converted_ub = Vector{Tuple{Any,Any}}()

            for (p, q) in ub
                pk = p <= m ? p + m : p - m
                qk = q <= m ? q + m : q - m
                # \bar{S}_{p,q} = (-1)^{p > m ? 1 : 0} * (-1)^{q > m ? 1 : 0} * S_{pk, qk}
                # sign = (-1)^( (p>m) + (q>m) )
                if (p <= m && q > m) || (p > m && q <= m)
                    sign_factor *= -1
                end
                push!(converted_ub, (pk, qk))
            end

            all_indices = [u; converted_ub]
            length(all_indices) % 2 != 0 ? 0 :
            (length(all_indices) == 0 ? 1 : sign_factor * integrate_indices_symplectic(all_indices, d))
        else
            # Symbolic d. Assume d is even.
            # Convert conjugates using symbolic partner indices pk = p + d/2, qk = q + d/2.
            # This is valid if p, q are "small" compared to d/2.
            m = d / 2
            sign_factor = 1
            converted_ub = Vector{Tuple{Any,Any}}()

            for (p, q) in ub
                pk = p + m
                qk = q + m
                push!(converted_ub, (pk, qk))
            end

            all_indices = [u; converted_ub]
            length(all_indices) % 2 != 0 ? 0 :
            (length(all_indices) == 0 ? 1 : sign_factor * integrate_indices_symplectic(all_indices, d))
        end
    end

INTEGRATION_RULES[:GUE] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_gue(all_indices, d))
    end

INTEGRATION_RULES[:GOE] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_goe(all_indices, d))
    end

INTEGRATION_RULES[:GSE] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_gse(all_indices, d))
    end

INTEGRATION_RULES[:COE] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        length(u) != length(ub) ? 0 : (length(u) == 0 ? 1 : integrate_indices_coe(u, ub, d))
    end

INTEGRATION_RULES[:CSE] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        d_un = Symbolics.unwrap(d)
        length(u) != length(ub) ? 0 :
        (
            begin
                phys_dim = mt isa Tuple ? mt[2] : d
                length(u) == 0 ? 1 : integrate_indices_cse(u, ub, d, phys_dim)
            end
        )
    end

INTEGRATION_RULES[:Perm] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        all_indices = [u; ub]
        length(all_indices) == 0 ? 1 : integrate_indices_permutation(all_indices, d)
    end

INTEGRATION_RULES[:CPerm] =
    (u, ub, d, mt) -> begin
        d = _ensure_symbolic_dim(d)
        all_indices = [u; ub]
        length(all_indices) == 0 ? 1 :
        integrate_indices_centered_permutation(all_indices, d)
    end

INTEGRATION_RULES[:Design] =
    (u, ub, d, mt) -> begin
        _, t_val = mt
        if length(u) != length(ub) || length(u) > t_val
            length(u) != length(ub) ? 0 :
            error("Integrand degree ($(length(u))) exceeds design order t=$t_val")
        end
        length(u) == 0 ? 1 : integrate_indices(u, ub, d)
    end

INTEGRATION_RULES[:DiagUnitary] =
    (u, ub, d, mt) -> begin
        if length(u) != length(ub) || any(x -> x[1] != x[2], u) || any(x -> x[1] != x[2], ub)
            return 0
        end
        length(u) == 0 ? 1 : integrate_indices_diagonal(u, ub, d)
    end

INTEGRATION_RULES[:GinUE] =
    (u, ub, d, mt) -> begin
        length(u) != length(ub) ? 0 : (length(u) == 0 ? 1 : integrate_indices_ginue(u, ub, d))
    end

INTEGRATION_RULES[:GinOE] =
    (u, ub, d, mt) -> begin
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_ginoe(all_indices, d))
    end

INTEGRATION_RULES[:GinSE] =
    (u, ub, d, mt) -> begin
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_ginse(all_indices, d))
    end


"""
    integrate_indices(U_idxs, U_bar_idxs, dim)

Low-level integration function using Weingarten calculus (Unitary).
"""
function integrate_indices(U_idxs::Vector{<:Tuple}, U_bar_idxs::Vector{<:Tuple}, dim)
    n = length(U_idxs)
    I = [x[1] for x in U_idxs]
    J = [x[2] for x in U_idxs]
    I_bar = [x[1] for x in U_bar_idxs]
    J_bar = [x[2] for x in U_bar_idxs]

    valid = get_matching_permutations(I, I_bar)
    valid_taus = get_matching_permutations(J, J_bar)

    if isempty(valid) || isempty(valid_taus)
        return 0
    end

    # Group terms by cycle type
    cycle_counts = Dict{Vector{Int},Int}()
    n_v = length(valid)
    n_vt = length(valid_taus)
    if n_v * n_vt > 100 # Only if there are enough terms to justify overhead
        p = Progress(n_v * n_vt; dt=10.0, desc="Calculating cycle types... ")
        for sigma in valid
            for tau in valid_taus
                inv_tau = invperm(tau)
                P = [sigma[inv_tau[i]] for i = 1:n]
                cycle_type = get_cycle_type(P)
                cycle_counts[cycle_type] = get(cycle_counts, cycle_type, 0) + 1
                next!(p)
            end
        end
    else
        for sigma in valid
            for tau in valid_taus
                inv_tau = invperm(tau)
                P = [sigma[inv_tau[i]] for i = 1:n]
                cycle_type = get_cycle_type(P)
                cycle_counts[cycle_type] = get(cycle_counts, cycle_type, 0) + 1
            end
        end
    end

    if dim isa Integer
        total = zero(Rational{BigInt})
        for (ct, count) in cycle_counts
            total += count * weingarten(ct, dim)
        end
        return total
    else
        # Symbolic dim: accumulate numerator/denominator from weingarten_raw
        total_num = zero(Num)
        total_den = one(Num)

        for (ct, count) in cycle_counts
            wnum, wden = IntU.weingarten_raw(ct, dim)
            if isequal(total_den, 1)
                total_num = count * wnum
                total_den = wden
            else
                if isequal(total_den, wden)
                    total_num += count * wnum
                else
                    total_num = total_num * wden + count * wnum * total_den
                    total_den = total_den * wden
                end
            end
        end

        return total_num / total_den
    end
end

function get_matching_permutations(target::AbstractVector, source::AbstractVector)
    n = length(target)
    if n != length(source)
        return Vector{Vector{Int}}()
    end


    source_groups = Dict{Any,Vector{Int}}()
    for (idx, val) in enumerate(source)
        push!(get!(source_groups, val, Int[]), idx)
    end

    target_groups = Dict{Any,Vector{Int}}()
    for (idx, val) in enumerate(target)
        push!(get!(target_groups, val, Int[]), idx)
    end


    if length(source_groups) != length(target_groups)
        return Vector{Vector{Int}}()
    end
    for (val, t_idxs) in target_groups
        if !haskey(source_groups, val) || length(source_groups[val]) != length(t_idxs)
            return Vector{Vector{Int}}()
        end
    end


    vals = collect(keys(target_groups))
    group_perms = [collect(permutations(source_groups[v])) for v in vals]

    res = Vector{Vector{Int}}()
    for combined_p in Iterators.product(group_perms...)
        p = zeros(Int, n)
        for (v_idx, v_perms) in enumerate(combined_p)
            t_idxs = target_groups[vals[v_idx]]
            for (i, t_idx) in enumerate(t_idxs)
                p[t_idx] = v_perms[i]
            end
        end
        push!(res, p)
    end

    return res
end

function get_cycle_type(p::Vector{Int})
    n = length(p)
    if n > 64
        visited = falses(n)
        lengths = Int[]
        for i = 1:n
            if !visited[i]
                curr = i
                len = 0
                while !visited[curr]
                    visited[curr] = true
                    curr = p[curr]
                    len += 1
                end
                push!(lengths, len)
            end
        end
        sort!(lengths, rev = true)
        return lengths
    end

    visited = UInt64(0)
    lengths = Int[]
    sizehint!(lengths, n)
    for i = 1:n
        if (visited & (UInt64(1) << (i - 1))) == 0
            curr = i
            len = 0
            while (visited & (UInt64(1) << (curr - 1))) == 0
                visited |= (UInt64(1) << (curr - 1))
                curr = p[curr]
                len += 1
            end
            push!(lengths, len)
        end
    end
    sort!(lengths, rev = true)
    return lengths
end

"""
    integrate_indices_orthogonal(indices, dim)
    
Low-level integration function using Orthogonal Weingarten calculus.
Indices are a list of (i, j) for O_{ij}.
Formula: sum_{pi, sigma in PairPartitions} delta_pi(i) * delta_sigma(j) * Wg(pi, sigma)
"""
function integrate_indices_orthogonal(indices::AbstractVector, dim)
    n = length(indices) # This is 2k
    if n % 2 != 0
        return 0
    end
    k = n ÷ 2

    I = [x[1] for x in indices]
    J = [x[2] for x in indices]

    # Shortcut: if all I are equal OR all J are equal
    is_I_uniform = all(x -> x == I[1], I)
    is_J_uniform = all(x -> x == J[1], J)

    if is_I_uniform || is_J_uniform
        # Row sum of Wg matrix is 1 / prod_{j=0}^{k-1} (dim + 2j)
        # Total sum = count_valid_partitions * row_sum
        if is_I_uniform && is_J_uniform
            # Both uniform: result is N_total / prod
            num = 1
            for i = 1:k
                num *= (2*i-1)
            end
            denom = one(Rational{BigInt})
            for j = 0:(k-1)
                denom *= (dim + 2*j)
            end
            return num / denom
        elseif is_I_uniform
            valid_sigma = get_matching_pair_partitions_filtered(J)
            denom = one(Rational{BigInt})
            for j = 0:(k-1)
                denom *= (dim + 2*j)
            end
            return length(valid_sigma) / denom
        else # is_J_uniform
            valid_pi = get_matching_pair_partitions_filtered(I)
            denom = one(Rational{BigInt})
            for j = 0:(k-1)
                denom *= (dim + 2*j)
            end
            return length(valid_pi) / denom
        end
    end

    valid_pi = get_matching_pair_partitions_filtered(I)
    valid_sigma = get_matching_pair_partitions_filtered(J)

    if isempty(valid_pi) || isempty(valid_sigma)
        return 0
    end

    # Group terms by canonicalized pair partitions
    pi_counts = Dict{Any,Int}()
    for pi in valid_pi
        c_pi = canonicalize_pair_partition(pi)
        pi_counts[c_pi] = get(pi_counts, c_pi, 0) + 1
    end

    sigma_counts = Dict{Any,Int}()
    for sigma in valid_sigma
        c_sigma = canonicalize_pair_partition(sigma)
        sigma_counts[c_sigma] = get(sigma_counts, c_sigma, 0) + 1
    end

    # Pre-aggregate counts by cycle type to minimize symbolic additions
    type_counts = Dict{Vector{Int}, BigInt}()
    for (c_pi, count_pi) in pi_counts
        for (c_sigma, count_sigma) in sigma_counts
            ct = get_full_cycle_type(c_pi, c_sigma)
            type_counts[ct] = get(type_counts, ct, zero(BigInt)) + BigInt(count_pi) * BigInt(count_sigma)
        end
    end

    if !(dim isa Integer)
        w, type_to_idx, _ = get_weingarten_reduced_data(k, dim, return_rationals=true)
        
        # Exact rational summation to bypass Symbolics.simplify bugs
        local_total_rat = nothing
        for (ct, count) in type_counts
            term = _rational_mul(w[type_to_idx[ct]], count)
            if local_total_rat === nothing
                local_total_rat = term
            else
                local_total_rat = _rational_add(local_total_rat, term)
            end
        end
        
        if local_total_rat === nothing
            return Num(0)
        end
        return from_rational(local_total_rat)
    end
    
    # Numeric case
    w, type_to_idx, _ = get_weingarten_reduced_data(k, dim)
    T_val = eltype(w)
    total = zero(T_val)
    for (ct, count) in type_counts
        val = w[type_to_idx[ct]]
        total += count * val
    end
    return total
end

"""
    get_matching_pair_partitions_filtered(indices)

Returns all pair partitions of 1..2k such that for every pair (a,b), indices[a] == indices[b].
"""
function get_matching_pair_partitions_filtered(indices::Vector{Int})
    n = length(indices)


    counts = Dict{Int,Vector{Int}}()
    for (pos, val) in enumerate(indices)
        push!(get!(counts, val, Int[]), pos)
    end

    for (val, pos_list) in counts
        if length(pos_list) % 2 != 0
            return Vector{Vector{Tuple{Int,Int}}}()
        end
    end



    group_partitions = Vector{Vector{Vector{Tuple{Int,Int}}}}()

    for (val, pos_list) in counts
        k_local = length(pos_list)
        # Generate partitions of 1..k_local
        local_parts_idx = get_pair_partitions(k_local)

        # Map to real positions
        real_parts = Vector{Vector{Tuple{Int,Int}}}()
        for p in local_parts_idx
            mapped = [(pos_list[a], pos_list[b]) for (a, b) in p]
            push!(real_parts, mapped)
        end
        push!(group_partitions, real_parts)
    end


    res = Vector{Vector{Tuple{Int,Int}}}()
    for combined in Iterators.product(group_partitions...)
        full_part = Vector{Tuple{Int,Int}}()
        for part in combined
            append!(full_part, part)
        end
        push!(res, full_part)
    end

    return res
end


"""
    integrate_indices_symplectic(indices, dim)

Low-level integration function using Symplectic Weingarten calculus.
Formula: sum_{pi, sigma} J_pi(i) * J_sigma(j) * Wg^Sp(pi, sigma)
"""
function integrate_indices_symplectic(indices::AbstractVector, dim)
    n = length(indices)
    if n % 2 != 0
        return 0
    end
    k = n ÷ 2
    partitions = get_pair_partitions(n)

    I = [x[1] for x in indices]
    J = [x[2] for x in indices]

    # Pre-calculate and group contractions
    pi_contractions = Dict{Any,Any}()
    for pi in partitions
        val_I = compute_symplectic_contraction(pi, I, dim)
        if !_iszero(val_I)
            c_pi = canonicalize_pair_partition(pi)
            pi_contractions[c_pi] = get(pi_contractions, c_pi, 0) + val_I
        end
    end

    if isempty(pi_contractions)
        return 0
    end

    sigma_contractions = Dict{Any,Any}()
    for sigma in partitions
        val_J = compute_symplectic_contraction(sigma, J, dim)
        if !_iszero(val_J)
            c_sigma = canonicalize_pair_partition(sigma)
            sigma_contractions[c_sigma] = get(sigma_contractions, c_sigma, 0) + val_J
        end
    end

    if isempty(sigma_contractions)
        return 0
    end

    w, type_to_idx, _ = get_weingarten_reduced_data(k, -dim)
    T = eltype(w)
    total = zero(T)

    n_pi = length(pi_contractions)
    n_sigma = length(sigma_contractions)
    if n_pi * n_sigma > 100
        p = Progress(n_pi * n_sigma; dt=10.0, desc="Calculating symplectic integrals... ")
        for (c_pi, val_pi) in pi_contractions
            for (c_sigma, val_sigma) in sigma_contractions
                ct = get_full_cycle_type(c_pi, c_sigma)
                loops = length(ct)
                wg = ((-1)^loops) * w[type_to_idx[ct]]
                total += val_pi * val_sigma * wg
                next!(p)
            end
        end
    else
        for (c_pi, val_pi) in pi_contractions
            for (c_sigma, val_sigma) in sigma_contractions
                ct = get_full_cycle_type(c_pi, c_sigma)
                loops = length(ct)
                wg = ((-1)^loops) * w[type_to_idx[ct]]
                total += val_pi * val_sigma * wg
            end
        end
    end

    if !(dim isa Integer)
        return Symbolics.simplify(Symbolics.wrap(total))
    end
    return total
end

"""
    integrate_indices_gue(indices, dim)
    
Low-level integration function using Wick's theorem for GUE.
Formula: sum_{pi in PairPartitions} prod_{(u, v) in pi} delta(i_u, j_v) * delta(j_u, i_v)
"""
function integrate_indices_gue(indices::AbstractVector, dim)
    n = length(indices) # Must be even
    partitions = get_pair_partitions(n)

    total = 0 // 1

    for pi in partitions
        term_val = 1
        possible = true

        for (u, v) in pi
            (i_u, j_u) = indices[u]
            (i_v, j_v) = indices[v]

            # Contraction < H_{i_u j_u} H_{i_v j_v} > = delta_{i_u j_v} * delta_{j_u i_v}
            # Check equalities
            if !_symbolic_isequal(i_u, j_v) || !_symbolic_isequal(j_u, i_v)
                possible = false
                break
            end
        end

        if possible
            total += 1
        end
    end

    return total
end

function integrate_indices_gue(indices::Vector{Any}, dim)
    return integrate_indices_gue(Vector{Tuple{Any,Any}}(indices), dim)
end

function integrate_indices_gue(indices::Vector{Tuple{Any,Any}}, dim)
    n = length(indices)
    partitions = get_pair_partitions(n)
    total = 0 // 1
    for pi in partitions
        possible = true
        for (u, v) in pi
            (i_u, j_u) = indices[u]
            (i_v, j_v) = indices[v]
            if !_symbolic_isequal(i_u, j_v) || !_symbolic_isequal(j_u, i_v)
                possible = false;
                break;
            end
        end
        if possible
            total += 1
        end
    end
    return total
end

"""
    integrate_indices_goe(indices, dim)
    
Low-level integration function using Wick's theorem for GOE.
Formula: sum_{pi in PairPartitions} prod_{(u, v) in pi} (delta(i_u, k_v)*delta(j_u, l_v) + delta(i_u, l_v)*delta(j_u, k_v))
where pair is H_{i_u j_u} and H_{k_v l_v} (indices re-labeled for clarity).
"""
function integrate_indices_goe(indices::AbstractVector, dim)
    n = length(indices) # Must be even
    partitions = get_pair_partitions(n)

    total = 0 // 1

    for pi in partitions
        term_val = 1
        possible = true

        for (u, v) in pi
            (i1, j1) = indices[u]
            (i2, j2) = indices[v]

            # GOE contraction: δ(i1,i2)δ(j1,j2) + δ(i1,j2)δ(j1,i2)
            val_pair = 0 // 1

            match1 = _symbolic_isequal(i1, i2) && _symbolic_isequal(j1, j2)
            match2 = _symbolic_isequal(i1, j2) && _symbolic_isequal(j1, i2)

            if match1
                val_pair += 1
            end
            if match2
                val_pair += 1
            end

            if val_pair == 0
                possible = false
                break
            end
            term_val *= val_pair
        end

        if possible
            total += term_val
        end
    end

    return total
end

# Removed helper methods that are no longer needed with relaxed signatures

function _get_J(i, j, d)
    # J = [0 I; -I 0], d even, n = d/2
    # J_{ij} = δ(i, j-n) - δ(i-n, j)
    if !(d isa Integer)
        return 0
    end
    n = d ÷ 2
    if i <= n && j > n && j == i + n
        return 1
    elseif i > n && j <= n && i == j + n
        return -1
    else
        return 0
    end
end

function integrate_indices_gse(indices::AbstractVector, dim)
    n = length(indices)
    partitions = get_pair_partitions(n)
    total = 0 // 1

    for pi in partitions
        # sum_{choices} (-1)^n2 * weight
        choice_combinations = collect(Iterators.product(fill([1, 2], n ÷ 2)...))
        for choices in choice_combinations
            term_val = 1 // 1
            possible = true
            n_type2 = 0
            for (p_idx, (u, v)) in enumerate(pi)
                (a, b) = indices[u]
                (c, d) = indices[v]

                choice = choices[p_idx]
                if choice == 1
                    # delta_ad delta_bc
                    if _symbolic_isequal(a, d) && _symbolic_isequal(b, c)
                        term_val *= 1
                    else
                        possible = false;
                        break
                    end
                else
                    # - J_ac J_bd
                    n_type2 += 1
                    jac = _get_J(a, c, dim)
                    jbd = _get_J(b, d, dim)
                    if jac == 0 || jbd == 0
                        possible = false;
                        break
                    end
                    term_val *= (jac * jbd)
                end
            end

            if possible
                total += term_val
            end
        end
    end
    return total
end

"""
    compute_symplectic_contraction(partition, indices, dim)

Computes prod_{(u, v) in partition} J(indices[u], indices[v]).
J = [0 I; -I 0].
"""
function compute_symplectic_contraction(partition, indices, dim)
    val = Num(1)

    for (u, v) in partition
        idx_u = indices[u]
        idx_v = indices[v]

        j_val = symplectic_form(idx_u, idx_v, dim)
        if _symbolic_isequal(j_val, 0)
            return 0
        end
        val *= j_val
    end
    return val
end

function symplectic_form(i, j, dim)
    # J matrix: J_{ij} = δ(j, i+m) - δ(j, i-m), m = dim/2

    u_dim = Symbolics.unwrap(dim)
    m = u_dim / 2

    if _symbolic_isequal(j, i + m)
        return 1
    elseif _symbolic_isequal(j, i - m)
        return -1
    else
        return 0
    end
end

"""
    integrate_indices_ginue(u_indices, ub_indices, dim)

Low-level integration for Complex Ginibre Ensemble.
Formula: sum_{sigma in S_n} prod_m delta(i_m, ib_sigma(m)) * delta(j_m, jb_sigma(m))
"""
function integrate_indices_ginue(
    u_indices::AbstractVector,
    ub_indices::AbstractVector,
    dim,
)
    n = length(u_indices)
    if n != length(ub_indices)
        return 0
    end
    perms = permutations(1:n)
    total = 0
    for p in perms
        possible = true
        for m = 1:n
            (i, j) = u_indices[m]
            (ib, jb) = ub_indices[p[m]]
            if !_symbolic_isequal(i, ib) || !_symbolic_isequal(j, jb)
                possible = false
                break
            end
        end
        if possible
            total += 1
        end
    end
    return total
end

"""
    integrate_indices_ginoe(indices, dim)

Low-level integration for Real Ginibre Ensemble.
Formula: sum_{pi in PairPartitions} prod_{(u, v) in pi} delta(i_u, i_v) * delta(j_u, j_v)
"""
function integrate_indices_ginoe(indices::AbstractVector, dim)
    n = length(indices)
    if isodd(n)
        return 0
    end

    partitions = get_pair_partitions(n)
    total = 0
    for pi in partitions
        possible = true
        for (u, v) in pi
            (i1, j1) = indices[u]
            (i2, j2) = indices[v]
            if !_symbolic_isequal(i1, i2) || !_symbolic_isequal(j1, j2)
                possible = false
                break
            end
        end
        if possible
            total += 1
        end
    end
    return total
end

"""
    integrate_indices_ginse(indices, dim)

Low-level integration for Symplectic Ginibre Ensemble.
Uses duality with GinOE.
"""
function integrate_indices_ginse(indices::AbstractVector, dim)
    n = length(indices)
    if isodd(n)
        return 0
    end

    res_oe = integrate_indices_ginoe(indices, dim)

    final_sign = ((-1)^(n ÷ 2 + 1))
    return final_sign * res_oe
end


"""
    _poly_degree(p, d)

Helper to get the degree of a polynomial `p` in variable `d`.
"""
function _poly_degree(p, d)
    deg = Symbolics.degree(Symbolics.wrap(p), Symbolics.wrap(d))
    return Int(Symbolics.unwrap(deg))
end


"""
    asymptotic(ex, d, order=1)

Generic asymptotic expansion of a rational function `ex` in powers of `1/d`.
"""
function asymptotic(ex, d, order = 1)
    return _expand_asymptotic(ex, d, order)
end

"""
    _expand_asymptotic(ex, d, order)

Helper to expand a rational function of `d` in powers of `1/d`.
Uses a series dictionary representation for robust recursive processing.
"""
function _expand_asymptotic(ex, d, order)
    ex_un = Symbolics.unwrap(ex)
    if ex_un isa AbstractArray
        return map(e -> _expand_asymptotic(e, d, order), ex_un)
    end
    if ex_un isa Complex
        re_part = _expand_asymptotic(real(ex_un), d, order)
        im_part = _expand_asymptotic(imag(ex_un), d, order)
        return re_part + im * im_part
    end
    if !(ex_un isa Symbolics.Num) &&
       !(ex_un isa SymbolicUtils.BasicSymbolic) &&
       !(ex_un isa Number)
        return ex_un
    end

    d_un = Symbolics.unwrap(d)
    d_num = Symbolics.wrap(d_un)

    series_dict = _asymptotic_series_dict(ex_un, d_un, order)

    if isempty(series_dict)
        return Num(0)
    end

    # Build Num from terms
    res = Num(0)
    for p in sort(collect(keys(series_dict)))
        coeff = series_dict[p]
        if !_iszero(coeff)
            if p == 0
                res = res + coeff
            elseif p > 0
                res = res + coeff * (1 / d_num)^p
            else
                res = res + coeff * d_num^(-p)
            end
        end
    end

    return Symbolics.expand(res)
end

function _asymptotic_series_dict(ex_un, d_un, order)
    if !(ex_un isa SymbolicUtils.BasicSymbolic)
        return Dict(0 => Symbolics.wrap(ex_un))
    end

    # Constants
    vars = Symbolics.get_variables(ex_un)
    if !any(v -> isequal(v, d_un), vars)
        return Dict(0 => Symbolics.wrap(ex_un))
    end

    if isequal(ex_un, d_un)
        return Dict(-1 => Num(1))
    end

    if SymbolicUtils.isadd(ex_un)
        res = Dict{Int,Any}()
        for arg in SymbolicUtils.arguments(ex_un)
            s_arg = _asymptotic_series_dict(arg, d_un, order)
            for (p, v) in s_arg
                if p <= order
                    res[p] = get(res, p, Num(0)) + v
                end
            end
        end
        return res
    end

    if SymbolicUtils.ismul(ex_un)
        res = Dict(0 => Num(1))
        for arg in SymbolicUtils.arguments(ex_un)
            s_arg = _asymptotic_series_dict(arg, d_un, order)
            new_res = Dict{Int,Any}()
            for (p1, v1) in res
                for (p2, v2) in s_arg
                    p = p1 + p2
                    if p <= order
                        val = v1 * v2
                        new_res[p] = get(new_res, p, Num(0)) + val
                    end
                end
            end
            res = new_res
            if isempty(res) break end
        end
        return res
    end

    if SymbolicUtils.ispow(ex_un)
        base, exp = SymbolicUtils.arguments(ex_un)
        if exp isa Integer && exp > 1
            # Simple integer power expansion
            res = Dict(0 => Num(1))
            s_base = _asymptotic_series_dict(Symbolics.unwrap(base), d_un, order)
            for _ = 1:exp
                new_res = Dict{Int,Any}()
                for (p1, v1) in res
                    for (p2, v2) in s_base
                        p = p1 + p2
                        if p <= order
                            new_res[p] = get(new_res, p, Num(0)) + v1 * v2
                        end
                    end
                end
                res = new_res
            end
            return res
        end
    end

    return _asymptotic_leaf_dict(ex_un, d_un, order)
end

function _asymptotic_leaf_dict(ex_un, d_un, order)
    d_num = Symbolics.wrap(d_un)
    # expand() is generally safer than simplify() for large mixed rational functions
    ex_sim = Symbolics.expand(Symbolics.wrap(ex_un))

    if ex_sim isa Complex
        re_dict = _asymptotic_series_dict(Symbolics.unwrap(real(ex_sim)), d_un, order)
        im_dict = _asymptotic_series_dict(Symbolics.unwrap(imag(ex_sim)), d_un, order)
        res = copy(re_dict)
        for (p, v) in im_dict
            res[p] = get(res, p, Num(0)) + im * v
        end
        return res
    end

    num = Symbolics.numerator(ex_sim)
    den = Symbolics.denominator(ex_sim)

    n = _poly_degree(num, d_un)
    m = _poly_degree(den, d_un)

    ϵ = Symbolics.variable(:ϵ)
    ϵ_un = Symbolics.unwrap(ϵ)

    p_eps = Symbolics.expand(Symbolics.substitute(num, Dict(d_num => 1 / ϵ)) * ϵ^n)
    q_eps = Symbolics.expand(Symbolics.substitute(den, Dict(d_num => 1 / ϵ)) * ϵ^m)

    f_analytic = Symbolics.simplify(p_eps / q_eps; expand = true)

    res_dict = Dict{Int,Any}()
    curr_deriv = f_analytic
    diff = m - n

    max_k = order - diff

    for k = 0:max_k
        v = Symbolics.substitute(curr_deriv, Dict(ϵ => 0))
        vu = Symbolics.unwrap(v)

        needs_sim = false
        if vu isa Number && (isnan(vu) || isinf(vu))
            needs_sim = true
        elseif !(vu isa Number)
            s_val = string(vu)
            if occursin("NaN", s_val) || occursin("Inf", s_val) || occursin("1//0", s_val)
                needs_sim = true
            end
        end

        if needs_sim
            curr_sim = Symbolics.expand(curr_deriv)
            val = Symbolics.substitute(curr_sim, Dict(ϵ => 0))
            if any(var -> isequal(var, ϵ_un), Symbolics.get_variables(val))
                val = Symbolics.wrap(
                    SymbolicUtils.Postwalk(x -> _symbolic_isequal(x, ϵ_un) ? 0 : x)(
                        Symbolics.unwrap(curr_deriv),
                    ),
                )
            end
        else
            val = v
        end

        if !_iszero(val)
            power = k + diff
            coeff = (1 // factorial(k))
            res_dict[power] = get(res_dict, power, Num(0)) + Symbolics.expand(val * coeff)
        end

        if k < max_k
            curr_deriv = Symbolics.derivative(curr_deriv, ϵ)
            if k % 2 == 0
                curr_deriv = Symbolics.expand(curr_deriv)
            end
        end
    end
    return res_dict
end

"""
    integrate_indices_diagonal(U_idxs, U_bar_idxs, dim)

Low-level integration function for the Diagonal Unitary group (Torus).
The integral is non-zero (equal to 1) only if the multisets of indices 
of U and U_bar are identical.
Indices are already checked to be diagonal (i == j) in process_term.
"""
function integrate_indices_diagonal(
    U_idxs::AbstractVector,
    U_bar_idxs::AbstractVector,
    dim,
)
    n = length(U_idxs)
    if n != length(U_bar_idxs)
        return 0
    end

    # Only rows matter since indices are diagonal
    rows_u = [x[1] for x in U_idxs]
    rows_ubar = [x[1] for x in U_bar_idxs]
    sort!(rows_u)
    sort!(rows_ubar)

    if rows_u == rows_ubar
        return 1
    else
        return 0
    end
end

function _standardize_sub(k)
    uk = Symbolics.unwrap(k)
    if uk isa LazyTrace
        if length(uk.cycles) == 1
            # This is a bit hacky as it depends on tr_val being in scope 
            # and implemented the way it is.
            return tr_val(uk.cycles[1])
        end
    end
    return k
end

function _get_denominators(ex)
    ex_un = Symbolics.unwrap(ex)
    dens = Num[]
    if SymbolicUtils.iscall(ex_un)
        op = SymbolicUtils.operation(ex_un)
        args = SymbolicUtils.arguments(ex_un)
        if op == (/)
            push!(dens, Symbolics.wrap(args[2]))
        end
        for a in args
            append!(dens, _get_denominators(a))
        end
    end
    return dens
end

"""
    evaluate(expr, dict)
    evaluate(expr, pair)

Shorthand for `Symbolics.substitute`. Useful for substituting symbolic dimensions 
with numeric values in integration results. Also handles substituting symbolic traces.

This function automatically handles removable singularities in fractions (e.g., \$0/0\$ forms).
If a denominator evaluates to zero after substitution, the expression is first simplified 
to attempt resolving the singularity before completing the evaluation.
"""
function evaluate(expr, dict)
    new_dict = if dict isa AbstractDict || dict isa AbstractVector
        Dict(_standardize_sub(k) => v for (k, v) in dict)
    else
        dict
    end

    # Handle removable singularities in fractions
    # Check if any denominator is zero after substitution
    dens = _get_denominators(expr)
    should_simplify = false
    for d in dens
        # d_val can be a Num or a Number
        d_val = Symbolics.substitute(d, new_dict)
        if _iszero(d_val)
            should_simplify = true
            break
        end
    end

    if should_simplify
        # Simplify first to resolve removable singularities like (x^2-1)/(x-1) -> x+1
        expr = Symbolics.simplify(expr)
    end

    res = Symbolics.substitute(expr, new_dict)

    # Try to return a number if the result is a Num wrapping a number
    if res isa Num
        val = Symbolics.unwrap(res)
        if val isa Number
            return val
        end

        if hasproperty(val, :val) ? val.val isa Number : false
            return val.val
        end
        if SymbolicUtils.iscall(val)
            op = SymbolicUtils.operation(val)
            args = SymbolicUtils.arguments(val)
            # Safe evaluation for basic arithmetic if all args are numeric
            if (op === (+) || op === (*) || op === (-) || op === (/)) &&
               all(x -> x isa Number, args)
                return op(args...)
            end
        end

        # Also try simplifying to see if it becomes a number
        # Note: simplify can return a Num, so we check again
        sim_res = Symbolics.simplify(res)
        val_sim = Symbolics.unwrap(sim_res)
        if val_sim isa Number
            return val_sim
        end
    end
    return res
end

function evaluate(expr, pair::Pair)
    return evaluate(expr, Dict(pair))
end
