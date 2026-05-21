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
    _concrete_numeric_value(x)

Return the concrete numeric payload of `x` when available, otherwise `nothing`.
Works for plain numbers and concrete symbolic constants (e.g. `Num(2.0)`).
"""
function _concrete_numeric_value(x)
    if x isa Num
        v = Symbolics.value(x)
        return v isa Number ? v : nothing
    end
    if x isa Number
        return x
    end
    v = try
        Symbolics.value(x)
    catch
        nothing
    end
    return v isa Number ? v : nothing
end

"""
    _assert_no_float_param(x, param_name, context="value")

Reject float-valued measure parameters. Accepts exact integer-like numeric values
and symbolic variables, but throws for concrete floats (including wrapped
constants like `Num(2.0)`).
"""
function _assert_no_float_param(
    x,
    param_name::AbstractString,
    context::AbstractString = "value",
)
    v = _concrete_numeric_value(x)
    if v isa AbstractFloat
        throw(ArgumentError(
            "$context requires `$param_name` to be an exact integer-like value " *
            "or symbolic variable; got float input $x."
        ))
    end
    return x
end

"""
    _try_extract_int(x)

Extract a plain `Int` from `x` if it wraps a concrete integer value.
Returns the `Int`, or `nothing` if `x` is symbolic or non-integer.
Handles: plain `Integer`, `Num`-wrapped integer constants (e.g. `Num(2)`).
"""
function _try_extract_int(x)
    if x isa Integer
        return typemin(Int) <= x <= typemax(Int) ? Int(x) : nothing
    end
    if x isa AbstractFloat && isinteger(x) && typemin(Int) <= x <= typemax(Int)
        return trunc(Int, x)
    end
    if x isa Rational && isinteger(x)
        n = numerator(x)
        return typemin(Int) <= n <= typemax(Int) ? Int(n) : nothing
    end
    if x isa Num
        u = Symbolics.unwrap(x)
        if !SymbolicUtils.issym(u) && !SymbolicUtils.iscall(u)
            st = SymbolicUtils.symtype(u)
            val = getfield(getfield(u, 1), 1)
            if st <: Integer
                return typemin(Int) <= val <= typemax(Int) ? Int(val) : nothing
            end
            if (st <: AbstractFloat || st <: Rational) && isinteger(val)
                return _try_extract_int(val)
            end
        end
    end
    # Fallback: try parsing the string representation (handles BasicSymbolic constants)
    str = string(x)
    return tryparse(Int, str)
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
        isfinite(v) || return v
        # Only rationalize if the result is exact (tol=0) and has a reasonable
        # denominator. This avoids silently zeroing tiny values like 1e-20.
        r = rationalize(v, tol = 0)
        if isfinite(r) && abs(denominator(r)) <= 10^15
            return r
        end
        return v
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

Abstract base type for strategies that identify random-matrix entries inside a
symbolic expression. Concrete subtypes must implement

    match_index(m::AbstractIndexMatcher, t) -> Union{Tuple, Nothing}

returning `(tag, i, j)` when `t` is a recognised random-matrix entry
(where `tag` is a `Symbol` like `:U` or `:U_bar`, and `i`, `j` are row/column
indices or `nothing` for metadata-only matching), or `nothing` otherwise.

The built-in subtype is [`MetadataMatcher`](@ref).
"""
abstract type AbstractIndexMatcher end
"""
    MetadataMatcher(type_tag::Symbol)

An [`AbstractIndexMatcher`](@ref) that recognises random-matrix entries by the
`special_type` metadata attached to a `SymbolicMatrix` symbol.

`type_tag` must match the `special_type` field of the target matrix, e.g.:
- `:U` for Haar-unitary matrices (`symbolic_unitary`)
- `:O` for orthogonal matrices (`symbolic_orthogonal`)
- `:Sp` for symplectic matrices (`symbolic_symplectic`)
- `:psi` for pure-state vectors (`symbolic_pure_state`)

Conjugate entries (`is_adj = true`) are tagged as `Symbol(type_tag, :_bar)`
(e.g. `:U_bar`).
"""
struct MetadataMatcher <: AbstractIndexMatcher
    type_tag::Symbol
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
