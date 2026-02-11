_to_Num(z::Complex) = Complex(Num(real(z)), Num(imag(z)))
_to_Num(z) = Num(z)

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


    try
        if a_un == b_un
            return true
        end
    catch
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
    try
        s = Symbolics.simplify(x)
        sv = Symbolics.value(s)
        if sv isa Number
            return iszero(sv)
        end
    catch
    end
    return _symbolic_isequal(v, 0)
end

function _robust_real(x)
    if x isa AbstractArray
        return map(_robust_real, x)
    end

    x_un = Symbolics.unwrap(x)


    v = Symbolics.value(x_un)
    if v isa AbstractFloat
        return rationalize(v, tol = 1e-13)
    end
    if v isa Real
        return v
    end
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


    try
        nx = _safe_Num(x_un)


        v = Symbolics.value(nx)
        if v isa AbstractFloat
            return rationalize(v, tol = 1e-13)
        end
        if v isa Real
            return v
        end


        nx = Symbolics.simplify(nx)
        v = Symbolics.value(nx)
        if v isa AbstractFloat
            return rationalize(v, tol = 1e-13)
        end
        if v isa Real
            return v
        end


        if _iszero(Symbolics.simplify(imag(nx)))
            rx = Symbolics.simplify(real(nx))
            vx = Symbolics.value(rx)
            if vx isa AbstractFloat
                return rationalize(vx, tol = 1e-13)
            end
            if vx isa Real
                return vx
            end
            return rx
        end
        return nx
    catch
    end

    return x_un
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


    try
        v = Symbolics.value(x)
        if v !== x
            return _is_manifestly_real(v)
        end
    catch
    end


    return false
end

"""
    AbstractIndexMatcher

Base type for index matching strategies. Subtypes must implement `match_index`.
"""
abstract type AbstractIndexMatcher end

"""
    LookupMatcher(U_lookup, U_bar_lookup)

A matcher that uses dictionaries to map symbolic variables to (row, column) indices.
Used primarily for GUE/GOE/GSE and list-based unitary integration.
"""
struct LookupMatcher <: AbstractIndexMatcher
    U_lookup::Dict{Any,Tuple}
    U_bar_lookup::Dict{Any,Tuple}
end

function match_index(m::LookupMatcher, t)
    t_un = Symbolics.unwrap(t)


    if haskey(m.U_lookup, t_un)
        v = m.U_lookup[t_un]

        if length(v) == 3 && v[3] == :conj
            return (:U_bar, v[1], v[2])
        else
            return (:U, v[1], v[2])
        end
    end


    if haskey(m.U_bar_lookup, t_un)
        v = m.U_bar_lookup[t_un]
        return (:U_bar, v[1], v[2])
    end

    return nothing
end

"""
    INTEGRATION_RULES

A dictionary mapping measure types (symbols) to their respective integration rule functions.
Each rule function should have the signature `(u_indices, u_bar_indices, dim, measure_type)`.
"""
const INTEGRATION_RULES = Dict{Any,Function}()

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
        return _robust_real(val_re + im * val_im)
    end

    expr_un = Symbolics.unwrap(expr)
    if expr_un isa Number && (expr_un isa Real || expr_un isa Complex)
        return expr_un
    end


    if expr_un isa Complex
        return _integrate_core(expr_un, dim, subs_dict, matcher, measure_type)
    end


    expr_num = _safe_Num(expr_un)
    try
        expr_num = Symbolics.expand(expr_num)
    catch
    end


    function is_add(t)
        Symbolics.iscall(t) && Symbolics.operation(t) == (+)
    end

    # Rewrite rules: expand abs², real, imag into polynomial form
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

    expr_unwrapped = Symbolics.unwrap(expr)


    chain = SymbolicUtils.Chain([
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
        r_pow_nested,
        r_hypot_pow,
        r_float_to_int_pow,
    ])
    expr_rewritten = SymbolicUtils.Postwalk(
        SymbolicUtils.PassThrough(chain);
        maketerm = (st, f, args, metadata; kwargs...) -> begin
            if f == complex || f == Base.complex
                return args[1] + im * args[2]
            end
            SymbolicUtils.maketerm(st, f, args, metadata)
        end,
    )(
        expr_unwrapped,
    )
    if expr_rewritten isa Complex
        return _integrate_core(expr_rewritten, dim, subs_dict, matcher, measure_type)
    end
    expr_num = _safe_Num(expr_rewritten)



    try
        expr_num = Symbolics.expand(_safe_Num(Symbolics.unwrap(expr_num)))
    catch
    end


    function robust_substitute(ex, dict)
        try
            res = Symbolics.substitute(Symbolics.wrap(ex), dict)

            if _symbolic_isequal(res, Symbolics.wrap(ex))
                throw(error("No change"))
            end
            return res
        catch
            # Fallback: manual Postwalk traversal with unwrapped keys
            unwrapped_dict =
                Dict(Symbolics.unwrap(k) => Symbolics.unwrap(v) for (k, v) in dict)

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
    end

    expr_subbed = if expr_num isa Complex
        robust_substitute(Symbolics.unwrap(real(expr_num)), subs_dict) +
        im * robust_substitute(Symbolics.unwrap(imag(expr_num)), subs_dict)
    else
        robust_substitute(Symbolics.unwrap(expr_num), subs_dict)
    end



    expanded_expr = try
        Symbolics.expand(_safe_Num(Symbolics.unwrap(expr_subbed)))
    catch
        _safe_Num(expr_subbed)
    end

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


    expanded_expr = try
        Symbolics.expand(expanded_expr)
    catch
        expanded_expr
    end


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

function integrate(expr::LazySum, measure)
    return sum(t -> integrate(t, measure), expr.terms)
end

"""
    integrate(expr::AbstractArray, measure)

Performs element-wise integration of a matrix or array of expressions.
Returns an array of the same shape containing integrated values.
"""
function integrate(expr::AbstractArray, measure)
    return map(t -> integrate(t, measure), expr)
end

"""
    integrate(expr, measure)

Top-level integration function. It first checks the [Pre-computed Integral Library](@ref) 
for instant results. If not found, it calls `fallback_integrate` for the specific measure.
"""
function integrate(expr, measure)
    lib_res = check_library(expr, measure)
    if lib_res !== nothing
        return lib_res
    end

    return fallback_integrate(expr, measure)
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
    if !(x isa Number) && !(x isa AbstractArray)
        return Num(x)
    end
    x_un = Symbolics.unwrap(x)
    return _to_Num(x_un)
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
    u_indices = Vector{Tuple{Int,Int}}()
    u_bar_indices = Vector{Tuple{Int,Int}}()

    function traverse(t)
        t_unwrapped = Symbolics.unwrap(t)

        match_res = match_index(matcher, t_unwrapped)
        if match_res !== nothing
            type, i, j = match_res
            if type == :U
                push!(u_indices, (i, j))
            else
                push!(u_bar_indices, (i, j))
            end
            return
        end

        if t isa Number && !(t isa Num) && !(t isa Complex{Num})
            coeff *= t
            return
        end

        if Symbolics.iscall(t_unwrapped)
            local op = Symbolics.operation(t_unwrapped)
            local args = Symbolics.arguments(t_unwrapped)

            if op == (*)
                for arg in args
                    traverse(arg)
                end
                return
            elseif op == (^)
                base = args[1]
                p_val = Symbolics.unwrap(args[2])
                p = try
                    parse(Int, string(p_val))
                catch
                    nothing
                end
                if p isa Integer
                    for _ = 1:p
                        traverse(base)
                    end
                    return
                end
                traverse(base)
                return
            elseif op == (/)
                traverse(args[1])
                coeff /= args[2]
                return
            elseif op == conj || op == Base.conj
                inner = Symbolics.unwrap(args[1])
                match_res_inner = match_index(matcher, inner)
                if match_res_inner !== nothing
                    type, i, j = match_res_inner
                    if type == :U
                        push!(u_bar_indices, (i, j))
                    else
                        push!(u_indices, (i, j))
                    end
                    return
                end
            elseif op == complex || op == Base.complex || op == (+) || op == (-)
                for arg in args
                    traverse(arg)
                end
                return
            end

        end

        coeff *= t
    end
    traverse(term)

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
        length(u) != length(ub) ? 0 : (length(u) == 0 ? 1 : integrate_indices(u, ub, d))
    end

INTEGRATION_RULES[:O] =
    (u, ub, d, mt) -> begin
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_orthogonal(all_indices, d))
    end

INTEGRATION_RULES[:Sp] =
    (u, ub, d, mt) -> begin
        # Map u_bar indices to u using symplectic metric relations
        # S_bar_{ij} = coeff * S_{p(i)p(j)}
        # This requires d to be an integer (or at least p(i) to be resolvable)
        
        # We need _j_pair_sign logic here. It is not exported from real_measures?
        # It is likely defined in IntU. Let's assume we can call it or reimplement it properly.
        # Check if we can find _j_pair_sign in integration_core or if we need to define it.
        # It's not in integration_core. It was in real_measures.
        # So we cannot call it inside integration_core if it's defined later/elsewhere without import.
        # But integration_core is included BEFORE real_measures usually?
        # Actually in partial viewing of real_measures.jl, it was defined there.
        # This means integration_core CANNOT depend on it if integration_core is loaded first.
        # I should define a helper here or check if it's available.
        # To be safe, I will implement the mapping logic locally here using a helper 
        # or just inline it since it's simple.
        
        mapped_ub = Vector{Tuple{Int,Int}}()
        total_coeff = 1
        
        if !isempty(ub) 
            if !(d isa Integer)
               # Cannot resolve symplectic pairs for symbolic dimension with explicit integer indices
               # Fallback to 0 or error? 
               # If indices were symbolic, maybe? But we have Int indices.
               return 0 
            end
            
            n_half = div(d, 2)
            
            for (i, j) in ub
                 # _j_pair_sign logic
                 p = (i <= n_half) ? i + n_half : i - n_half
                 sign_i = (i <= n_half) ? 1 : -1
                 
                 q = (j <= n_half) ? j + n_half : j - n_half
                 sign_j = (j <= n_half) ? 1 : -1
                 
                 push!(mapped_ub, (p, q))
                 total_coeff *= (sign_i * sign_j)
            end
        end

        all_indices = [u; mapped_ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : total_coeff * integrate_indices_symplectic(all_indices, d))
    end

INTEGRATION_RULES[:GUE] =
    (u, ub, d, mt) -> begin
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_gue(all_indices, d))
    end

INTEGRATION_RULES[:GOE] =
    (u, ub, d, mt) -> begin
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_goe(all_indices, d))
    end

INTEGRATION_RULES[:GSE] =
    (u, ub, d, mt) -> begin
        all_indices = [u; ub]
        length(all_indices) % 2 != 0 ? 0 :
        (length(all_indices) == 0 ? 1 : integrate_indices_gse(all_indices, d))
    end

INTEGRATION_RULES[:COE] =
    (u, ub, d, mt) -> begin
        length(u) != length(ub) ? 0 : (length(u) == 0 ? 1 : integrate_indices_coe(u, ub, d))
    end

INTEGRATION_RULES[:CSE] =
    (u, ub, d, mt) -> begin
        length(u) != length(ub) ? 0 :
        (
            begin
                phys_dim = mt isa Tuple ? mt[2] : d
                length(u) == 0 ? 1 : integrate_indices_cse(u, ub, d, phys_dim)
            end
        )
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

INTEGRATION_RULES[:Perm] =
    (u, ub, d, mt) -> begin
        all_indices = [u; ub]
        length(all_indices) == 0 ? 1 : integrate_indices_permutation(all_indices, d)
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
function integrate_indices(
    U_idxs::Vector{Tuple{Int,Int}},
    U_bar_idxs::Vector{Tuple{Int,Int}},
    dim,
)
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
    for sigma in valid
        for tau in valid_taus
            inv_tau = invperm(tau)
            P = [sigma[inv_tau[i]] for i = 1:n]
            cycle_type = get_cycle_type(P)
            cycle_counts[cycle_type] = get(cycle_counts, cycle_type, 0) + 1
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

function get_matching_permutations(target::Vector{Int}, source::Vector{Int})
    n = length(target)
    if n != length(source)
        return Vector{Vector{Int}}()
    end


    source_groups = Dict{Int,Vector{Int}}()
    for (idx, val) in enumerate(source)
        push!(get!(source_groups, val, Int[]), idx)
    end

    target_groups = Dict{Int,Vector{Int}}()
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
    visited = falses(length(p))
    lengths = Int[]
    for i = 1:length(p)
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

"""
    integrate_indices_orthogonal(indices, dim)
    
Low-level integration function using Orthogonal Weingarten calculus.
Indices are a list of (i, j) for O_{ij}.
Formula: sum_{pi, sigma in PairPartitions} delta_pi(i) * delta_sigma(j) * Wg(pi, sigma)
"""
function integrate_indices_orthogonal(indices::Vector{Tuple{Int,Int}}, dim)
    n = length(indices) # This is 2k
    if n % 2 != 0
        return 0
    end

    I = [x[1] for x in indices]
    J = [x[2] for x in indices]

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

    val_mat, lookup = get_weingarten_orthogonal_data(n ÷ 2, dim)

    total = 0 // 1
    for (c_pi, count_pi) in pi_counts
        idx_pi = get(lookup, c_pi, nothing)
        idx_pi === nothing && continue

        for (c_sigma, count_sigma) in sigma_counts
            idx_sigma = get(lookup, c_sigma, nothing)
            idx_sigma === nothing && continue

            val = val_mat[idx_pi, idx_sigma]
            total += (count_pi * count_sigma) * val
        end
    end
    if !(dim isa Integer)
        try
            return Symbolics.simplify(Symbolics.wrap(total))
        catch
        end
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
function integrate_indices_symplectic(indices::Vector{Tuple{Int,Int}}, dim)
    n = length(indices)
    k = div(n, 2)
    partitions = get_pair_partitions(n)

    I = [x[1] for x in indices]
    J_idx = [x[2] for x in indices]

    dim_int = dim isa Integer ? dim : nothing

    total = 0 // 1

    pi_contractions = Dict{Any,Any}()
    sigma_contractions = Dict{Any,Any}()

    for pi in partitions
        val_I = compute_symplectic_contraction(pi, I, dim)
        if !_symbolic_isequal(val_I, 0)
            pi_contractions[pi] = val_I
        end

        val_J = compute_symplectic_contraction(pi, J_idx, dim)
        if !_symbolic_isequal(val_J, 0)
            sigma_contractions[pi] = val_J
        end
    end

    if isempty(pi_contractions) || isempty(sigma_contractions)
        return 0
    end


    total = 0 // 1

    for (pi, val_pi) in pi_contractions
        for (sigma, val_sigma) in sigma_contractions
            wg = weingarten_symplectic_val(pi, sigma, dim)
            total += val_pi * val_sigma * wg
        end
    end

    if !(dim isa Integer)
        try
            return Symbolics.simplify(Symbolics.wrap(total))
        catch
        end
    end
    return total
end

"""
    integrate_indices_gue(indices, dim)
    
Low-level integration function using Wick's theorem for GUE.
Formula: sum_{pi in PairPartitions} prod_{(u, v) in pi} delta(i_u, j_v) * delta(j_u, i_v)
"""
function integrate_indices_gue(indices::Vector{Tuple{Int,Int}}, dim)
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
function integrate_indices_goe(indices::Vector{Tuple{Int,Int}}, dim)
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

function integrate_indices_goe(indices::Vector{Any}, dim)
    return integrate_indices_goe(Vector{Tuple{Any,Any}}(indices), dim)
end

function integrate_indices_goe(indices::Vector{Tuple{Any,Any}}, dim)
    n = length(indices) # Must be even
    partitions = get_pair_partitions(n)

    total = 0 // 1

    for pi in partitions
        term_val = 1
        possible = true

        for (u, v) in pi
            (i1, j1) = indices[u]
            (i2, j2) = indices[v]

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

function integrate_indices_gse(indices::Vector{Any}, dim)
    return integrate_indices_gse(Vector{Tuple{Int,Int}}(indices), dim)
end

function integrate_indices_gse(indices::Vector{Tuple{Int,Int}}, dim)
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
    val = 1

    u_unwrapped = Symbolics.unwrap(dim)
    is_d_numeric = (u_unwrapped isa Number)

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

    if !(i isa Integer) || !(j isa Integer)
        # Symbolic indices not fully supported yet for contraction
        return 0
    end


    u_dim = Symbolics.unwrap(dim)
    if !(u_dim isa Number)
        return 0
    end

    dim_val = Int(u_dim)
    m = div(dim_val, 2)


    if i < 1 || i > 2*m || j < 1 || j > 2*m
        return 0
    end

    if j == i + m
        return 1
    elseif j == i - m
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
    u_indices::Vector{Tuple{Int,Int}},
    ub_indices::Vector{Tuple{Int,Int}},
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
function integrate_indices_ginoe(indices::Vector{Tuple{Int,Int}}, dim)
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
function integrate_indices_ginse(indices::Vector{Tuple{Int,Int}}, dim)
    n = length(indices)
    if isodd(n)
        return 0
    end

    res_oe = integrate_indices_ginoe(indices, dim)

    final_sign = ((-1)^(n ÷ 2 + 1))
    return final_sign * res_oe
end


# Overloads for symbolic arrays
"""
    tr(A)

Compute the trace of a matrix. Works for both standard matrices (via `LinearAlgebra.tr`) 
and symbolic arrays (by summing diagonal elements).
"""
function tr(A::Symbolics.Arr{T,2}) where {T}
    return sum(A[i, i] for i = 1:size(A, 1))
end

function tr(A)
    return LinearAlgebra.tr(A)
end


"""
    _poly_degree(p, d)

Helper to get the degree of a polynomial `p` in variable `d`.
"""
function _poly_degree(p, d)
    try
        deg = Symbolics.degree(Symbolics.wrap(p), Symbolics.wrap(d))
        return Int(Symbolics.unwrap(deg))
    catch
        # Fallback: basic manual check
        p_un = Symbolics.unwrap(p)
        d_un = Symbolics.unwrap(d)
        if _symbolic_isequal(p_un, d_un)
            return 1
        end
        return 0
    end
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
Uses Taylor expansion in ε = 1/d around ε = 0.
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

    # Substitute d -> 1/ε and clear denominators
    ex_sim = try
        Symbolics.simplify(Symbolics.wrap(ex_un))
    catch
        Symbolics.wrap(ex_un)
    end


    if ex_sim isa Complex
        re_part = _expand_asymptotic(real(ex_sim), d, order)
        im_part = _expand_asymptotic(imag(ex_sim), d, order)
        return re_part + im * im_part
    end

    num = Symbolics.numerator(ex_sim)
    den = Symbolics.denominator(ex_sim)


    n = _poly_degree(num, d_un)
    m = _poly_degree(den, d_un)

    ϵ = Symbolics.variable(:ϵ)
    ϵ_un = Symbolics.unwrap(ϵ)


    p_eps = try
        Symbolics.simplify(Symbolics.substitute(num, Dict(d_num => 1/ϵ)) * ϵ^n; expand = true)
    catch
        Symbolics.substitute(num, Dict(d_num => 1/ϵ)) * ϵ^n
    end
    q_eps = try
        Symbolics.simplify(Symbolics.substitute(den, Dict(d_num => 1/ϵ)) * ϵ^m; expand = true)
    catch
        Symbolics.substitute(den, Dict(d_num => 1/ϵ)) * ϵ^m
    end


    f_analytic = try
        Symbolics.simplify(p_eps / q_eps; expand = true)
    catch
        p_eps / q_eps
    end

    total = Num(0)
    curr_deriv = f_analytic
    diff = m - n

    max_k = order - diff

    for k = 0:max_k
        val = try
            v = Symbolics.substitute(curr_deriv, Dict(ϵ => 0))
            vu = Symbolics.unwrap(v)
            if isnan(vu) || isinf(vu)
                curr_sim = Symbolics.simplify(curr_deriv; expand = true)
                Symbolics.substitute(curr_sim, Dict(ϵ => 0))
            else
                v
            end
        catch
            try
                curr_sim = Symbolics.simplify(curr_deriv; expand = true)
                Symbolics.substitute(curr_sim, Dict(ϵ => 0))
            catch
                # Last resort: postwalk replacement
                Symbolics.wrap(
                    SymbolicUtils.Postwalk(x -> _symbolic_isequal(x, ϵ_un) ? 0 : x)(
                        Symbolics.unwrap(curr_deriv),
                    ),
                )
            end
        end

        if !_iszero(val)
            power = k + diff
            term = (val * (1//factorial(k))) * (1/d_num)^power
            total += term
        end

        if k < max_k
            curr_deriv = Symbolics.derivative(curr_deriv, ϵ)
            # Periodic simplification to manage expression size
            if k % 2 == 0
                try
                    curr_deriv = Symbolics.simplify(curr_deriv; expand = true)
                catch
                end
            end
        end
    end

    return Symbolics.simplify(total)
end

"""
    integrate_indices_diagonal(U_idxs, U_bar_idxs, dim)

Low-level integration function for the Diagonal Unitary group (Torus).
The integral is non-zero (equal to 1) only if the multisets of indices 
of U and U_bar are identical.
Indices are already checked to be diagonal (i == j) in process_term.
"""
function integrate_indices_diagonal(
    U_idxs::Vector{Tuple{Int,Int}},
    U_bar_idxs::Vector{Tuple{Int,Int}},
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
