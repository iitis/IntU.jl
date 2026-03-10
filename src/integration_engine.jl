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
    end)(Symbolics.unwrap(ex))

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
                    p = Progress(n_terms; dt = 10.0, desc = "Integrating terms... ")
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
        throw(
            ArgumentError(
                "Graphical integration (LazyTrace) not implemented for measure type $measure_type. Try expanding traces element-wise.",
            ),
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
        throw(ArgumentError("Unknown measure type: $measure_type"))
    end
end
