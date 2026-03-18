
"""
    _poly_degree(p, d)

Helper to get the degree of a polynomial `p` in variable `d`.
"""
function _poly_degree(p, d)
    deg = Symbolics.degree(Symbolics.wrap(p), Symbolics.wrap(d))
    return Int(Symbolics.unwrap(deg))
end
"""
    asymptotic(expr, measure::AbstractMeasure, order=1)

Returns the series expansion of the integral in powers of `1/d`.
Works for any measure that implements `_reconstruct_symbolic`.
"""
function asymptotic(expr, measure::AbstractMeasure, order = 1)
    d = measure.dim
    if d isa SymbolicMatrix
        d = d.dim
    end

    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = _reconstruct_symbolic(measure, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
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
            if isempty(res)
                break
            end
        end
        return res
    end

    if SymbolicUtils.ispow(ex_un)
        base, exp = SymbolicUtils.arguments(ex_un)
        if exp isa Integer && exp > 1
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

function _standardize_sub(k)
    uk = Symbolics.unwrap(k)
    if uk isa LazyTrace
        if length(uk.cycles) == 1
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

    dens = _get_denominators(expr)
    should_simplify = false
    for d in dens
        d_val = Symbolics.substitute(d, new_dict)
        if _iszero(d_val)
            should_simplify = true
            break
        end
    end

    if should_simplify
        expr = Symbolics.simplify(expr)
    end

    res = Symbolics.substitute(expr, new_dict)

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
            if (op === (+) || op === (*) || op === (-) || op === (/)) &&
               all(x -> x isa Number, args)
                return op(args...)
            end
        end

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
