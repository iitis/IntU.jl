
"""
    check_library(expr, measure)

Check if the integral of `expr` over `measure` is available in the pre-computed library.
Returns the symbolic result if found, otherwise `nothing`.
"""
check_library(expr, measure) = nothing
check_library(expr, measure::HaarMeasure) = check_haar_library(expr, measure)
check_library(expr, measure::GUEMeasure) = check_gaussian_library(expr, measure, :GUE)
check_library(expr, measure::GOEMeasure) = check_gaussian_library(expr, measure, :GOE)
check_library(expr, measure::GSEMeasure) = check_gaussian_library(expr, measure, :GSE)
check_library(expr, measure::GinUEMeasure) = check_ginibre_library(expr, measure, :GinUE)
check_library(expr, measure::GinOEMeasure) = check_ginibre_library(expr, measure, :GinOE)
check_library(expr, measure::GinSEMeasure) = check_ginibre_library(expr, measure, :GinSE)
check_library(expr, measure::OrthogonalMeasure) = check_orthogonal_library(expr, measure)
check_library(expr, measure::SymplecticMeasure) = check_symplectic_library(expr, measure)
check_library(expr, measure::COEMeasure) = check_coe_library(expr, measure)
check_library(expr, measure::CSEMeasure) = check_cse_library(expr, measure)
check_library(expr, measure::PureStateMeasure) = check_pure_library(expr, measure)


function _haar_trace_moment_value(k::Int, d::Integer)
    total = zero(Rational{BigInt})
    for part in partitions(k)
        if length(part) <= d
            f_lambda = character_at_id(part)
            total += f_lambda^2
        end
    end
    return total
end

function _match_haar_pure_trace_moment(expr::LazyTrace, measure::HaarMeasure)
    u_names = Set{Symbol}()
    n_u = 0
    n_u_dag = 0

    for cycle in expr.cycles
        for f in cycle
            if !(f isa SymbolicMatrix) || f.special_type != :U
                return nothing
            end
            if !isequal(f.dim, measure.dim)
                return nothing
            end
            push!(u_names, f.name)
            if f.is_adj
                n_u_dag += 1
            else
                n_u += 1
            end
        end
    end

    if length(u_names) != 1 || n_u == 0 || n_u != n_u_dag
        return nothing
    end

    return n_u
end

function check_haar_library(expr, measure)
    if expr isa LazyTrace
        prefactor = expr.prefactor

        k = _match_haar_pure_trace_moment(expr, measure)
        if k !== nothing
            d = measure.dim
            if d isa Integer
                return prefactor * _haar_trace_moment_value(k, d)
            elseif k <= 1
                return prefactor * one(BigInt)
            else
                throw(
                    ArgumentError(
                        "|tr(U)|^$(2k) requires a concrete integer dimension " *
                        "(result depends on d in a non-polynomial way).",
                    ),
                )
            end
        end

        if length(expr.cycles) != 1
            return nothing
        end

        factors = expr.cycles[1]

        if length(factors) == 4
            for i = 1:4
                shifted = circshift(factors, -i+1)
                U_cand = shifted[1]
                U_dag_cand = shifted[3]

                if !(U_cand isa SymbolicMatrix) || U_cand.special_type != :U
                    continue
                end

                if !isequal(U_cand.dim, measure.dim)
                    continue
                end

                if !(U_dag_cand isa SymbolicMatrix) ||
                   U_dag_cand.special_type != :U ||
                   U_dag_cand.name != U_cand.name ||
                   U_dag_cand.is_adj == U_cand.is_adj
                    continue
                end

                A = shifted[2]
                B = shifted[4]

                if (A isa SymbolicMatrix && A.name == U_cand.name) ||
                   (B isa SymbolicMatrix && B.name == U_cand.name)
                    continue
                end

                return prefactor * (tr_val([A]) * tr_val([B])) / measure.dim
            end
        end
    end
    return nothing
end


function _collect_scalar_monomial_indices(expr, type_tag::Symbol)
    matcher = MetadataMatcher(type_tag)
    coeff = 1 // 1
    u_indices = Tuple{Any,Any}[]
    u_bar_indices = Tuple{Any,Any}[]

    function _push_match!(match_res, conjugated::Bool)
        tag, i, j = match_res
        if i === nothing || j === nothing
            return false
        end
        is_bar = endswith(String(tag), "_bar")
        final_is_bar = is_bar ⊻ conjugated
        push!(final_is_bar ? u_bar_indices : u_indices, (i, j))
        return true
    end

    function traverse(t, conjugated::Bool = false)
        t_un = Symbolics.unwrap(t)

        m = match_index(matcher, t_un)
        if m !== nothing
            return _push_match!(m, conjugated)
        end

        if t_un isa Number
            coeff *= conjugated ? conj(t_un) : t_un
            return true
        end

        if Symbolics.iscall(t_un)
            op = Symbolics.operation(t_un)
            args = Symbolics.arguments(t_un)

            if op == (*)
                for arg in args
                    traverse(arg, conjugated) || return false
                end
                return true
            elseif op == (/)
                traverse(args[1], conjugated) || return false
                denom = conjugated ? conj(args[2]) : args[2]
                coeff /= denom
                return true
            elseif op == (^)
                base = args[1]
                p = _try_extract_int(Symbolics.unwrap(args[2]))
                if !(p isa Integer) || p < 0
                    return false
                end

                base_un = Symbolics.unwrap(base)
                if Symbolics.iscall(base_un)
                    bop = Symbolics.operation(base_un)
                    bargs = Symbolics.arguments(base_un)
                    if (bop == abs || bop == Base.abs) && iseven(p)
                        k = div(p, 2)
                        x = bargs[1]
                        for _ = 1:k
                            traverse(x, conjugated) || return false
                            traverse(x, !conjugated) || return false
                        end
                        return true
                    elseif bop == abs2 || bop == Base.abs2
                        x = bargs[1]
                        for _ = 1:p
                            traverse(x, conjugated) || return false
                            traverse(x, !conjugated) || return false
                        end
                        return true
                    end
                end

                for _ = 1:p
                    traverse(base, conjugated) || return false
                end
                return true
            elseif op == conj || op == Base.conj
                return traverse(args[1], !conjugated)
            elseif op == (-)
                if length(args) == 1
                    coeff *= -1
                    return traverse(args[1], conjugated)
                end
                return false
            elseif op == (+) ||
                   op == abs ||
                   op == Base.abs ||
                   op == abs2 ||
                   op == Base.abs2 ||
                   op == real ||
                   op == Base.real ||
                   op == imag ||
                   op == Base.imag ||
                   op == hypot ||
                   op == Base.hypot ||
                   op == complex ||
                   op == Base.complex
                return false
            end
        end

        coeff *= conjugated ? conj(t) : t
        return true
    end

    ok = traverse(expr, false)
    if !ok
        return nothing
    end
    return coeff, u_indices, u_bar_indices
end

function _pair_counts(indices::AbstractVector{<:Tuple})
    counts = Dict{Tuple{Any,Any},Int}()
    for idx in indices
        counts[idx] = get(counts, idx, 0) + 1
    end
    return counts
end

function _axis_counts(pair_counts::Dict{Tuple{Any,Any},Int}, axis::Int)
    counts = Dict{Any,Int}()
    for ((i, j), c) in pair_counts
        key = axis == 1 ? i : j
        counts[key] = get(counts, key, 0) + c
    end
    return counts
end

function _all_same_row(pair_counts::Dict{Tuple{Any,Any},Int})
    rows = Set(first(idx) for idx in keys(pair_counts))
    return length(rows) == 1
end

function _match_orthogonal_pattern(pair_counts::Dict{Tuple{Any,Any},Int})
    mults = sort(collect(values(pair_counts)))
    if length(pair_counts) == 1 && (mults[1] == 2 || mults[1] == 4)
        return true
    end

    if length(pair_counts) == 3 && mults == [2, 4, 6]
        if _all_same_row(pair_counts)
            return true
        end
        row_mults = sort(collect(values(_axis_counts(pair_counts, 1))))
        col_mults = sort(collect(values(_axis_counts(pair_counts, 2))))
        return row_mults == [4, 8] && col_mults == [2, 4, 6]
    end

    return false
end

function _match_symplectic_pattern(pair_counts::Dict{Tuple{Any,Any},Int})
    mults = sort(collect(values(pair_counts)))
    if length(pair_counts) == 1 && (mults[1] == 2 || mults[1] == 4)
        return true
    end
    if length(pair_counts) == 2 && mults == [2, 2]
        return true
    end
    return length(pair_counts) == 3 && mults == [2, 4, 6] && _all_same_row(pair_counts)
end

function _match_coe_pattern(pair_counts::Dict{Tuple{Any,Any},Int})
    mults = sort(collect(values(pair_counts)))
    if length(pair_counts) == 1
        return mults[1] == 2 || mults[1] == 4
    end
    return length(pair_counts) == 2 && mults == [2, 2]
end

function _match_cse_pattern(pair_counts::Dict{Tuple{Any,Any},Int})
    mults = sort(collect(values(pair_counts)))
    return length(pair_counts) == 1 && (mults[1] == 2 || mults[1] == 4)
end

function _check_scalar_rule_library(expr, measure, tag::Symbol, pattern_matcher)
    parsed = _collect_scalar_monomial_indices(expr, tag)
    parsed === nothing && return nothing
    coeff, u_indices, u_bar_indices = parsed
    if isempty(u_indices) && isempty(u_bar_indices)
        return nothing
    end

    pair_counts = _pair_counts(vcat(u_indices, u_bar_indices))
    pattern_matcher(pair_counts) || return nothing

    val = INTEGRATION_RULES[tag](u_indices, u_bar_indices, measure.dim, measure.dim)
    return coeff * val
end

function check_orthogonal_library(expr, measure::OrthogonalMeasure)
    return _check_scalar_rule_library(expr, measure, :O, _match_orthogonal_pattern)
end

function check_symplectic_library(expr, measure::SymplecticMeasure)
    return _check_scalar_rule_library(expr, measure, :Sp, _match_symplectic_pattern)
end

function check_coe_library(expr, measure::COEMeasure)
    return _check_scalar_rule_library(expr, measure, :COE, _match_coe_pattern)
end

function check_cse_library(expr, measure::CSEMeasure)
    return _check_scalar_rule_library(expr, measure, :CSE, _match_cse_pattern)
end

function _check_gaussian_scalar_second_moment(expr, measure, type::Symbol)
    parsed = _collect_scalar_monomial_indices(expr, type)
    parsed === nothing && return nothing
    coeff, u_indices, u_bar_indices = parsed
    n = length(u_indices) + length(u_bar_indices)
    n == 2 || return nothing
    val = INTEGRATION_RULES[type](u_indices, u_bar_indices, measure.dim, measure.dim)
    return coeff * val
end

function _same_tagged_name(f::SymbolicMatrix, g::SymbolicMatrix, tag::Symbol)
    return f.special_type == tag && g.special_type == tag && f.name == g.name
end

function _is_adjoint_pair(f::SymbolicMatrix, g::SymbolicMatrix, tag::Symbol)
    _same_tagged_name(f, g, tag) || return false
    return (f.is_adj != g.is_adj) && (f.is_trans != g.is_trans)
end

function _is_transpose_pair(f::SymbolicMatrix, g::SymbolicMatrix, tag::Symbol)
    _same_tagged_name(f, g, tag) || return false
    return (f.is_adj == g.is_adj) && (f.is_trans != g.is_trans)
end

function _is_len2_cycle_pair(cycle, rel)
    length(cycle) == 2 || return false
    a, b = cycle
    return (a isa SymbolicMatrix) && (b isa SymbolicMatrix) && rel(a, b)
end

function _is_len4_alternating_cycle(cycle, rel)
    length(cycle) == 4 || return false
    for i = 1:4
        a = cycle[i]
        b = cycle[mod1(i + 1, 4)]
        if !(a isa SymbolicMatrix && b isa SymbolicMatrix && rel(a, b))
            return false
        end
    end
    return true
end

function check_ginibre_library(expr, measure, type::Symbol)
    if !(expr isa LazyTrace)
        return nothing
    end

    prefactor = expr.prefactor
    d = measure.dim
    cycles = expr.cycles

    if type == :GinUE
        rel = (a, b) -> _is_adjoint_pair(a, b, :GinUE)
        if length(cycles) == 1
            c = cycles[1]
            if _is_len2_cycle_pair(c, rel)
                return prefactor * d^2
            elseif _is_len4_alternating_cycle(c, rel)
                return prefactor * (2d^3)
            end
        elseif length(cycles) == 2 && all(c -> _is_len2_cycle_pair(c, rel), cycles)
            return prefactor * (d^4 + d^2)
        end
    elseif type == :GinOE
        rel = (a, b) -> _is_transpose_pair(a, b, :GinOE)
        if length(cycles) == 1 && _is_len2_cycle_pair(cycles[1], rel)
            return prefactor * d^2
        end
    elseif type == :GinSE
        rel = (a, b) -> _is_adjoint_pair(a, b, :GinSE)
        if length(cycles) == 1 && _is_len2_cycle_pair(cycles[1], rel)
            return prefactor * d^2
        end
    end

    return nothing
end

function check_gaussian_library(expr, measure, type)
    if expr isa LazyTrace
        if length(expr.cycles) != 1
            return nothing
        end

        factors = expr.cycles[1]
        prefactor = expr.prefactor

        expected_tag = type

        if !all(f -> (f isa SymbolicMatrix) && f.special_type == expected_tag, factors)
            return nothing
        end

        k = length(factors)
        d = measure.dim

        val = nothing
        if type == :GUE
            if k == 2
                val = d^2
            elseif k == 4
                val = 2d^3 + d
            elseif k == 6
                val = 5d^4 + 10d^2
            end
        elseif type == :GOE
            if k == 2
                val = d^2 + d
            elseif k == 4
                val = 2d^3 + 5d^2 + 5d
            end
        elseif type == :GSE
            if k == 2
                val = d^2 - d
            elseif k == 4
                val = 2d^3 - 5d^2 + 5d
            end
        end

        if val !== nothing
            return prefactor * val
        end
    end

    return _check_gaussian_scalar_second_moment(expr, measure, type)
end


function check_pure_library(expr, measure)
    return nothing
end
