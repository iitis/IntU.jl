
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
            moment = d isa Integer ? _haar_trace_moment_value(k, d) : factorial(k)
            return prefactor * moment
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


function check_gaussian_library(expr, measure, type)
    if !(expr isa LazyTrace)
        return nothing
    end

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

    return nothing
end


function check_pure_library(expr, measure)
    return nothing
end
