# src/library.jl

"""
    check_library(expr, measure)

Check if the integral of `expr` over `measure` is available in the pre-computed library.
Returns the symbolic result if found, otherwise `nothing`.
"""
function check_library(expr, measure)
    if measure isa HaarMeasure
        return check_haar_library(expr, measure)
    elseif measure isa GUEMeasure
        return check_gaussian_library(expr, measure, :GUE)
    elseif measure isa GOEMeasure
        return check_gaussian_library(expr, measure, :GOE)
    elseif measure isa GSEMeasure
        return check_gaussian_library(expr, measure, :GSE)
    elseif measure isa PureStateMeasure
        return check_pure_library(expr, measure)
    end
    return nothing
end

# --- Haar Unitary Library ---

function check_haar_library(expr, measure)
    # Check for tr(U A U' B) where U is the integration variable
    # This matches LazyTrace of [U, A, U', B]
    if expr isa LazyTrace
        if length(expr.cycles) != 1
            return nothing
        end

        factors = expr.cycles[1]
        prefactor = expr.prefactor

        if length(factors) == 4
            # Normalize cycle
            # We want to find a cyclic shift that is [U, A, U', B]
            for i = 1:4
                shifted = circshift(factors, -i+1)
                U_cand = shifted[1]
                U_dag_cand = shifted[3]
                
                # Check 1: U_cand is a Unitary (:U)
                if U_cand.special_type != :U
                    continue
                end
                
                # Check 2: Dimensions match measure
                if !isequal(U_cand.dim, measure.dim)
                    continue
                end

                # Check 3: U_dag_cand matches U_cand name and is adjoint pair
                # We require them to be the same matrix, but one adjoint and one not (or opposite adjointness)
                if U_dag_cand.special_type != :U ||
                   U_dag_cand.name != U_cand.name ||
                   U_dag_cand.is_adj == U_cand.is_adj
                   continue
                end
                
                # Check 4: A and B do not depend on U
                A = shifted[2]
                B = shifted[4]
                
                if A.name == U_cand.name || B.name == U_cand.name
                    continue
                end
                
                return prefactor * (tr_val([A]) * tr_val([B])) / measure.dim
            end
        end
    end
    return nothing
end

# --- Gaussian Ensembles Library ---

function check_gaussian_library(expr, measure, type)
    if !(expr isa LazyTrace)
        return nothing
    end

    if length(expr.cycles) != 1
        return nothing
    end

    factors = expr.cycles[1]
    prefactor = expr.prefactor

    # For the library, we check if all factors are of the expected special_type
    expected_tag = (type == :GUE || type == :GOE || type == :GSE) ? :H : :G
    
    if !all(f -> f.special_type == expected_tag, factors)
        return nothing
    end

    k = length(factors) # tr(H^k)
    d = measure.dim

    val = nothing
    if type == :GUE
        if k == 2
            return d^2
        elseif k == 4
            return 2d^3 + d
        elseif k == 6
            return 5d^4 + 10d^2
        end
    elseif type == :GOE
        if k == 2
            return d^2 + d
        elseif k == 4
            return 2d^3 + 5d^2 + 5d
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

# --- Pure States Library ---

function check_pure_library(expr, measure)
    # Most pure state integrals are on indices, but we could add matrix-level later
    # For now, if someone does tr(|psi><psi|), it's 1.
    return nothing
end
