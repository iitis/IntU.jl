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
        factors = expr.factors
        if length(factors) == 4
            # Normalize cycle
            # We want to find a cyclic shift that is [U, A, U', B]
            # where factors[1] == measure.U and factors[3] == measure.U'
            for i in 1:4
                shifted = circshift(factors, -i+1)
                if shifted[1].special_type == :U && shifted[1].name == (measure.U isa SymbolicMatrix ? measure.U.name : :nothing) &&
                   shifted[3].special_type == :U_dag && shifted[3].name == (measure.U isa SymbolicMatrix ? measure.U.name : :nothing)
                   
                   A = shifted[2]
                   B = shifted[4]
                   return (tr_val([A]) * tr_val([B])) / measure.dim
                end
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
    
    factors = expr.factors
    H_name = measure.H isa SymbolicMatrix ? measure.H.name : :H
    
    # Check if all factors are H
    if !all(f -> f.name == H_name, factors)
        return nothing
    end
    
    k = length(factors) # tr(H^k)
    d = measure.dim
    
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
            return d^2 - d
        elseif k == 4
            return 2d^3 - 5d^2 + 5d
        end
    end
    
    return nothing
end

# --- Pure States Library ---

function check_pure_library(expr, measure)
    # Most pure state integrals are on indices, but we could add matrix-level later
    # For now, if someone does tr(|psi><psi|), it's 1.
    return nothing
end
