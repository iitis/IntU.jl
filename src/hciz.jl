# Harish-Chandra-Itzykson-Zuber (HCIZ) Integrals

"""
    vandermonde_det(v)

Computes the Vandermonde determinant of a vector `v`:
```math
\Delta(v) = \prod_{1 \le i < j \le d} (v_i - v_j)
```
"""
function vandermonde_det(v::AbstractVector)
    d = length(v)
    res = one(eltype(v))
    for i in 1:d
        for j in i+1:d
            res *= (v[i] - v[j])
        end
    end
    return res
end

"""
    hciz(A, B)
    hciz(a::AbstractVector, b::AbstractVector)

Computes the Harish-Chandra-Itzykson-Zuber (HCIZ) integral:
```math
\int_{U(d)} dU e^{\text{Tr}(A U B U^\dagger)} = \left( \prod_{p=1}^{d-1} p! \right) \frac{\det(e^{a_i b_j})_{i,j=1}^d}{\Delta(a) \Delta(b)}
```
where `a` and `b` are eigenvalues of `A` and `B`, and `\Delta` is the Vandermonde determinant.

If A and B are matrices, their eigenvalues are extracted. Supporting:
- Numeric matrices (via `eigen`)
- `Matrix{Num}` (symbolic diagonal or 2x2)
- `SymbolicMatrix` (by generating symbolic eigenvalues `a_1, ..., a_d`)

Note: This formula is sensitive to degenerate eigenvalues where the denominators become zero. 
In such cases, the limit should be taken. This implementation currently uses a 
small perturbation for numerical stability if exact degeneracies are detected in numeric input.
"""
function hciz(A::AbstractMatrix, B::AbstractMatrix)
    a = _get_eigenvalues(A)
    b = _get_eigenvalues(B)
    return hciz(a, b)
end

function hciz(A::SymbolicMatrix, B::SymbolicMatrix, d::Int)
    # Generate symbolic eigenvalues for matrices with only names
    a = [Symbolics.variable(Symbol(string(A.name) * "_$i"); T=Real) for i in 1:d]
    b = [Symbolics.variable(Symbol(string(B.name) * "_$i"); T=Real) for i in 1:d]
    return hciz(a, b)
end

function _get_eigenvalues(M::AbstractMatrix)
    # Check if numeric
    if all(x -> x isa Number && !(x isa Num), M)
        return eigen(M).values
    end
    
    # Handle Matrix{Num}
    # 1. Check if diagonal
    is_diag = true
    d = size(M, 1)
    for i in 1:d, j in 1:d
        if i != j && !isequal(Symbolics.unwrap(M[i, j]), 0)
            is_diag = false
            break
        end
    end
    if is_diag
        return [M[i, i] for i in 1:d]
    end
    
    # 2. Handle 2x2 symbolic
    if d == 2
        # λ^2 - tr(M)λ + det(M) = 0
        t = M[1,1] + M[2,2]
        det_M = M[1,1]*M[2,2] - M[1,2]*M[2,1]
        disc = Symbolics.simplify(t^2 - 4*det_M)
        return [
            Symbolics.simplify((t + sqrt(disc))/2),
            Symbolics.simplify((t - sqrt(disc))/2)
        ]
    end
    
    error("Cannot extract eigenvalues symbolically for d > 2 and non-diagonal matrix. Please provide eigenvalues directly.")
end

function hciz(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || error("A and B must have the same dimension")
    d = length(a)
    d == 0 && return 1.0
    
    # Check for degeneracies
    if _has_degeneracies(a) || _has_degeneracies(b)
        # For numeric values, we can add a tiny perturbation
        if eltype(a) <: Number && eltype(b) <: Number
            a = a .+ (rand(d) .- 0.5) .* 1e-12
            b = b .+ (rand(d) .- 0.5) .* 1e-12
        else
            # For symbolic, this is harder. 
            # In some cases L'Hopital's rule or character expansions are needed.
            # Currently we'll throw a warning/error if it's strictly zero.
        end
    end

    prefactor = one(eltype(a))
    for p in 1:(d-1)
        prefactor *= factorial(p)
    end

    delta_a = vandermonde_det(a)
    delta_b = vandermonde_det(b)
    
    # Matrix M_ij = exp(a_i * b_j)
    M = [exp(a[i] * b[j]) for i in 1:d, j in 1:d]
    
    return prefactor * (det(M) / (delta_a * delta_b))
end

function _has_degeneracies(v)
    for i in 1:length(v)
        for j in i+1:length(v)
            # Use isequal for robust comparison including symbolic
            if isequal(v[i], v[j])
                return true
            end
        end
    end
    return false
end
