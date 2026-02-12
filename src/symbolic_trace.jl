# src/symbolic_trace.jl

"""
    SymbolicMatrix(name::Symbol)
    SymbolicMatrix(name::Symbol, special_type::Symbol)
    SymbolicMatrix(name::Symbol, special_type::Symbol, dim)
    SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol, dim)

A wrapper associated with a symbolic name to represent a matrix in a coordinate-free way.
Used in the symbolic trace logic (via `tr_lazy`) and for metadata-driven element-wise integration.

# Special Types
The `special_type` field determines how the matrix entries are tagged via metadata for the integration engine:
- `:Constant` (default): Treated as a fixed constant matrix.
- `:U`: Marks the matrix as a Haar-random unitary (used by `dU`, `dSU`, `dStiefel`, `dCUE`).
- `:O`: Marks the matrix as a Haar-random real orthogonal matrix (used by `dO`, `dGOE`, `dCOE`).
- `:Sp`: Marks the matrix as a Haar-random symplectic matrix (used by `dSp`, `dGSE`, `dCSE`).
- `:Perm`: Marks the matrix as a random permutation matrix (used by `dPerm`, `dCPerm`).
- `:DiagUnitary`: Marks the matrix as a random diagonal unitary matrix (used by `dDiagUnitary`).
- `:GUE` / `:GOE` / `:GSE`: Marks the matrix as a member of the respective Gaussian ensemble.
"""
struct SymbolicMatrix <: AbstractMatrix{Any}
    name::Symbol
    is_adj::Bool
    special_type::Symbol 
    dim::Any

    function SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol, dim::Any)
        new(name, is_adj, special_type, dim)
    end
end

"""
    MatrixMetadata

Internal type used as a key for Symbolics/SymbolicUtils metadata.
"""
struct MatrixMetadata end

SymbolicMatrix(name::Symbol) = SymbolicMatrix(name, false, :Constant, nothing)
SymbolicMatrix(name::Symbol, special_type::Symbol) = SymbolicMatrix(name, false, special_type, nothing)
SymbolicMatrix(name::Symbol, special_type::Symbol, dim) = SymbolicMatrix(name, false, special_type, dim)

import Base: *, adjoint, show, ^, size, getindex
import LinearAlgebra: tr

# AbstractMatrix implementation
function Base.size(A::SymbolicMatrix)
    d = A.dim === nothing ? typemax(Int) : A.dim
    return (d, d)
end
function Base.getindex(A::SymbolicMatrix, i::Integer, j::Integer)
    # Generate symbol name
    s_name = Symbol(A.name, :_, i, :_, j)
    if A.is_adj
        s_name = Symbol(A.name, :_, j, :_, i)
    end
    
    # Create variable and attach metadata
    # We use T=Number to prevent it from being treated as Real and simplifying conj(v) -> v
    # This keeps v as an opaque complex-like symbol.
    v = Symbolics.variable(s_name, T = Complex{Real})
    
    meta = Dict(
        :name => A.name,
        :type => A.special_type,
        :indices => (A.is_adj ? (j, i) : (i, j)),
        :is_adj => A.is_adj
    )
    
    v = Num(Symbolics.setmetadata(v, MatrixMetadata, meta))
    
    return A.is_adj ? conj(v) : v
end

function adjoint(A::SymbolicMatrix)
    return SymbolicMatrix(A.name, !A.is_adj, A.special_type, A.dim)
end

function Base.show(io::IO, A::SymbolicMatrix)
    print(io, A.name)
    if A.is_adj
        print(io, "'")
    end
end

function Base.show(io::IO, ::MIME"text/plain", A::SymbolicMatrix)
    adj_str = A.is_adj ? "'" : ""
    print(io, "SymbolicMatrix (type=$(A.special_type)): $(A.name)$adj_str")
end

"""
    LazyTrace

Represents the product of traces of sequences of SymbolicMatrices.
e.g. `tr(A B) * tr(C)` is represented as one LazyTrace with two cycles.
"""
struct LazyTrace
    cycles::Vector{Vector{SymbolicMatrix}} # list of cycles. Each cycle is a product trace.
    prefactor::Any                         # Scalar prefactor (Number or Num)
end

"""
    LazySum

Represents a sum of LazyTraces.
"""
struct LazySum
    terms::Vector{LazyTrace}
end

"""
    SymbolicMatrixProduct

Represents a product of SymbolicMatrices.
"""
struct SymbolicMatrixProduct <: AbstractMatrix{Num}
    factors::Vector{Any} # Can be SymbolicMatrix, SymbolicKron, or AbstractMatrix
end

function Base.size(P::SymbolicMatrixProduct)
    if isempty(P.factors)
        return (0, 0)
    end
    return (size(P.factors[1], 1), size(P.factors[end], 2))
end

function Base.getindex(P::SymbolicMatrixProduct, i::Integer, j::Integer)
    # Recursively handle products
    # For a product A * B, (AB)_ij = sum_k A_ik B_kj
    # We only support small fixed dimensions for explicit getindex if not using tr_lazy
    factors = P.factors
    if length(factors) == 1
        return factors[1][i, j]
    end
    
    A = factors[1]
    # Rest of the product
    B = length(factors) == 2 ? factors[2] : SymbolicMatrixProduct(factors[2:end])
    
    dim = size(A, 2)
    if dim isa Integer
        return sum(A[i, k] * B[k, j] for k = 1:dim)
    else
        # Symbolic dimension sum. For integration, we might need a better way.
        # But let's try to return a symbolic sum.
        k = Symbolics.variable(:k, T = Int)
        # return Symbolics.sum(A[i, k] * B[k, j], k => 1:dim) # Symbolics.sum is still experimental or in SymbolicUtils
        # For now, let's use a dummy name to avoid infinite sums if it's evaluated
        return Num(Symbolics.variable(Symbol("sum_$(A)_$(B)_$(i)_$(j)"); T = Number))
    end
end

function *(A::SymbolicMatrix, B::SymbolicMatrix)
    return SymbolicMatrixProduct([A, B])
end

function *(A::SymbolicMatrixProduct, B::SymbolicMatrix)
    return SymbolicMatrixProduct(vcat(A.factors, [B]))
end

function *(A::SymbolicMatrix, B::SymbolicMatrixProduct)
    return SymbolicMatrixProduct(vcat([A], B.factors))
end

function *(A::SymbolicMatrixProduct, B::SymbolicMatrixProduct)
    return SymbolicMatrixProduct(vcat(A.factors, B.factors))
end

function ^(A::SymbolicMatrix, n::Integer)
    return SymbolicMatrixProduct(fill(A, n))
end

"""
    SymbolicKron
    
Represents the tensor product of two symbolic matrices.
"""
struct SymbolicKron <: AbstractMatrix{Num}
    A::Any # SymbolicMatrix or SymbolicKron or SymbolicMatrixProduct
    B::Any
end

function Base.size(K::SymbolicKron)
    sA = size(K.A)
    sB = size(K.B)
    return (sA[1] * sB[1], sA[2] * sB[2])
end

function Base.getindex(K::SymbolicKron, i::Integer, j::Integer)
    mA, nA = size(K.A)
    mB, nB = size(K.B)
    iA, iB = divrem(i - 1, mB) .+ 1
    jA, jB = divrem(j - 1, nB) .+ 1
    return K.A[iA, jA] * K.B[iB, jB]
end

function Base.kron(A::SymbolicMatrix, B::SymbolicMatrix)
    return SymbolicKron(A, B)
end

function Base.kron(A::SymbolicMatrix, B::AbstractMatrix)
    return SymbolicKron(A, B)
end

function Base.kron(A::AbstractMatrix, B::SymbolicMatrix)
    return SymbolicKron(A, B)
end

function adjoint(P::SymbolicMatrixProduct)
    return SymbolicMatrixProduct(reverse([adjoint(f) for f in P.factors]))
end

function adjoint(K::SymbolicKron)
    return SymbolicKron(adjoint(K.A), adjoint(K.B))
end

function tr(A::SymbolicMatrix)
    return tr_lazy(A)
end

function tr(A::SymbolicMatrixProduct)
    return tr_lazy(A.factors)
end

function tr(A::Symbolics.Arr{T,2}) where {T}
    return sum(A[i, i] for i = 1:size(A, 1))
end

"""
    tr_lazy(product)

Create a LazyTrace from a product of SymbolicMatrices (or a single one).
"""
function tr_lazy(product::AbstractVector)
    return LazyTrace([collect(SymbolicMatrix, product)], 1)
end

function tr_lazy(product::SymbolicMatrix)
    return LazyTrace([[product]], 1)
end

function tr_lazy(product::SymbolicMatrixProduct)
    return LazyTrace([product.factors], 1)
end

# Arithmetic Operations

# Multiplication: LazyTrace * LazyTrace -> LazyTrace (merge cycles)
function Base.:*(a::LazyTrace, b::LazyTrace)
    return LazyTrace(vcat(a.cycles, b.cycles), a.prefactor * b.prefactor)
end

# Multiplication: LazyTrace * Number -> LazyTrace
function Base.:*(a::LazyTrace, b::Number)
    return LazyTrace(a.cycles, a.prefactor * b)
end
function Base.:*(b::Number, a::LazyTrace)
    return LazyTrace(a.cycles, a.prefactor * b)
end

# Multiplication: LazyTrace * Num -> LazyTrace
function Base.:*(a::LazyTrace, b::Num)
    return LazyTrace(a.cycles, a.prefactor * b)
end
function Base.:*(b::Num, a::LazyTrace)
    return LazyTrace(a.cycles, a.prefactor * b)
end

# Complex Conjugation: conj(LazyTrace)
# conj(Tr(A B ...)) = Tr( (A B ...)' ) = Tr( ...' B' A' )
function Base.conj(t::LazyTrace)
    new_cycles = Vector{Vector{SymbolicMatrix}}()
    for cycle in t.cycles
        # Reverse order and adjoint each factor
        new_cycle = reverse([adjoint(f) for f in cycle])
        push!(new_cycles, new_cycle)
    end
    return LazyTrace(new_cycles, Symbolics.conj(t.prefactor))
end

# Absolute value squared: |Tr(U)|^2 = Tr(U) * conj(Tr(U))
function Base.abs2(t::LazyTrace)
    return t * conj(t)
end

# Addition: LazyTrace + LazyTrace -> LazySum
function Base.:+(a::LazyTrace, b::LazyTrace)
    return LazySum([a, b])
end

# Addition: LazySum + LazyTrace -> LazySum
function Base.:+(a::LazySum, b::LazyTrace)
    return LazySum(vcat(a.terms, b))
end
function Base.:+(b::LazyTrace, a::LazySum)
    return LazySum(vcat(b, a.terms))
end

# Addition: LazySum + LazySum -> LazySum
function Base.:+(a::LazySum, b::LazySum)
    return LazySum(vcat(a.terms, b.terms))
end

# Pow: LazyTrace ^ Integer
function Base.:^(a::LazyTrace, n::Integer)
    if n == 0
        return LazyTrace([], 1)
    end
    if n == 1
        return a
    end
    # Repeat cycles n times
    new_cycles = Vector{Vector{SymbolicMatrix}}()
    for _ = 1:n
        append!(new_cycles, a.cycles)
    end
    return LazyTrace(new_cycles, a.prefactor^n)
end

# Show methods
function show(io::IO, t::LazyTrace)
    if t.prefactor != 1
        print(io, t.prefactor, "*")
    end
    if isempty(t.cycles)
        print(io, "1")
        return
    end
    for (i, cycle) in enumerate(t.cycles)
        print(io, "tr(")
        for (j, f) in enumerate(cycle)
            print(io, f)
            if j < length(cycle)
                print(io, " * ")
            end
        end
        print(io, ")")
    end
end

function show(io::IO, s::LazySum)
    for (i, t) in enumerate(s.terms)
        print(io, t)
        if i < length(s.terms)
            print(io, " + ")
        end
    end
end




"""
    tr_val(p::SymbolicMatrixProduct)

Symbolic function representing the trace of a matrix product.
"""
function tr_val(factors::Vector{SymbolicMatrix})
    # Create a symbolic term representing this trace
    # This avoids issues with Term multiplication and simplification
    if isempty(factors)
        return 1
    end

    s_parts = String[]
    for (i, f) in enumerate(factors)
        push!(s_parts, string(f))
    end

    inner_content_name = join(s_parts, "*")

    name = "tr(" * inner_content_name * ")"
    return Num(Symbolics.variable(Symbol(name); T = Real))
end
