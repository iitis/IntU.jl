# src/symbolic_trace.jl

"""
    SymbolicMatrix(name::Symbol)
    SymbolicMatrix(name::Symbol, special_type::Symbol)
    SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol)

A wrapper associated with a symbolic name to represent a matrix in a coordinate-free way.
Used in the symbolic trace logic (via `tr_lazy`).

# Special Types
The `special_type` field determines how the matrix is handled during trace-based integration:
- `:Constant` (default): Treated as a fixed constant matrix.
- `:U`: Marks the matrix as a Haar-random unitary for integration under `dU`.
- `:U_dag`: Marks the matrix as the adjoint of a Haar-random unitary.

For Gaussian measures (GUE, GinUE, etc.), matrices are primarily identified by their name, 
but using `:Constant` for fixed matrices is recommended.
"""
struct SymbolicMatrix <: AbstractMatrix{Num}
    name::Symbol
    is_adj::Bool
    special_type::Symbol  # :U, :U_dag, :Constant
end

SymbolicMatrix(name::Symbol) = SymbolicMatrix(name, false, :Constant)
SymbolicMatrix(name::Symbol, special_type::Symbol) = SymbolicMatrix(name, false, special_type)

import Base: *, adjoint, show, ^, size, getindex
import LinearAlgebra: tr

# AbstractMatrix implementation
Base.size(::SymbolicMatrix) = (typemax(Int), typemax(Int))
function Base.getindex(A::SymbolicMatrix, i::Integer, j::Integer)
    if A.is_adj
        return conj(Symbolics.variable(Symbol(A.name, :_, j, :_, i), T = Complex{Num}))
    end
    return Symbolics.variable(Symbol(A.name, :_, i, :_, j), T = Complex{Num})
end

function adjoint(A::SymbolicMatrix)
    new_type = A.special_type
    if A.special_type == :U
        new_type = :U_dag
    elseif A.special_type == :U_dag
        new_type = :U
    end
    return SymbolicMatrix(A.name, !A.is_adj, new_type)
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
struct SymbolicMatrixProduct
    factors::Vector{SymbolicMatrix}
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
