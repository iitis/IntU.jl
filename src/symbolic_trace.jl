# src/symbolic_trace.jl

"""
    SymbolicMatrix(name::Symbol)
    SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol)

A wrapper associated with a symbolic name to represent a matrix in a coordinate-free way.
`special_type` can be `:U`, `:U_dag`, or `:Constant`.
"""
struct SymbolicMatrix
    name::Symbol
    is_adj::Bool
    special_type::Symbol  # :U, :U_dag, :Constant
end

SymbolicMatrix(name::Symbol) = SymbolicMatrix(name, false, :Constant)

import Base: *, adjoint, show, ^

function adjoint(A::SymbolicMatrix)
    new_type = A.special_type
    if A.special_type == :U
        new_type = :U_dag
    elseif A.special_type == :U_dag
        new_type = :U
    end
    # For constants, new_type remains :Constant (or we could have :Constant_dag, but :Constant handles both usually)
    return SymbolicMatrix(A.name, !A.is_adj, new_type)
end

function show(io::IO, A::SymbolicMatrix)
    print(io, A.name)
    if A.is_adj
        print(io, "'")
    end
end

"""
    LazyTrace

Represents the trace of a sequence of SymbolicMatrices.
"""
struct LazyTrace
    factors::Vector{SymbolicMatrix}
end

function *(A::SymbolicMatrix, B::SymbolicMatrix)
    return [A, B]
end

function *(A::Vector{SymbolicMatrix}, B::SymbolicMatrix)
    return push!(copy(A), B)
end

function *(A::SymbolicMatrix, B::Vector{SymbolicMatrix})
    return vcat([A], B)
end

function *(A::Vector{SymbolicMatrix}, B::Vector{SymbolicMatrix})
    return vcat(A, B)
end

function ^(A::SymbolicMatrix, n::Integer)
    return fill(A, n)
end

function tr(A::SymbolicMatrix)
    return tr_lazy(A)
end

function tr(A::Vector{SymbolicMatrix})
    return tr_lazy(A)
end

"""
    tr_lazy(product)

Create a LazyTrace from a product of SymbolicMatrices (or a single one).
"""
function tr_lazy(product::AbstractVector)
    return LazyTrace(collect(SymbolicMatrix, product))
end

function tr_lazy(product::SymbolicMatrix)
    return LazyTrace([product])
end

function show(io::IO, t::LazyTrace)
    print(io, "tr(")
    for (i, f) in enumerate(t.factors)
        print(io, f)
        if i < length(t.factors)
            print(io, " * ")
        end
    end
    print(io, ")")
end




"""
    tr_val(p::SymbolicMatrixProduct)

Symbolic function representing the trace of a matrix product.
"""
function tr_val(factors::Vector{SymbolicMatrix})
    # Create a symbolic variable representing this trace
    # This avoids issues with Term multiplication and simplification
    if isempty(factors)
        return 1 # Should handle dim separately, but trace of empty is not passed here usually
    end
    
    # Construct a nice string representation
    # e.g. "tr(A * B)"
    s_parts = String[]
    for (i, f) in enumerate(factors)
        push!(s_parts, string(f))
    end
    name = "tr_val(" * join(s_parts, "*") * ")"
    return Symbolics.variable(Symbol(name); T=Real)
end
# Symbolics metadata might go here if needed.
# Actually, we rely on Term wrapping it.

