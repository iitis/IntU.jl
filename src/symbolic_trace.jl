
"""
    SymbolicMatrix(name::Symbol)
    SymbolicMatrix(name::Symbol, special_type::Symbol)
    SymbolicMatrix(name::Symbol, special_type::Symbol, dim)
    SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol, dim)

A wrapper associated with a symbolic name to represent a matrix in a coordinate-free way.
Used in the symbolic trace logic (via `tr_lazy`) and for metadata-driven element-wise integration.
"""
struct SymbolicMatrix <: AbstractMatrix{Num}
    name::Symbol
    is_adj::Bool
    special_type::Symbol 
    dim::Union{Nothing, Integer, Num, Tuple{Union{Integer, Num}, Union{Integer, Num}}, Any}

    function SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol, dim::Any)
        new(name, is_adj, special_type, dim)
    end
end

struct MatrixMetadata end

SymbolicMatrix(name::Symbol) = SymbolicMatrix(name, false, :Constant, nothing)
SymbolicMatrix(name::Symbol, special_type::Symbol) = SymbolicMatrix(name, false, special_type, nothing)
SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol) = SymbolicMatrix(name, is_adj, special_type, nothing)
SymbolicMatrix(name::Symbol, special_type::Symbol, dim) = SymbolicMatrix(name, false, special_type, dim)

import Base: *, adjoint, transpose, show, ^, size, getindex
import LinearAlgebra: tr

function Base.size(A::SymbolicMatrix)
    d = A.dim === nothing ? (typemax(Int), typemax(Int)) : A.dim
    if d isa Tuple
        rows, cols = d
        return A.is_adj ? (cols, rows) : (rows, cols)
    else
        return (d, d)
    end
end

function Base.getindex(A::SymbolicMatrix, i::Integer, j::Integer)
    s_name = Symbol(A.name, :_, i, :_, j)
    if A.is_adj
        s_name = Symbol(A.name, :_, j, :_, i)
    end
    
    meta = Dict(
        :name => A.name,
        :type => A.special_type,
        :indices => (A.is_adj ? (j, i) : (i, j)),
        :is_adj => A.is_adj
    )

    # Use Matrix{Any} in integration to avoid symtype errors
    v = Symbolics.variable(s_name, T = Number)
    v_un = Symbolics.unwrap(v)
    v_meta_un = SymbolicUtils.setmetadata(v_un, MatrixMetadata, meta)
    v_meta = Symbolics.wrap(v_meta_un)
    
    return A.is_adj ? conj(v_meta) : v_meta
end

function Base.getindex(A::SymbolicMatrix, i::Union{Integer, AbstractVector, Colon}, j::Union{Integer, AbstractVector, Colon})
    rows = (i isa Colon) ? (1:size(A, 1)) : i
    cols = (j isa Colon) ? (1:size(A, 2)) : j
    
    # Handle single element access by dispatching back to (Integer, Integer)
    if rows isa Integer && cols isa Integer
        return invoke(getindex, Tuple{SymbolicMatrix, Integer, Integer}, A, rows, cols)
    end
    
    res = Matrix{Any}(undef, length(rows), length(cols))
    for (r_idx, r) in enumerate(rows)
        for (c_idx, c) in enumerate(cols)
            res[r_idx, c_idx] = A[r, c]
        end
    end
    # Return vector if single column/row was requested (standard Julia behavior)
    if length(cols) == 1 && (j isa Integer)
        return res[:, 1]
    elseif length(rows) == 1 && (i isa Integer)
        return res[1, :]
    end
    return SymbolicMatrixProduct([res])
end

function Base.adjoint(A::SymbolicMatrix)
    return SymbolicMatrix(A.name, !A.is_adj, A.special_type, A.dim)
end

function Base.transpose(A::SymbolicMatrix)
    return SymbolicMatrix(A.name, !A.is_adj, A.special_type, A.dim)
end

function Base.show(io::IO, A::SymbolicMatrix)
    print(io, A.name)
    if A.is_adj
        print(io, "'")
    end
end

struct LazyTrace
    cycles::Vector{Vector{AbstractMatrix}}
    prefactor::Union{Num, Number}
end

function Base.getproperty(t::LazyTrace, s::Symbol)
    if s == :factors
        return vcat(t.cycles...)
    end
    return getfield(t, s)
end

struct LazySum
    terms::Vector{LazyTrace}
end

struct SymbolicMatrixProduct <: AbstractMatrix{Num}
    factors::Vector{AbstractMatrix}
end

const SymbolicAny = Union{SymbolicMatrix, SymbolicMatrixProduct}

function Base.size(P::SymbolicMatrixProduct)
    if isempty(P.factors) return (0, 0) end
    return (size(P.factors[1], 1), size(P.factors[end], 2))
end

function Base.size(P::SymbolicMatrixProduct, i::Integer)
    if i == 1
        return size(P.factors[1], 1)
    elseif i == 2
        return size(P.factors[end], 2)
    else
        return 1
    end
end

function Base.getindex(P::SymbolicMatrixProduct, i::Integer, j::Integer)
    factors = P.factors
    if length(factors) == 1 return factors[1][i, j] end
    A = factors[1]
    B = length(factors) == 2 ? factors[2] : SymbolicMatrixProduct(factors[2:end])
    dim = size(A, 2)
    if dim isa Integer
        return Symbolics.wrap(sum(A[i, k] * B[k, j] for k = 1:dim))
    else
        return Num(Symbolics.variable(Symbol("sum_$(A)_$(B)_$(i)_$(j)"); T = Number))
    end
end

# Multiplication logic
_factors(A::SymbolicMatrix) = Any[A]
_factors(P::SymbolicMatrixProduct) = P.factors
_factors(A::AbstractMatrix) = Any[A]

function *(A::SymbolicAny, B::SymbolicAny)
    return SymbolicMatrixProduct(vcat(_factors(A), _factors(B)))
end
function *(A::SymbolicAny, B::AbstractMatrix)
    return SymbolicMatrixProduct(vcat(_factors(A), Any[B]))
end
function *(A::AbstractMatrix, B::SymbolicAny)
    return SymbolicMatrixProduct(vcat(Any[A], _factors(B)))
end

function *(A::SymbolicAny, B::SymbolicAny, Cs::SymbolicAny...)
    res_factors = vcat(_factors(A), _factors(B))
    for C in Cs
        append!(res_factors, _factors(C))
    end
    return SymbolicMatrixProduct(res_factors)
end

# Resolve ambiguities with LinearAlgebra and Symbolics
using LinearAlgebra
for T in [Adjoint, Transpose]
    @eval function *(A::$T{<:Any, <:AbstractVector}, B::SymbolicAny)
        return SymbolicMatrixProduct(vcat(Any[A], _factors(B)))
    end
    @eval function *(A::SymbolicAny, B::$T{<:Any, <:AbstractVector})
        return SymbolicMatrixProduct(vcat(_factors(A), Any[B]))
    end
end

function *(A::Symbolics.Arr, B::SymbolicAny)
    return SymbolicMatrixProduct(vcat(Any[A], _factors(B)))
end
function *(A::SymbolicAny, B::Symbolics.Arr)
    return SymbolicMatrixProduct(vcat(_factors(A), Any[B]))
end

# Resolve specific ambiguities discovered during tests
for T in [Adjoint, Transpose]
    # 3-arg
    @eval function *(A::$T{<:Any, <:AbstractVector}, B::SymbolicAny, C::SymbolicAny)
        return (A * B) * C
    end
    @eval function *(A::$T{<:Any, <:AbstractVector}, B::AbstractMatrix, C::SymbolicAny)
        return (A * B) * C
    end
    @eval function *(A::$T{<:Any, <:AbstractVector}, B::SymbolicAny, C::AbstractMatrix)
        return (A * B) * C
    end
    
    # 4-arg (to prevent LinearAlgebra from catching it)
    @eval function *(A::$T{<:Any, <:AbstractVector}, B::AbstractMatrix, C::AbstractMatrix, D::SymbolicAny)
        return (A * B) * (C * D)
    end
    @eval function *(A::$T{<:Any, <:AbstractVector}, B::AbstractMatrix, C::SymbolicAny, D::AbstractMatrix)
        return (A * B) * (C * D)
    end

    # Disambiguate with internal overlaps
    @eval function *(A::$T{<:Any, <:AbstractVector}, B::AbstractMatrix, C::SymbolicAny, D::SymbolicAny)
        return (A * B) * (C * D)
    end

    # Disambiguate with Symbolics.Arr (3-arg)
    @eval function *(A::$T{W, <:AbstractVector}, B::Symbolics.Arr, C::SymbolicAny) where W
        return (A * B) * C
    end
    @eval function *(A::$T{W, <:AbstractVector}, B::Symbolics.Arr, C::SymbolicAny) where W<:Number
        return (A * B) * C
    end

    # Disambiguate with Symbolics.Arr (4-arg)
    for W_type in [Any, Number]
        @eval function *(A::$T{W, <:AbstractVector}, B::Symbolics.Arr, C::AbstractMatrix, D::SymbolicAny) where W<:$W_type
            return (A * B) * (C * D)
        end
        @eval function *(A::$T{W, <:AbstractVector}, B::Symbolics.Arr, C::SymbolicAny, D::AbstractMatrix) where W<:$W_type
            return (A * B) * (C * D)
        end
        @eval function *(A::$T{W, <:AbstractVector}, B::Symbolics.Arr, C::SymbolicAny, D::SymbolicAny) where W<:$W_type
            return (A * B) * (C * D)
        end
    end
end

function ^(A::SymbolicMatrix, n::Integer)
    return SymbolicMatrixProduct(fill(A, n))
end

function adjoint(P::SymbolicMatrixProduct)
    return SymbolicMatrixProduct(reverse([adjoint(f) for f in P.factors]))
end

function tr(A::SymbolicMatrix)
    return tr_lazy(A)
end
function tr(A::SymbolicMatrixProduct)
    return tr_lazy(A.factors)
end

"""
    tr_lazy(product)

Creates a `LazyTrace` representing the symbolic trace of a matrix product.
The product can be a `SymbolicMatrix`, `SymbolicMatrixProduct`, or a vector of matrices.
"""
function tr_lazy(product::AbstractVector)
    return LazyTrace([collect(Any, product)], 1)
end
function tr_lazy(product::SymbolicMatrix)
    return LazyTrace([[product]], 1)
end
function tr_lazy(product::SymbolicMatrixProduct)
    return LazyTrace([product.factors], 1)
end

function Base.:*(a::LazyTrace, b::LazyTrace)
    return LazyTrace(vcat(a.cycles, b.cycles), a.prefactor * b.prefactor)
end
function Base.:*(a::LazyTrace, b::Number)
    return LazyTrace(a.cycles, a.prefactor * b)
end
function Base.:*(b::Number, a::LazyTrace)
    return LazyTrace(a.cycles, a.prefactor * b)
end
function Base.:*(a::LazyTrace, b::Num)
    return LazyTrace(a.cycles, a.prefactor * b)
end
function Base.:*(b::Num, a::LazyTrace)
    return LazyTrace(a.cycles, a.prefactor * b)
end

function Base.conj(t::LazyTrace)
    new_cycles = Vector{Vector{SymbolicMatrix}}()
    for cycle in t.cycles
        new_cycle = reverse([adjoint(f) for f in cycle])
        push!(new_cycles, new_cycle)
    end
    return LazyTrace(new_cycles, Symbolics.conj(t.prefactor))
end

function Base.abs2(t::LazyTrace)
    return t * conj(t)
end

function Base.:+(a::LazyTrace, b::LazyTrace)
    return LazySum([a, b])
end
function Base.:+(a::LazySum, b::LazyTrace)
    return LazySum(vcat(a.terms, b))
end
function Base.:+(b::LazyTrace, a::LazySum)
    return LazySum(vcat(b, a.terms))
end
function Base.:+(a::LazySum, b::LazySum)
    return LazySum(vcat(a.terms, b.terms))
end

function Base.:^(a::LazyTrace, n::Integer)
    if n == 0 return LazyTrace([], 1) end
    if n == 1 return a end
    new_cycles = Vector{Vector{AbstractMatrix}}()
    for _ = 1:n append!(new_cycles, a.cycles) end
    return LazyTrace(new_cycles, a.prefactor^n)
end

function show(io::IO, t::LazyTrace)
    if t.prefactor != 1 print(io, t.prefactor, "*") end
    if isempty(t.cycles) print(io, "1"); return end
    for (i, cycle) in enumerate(t.cycles)
        print(io, "tr(")
        for (j, f) in enumerate(cycle)
            print(io, f)
            if j < length(cycle) print(io, " * ") end
        end
        print(io, ")")
    end
end

"""
    tr_val(factors::AbstractVector)

Evaluates the trace of a product of matrices. If all factors are concrete, returns 
the numeric trace. If any factor is symbolic, returns a symbolic representation 
normalized by circular shifts and adjoints to ensure unique naming.
"""
function tr_val(factors::AbstractVector)
    if isempty(factors) return 1 end
    
    # Try to evaluate if all factors are concrete matrices
    # Check if any factor is symbolic
    is_symbolic = any(f -> f isa SymbolicMatrix || f isa SymbolicMatrixProduct || f isa LazyTrace, factors)
    
    if !is_symbolic
        try
            prod_val = prod(factors)
            if prod_val isa AbstractMatrix && eltype(prod_val) <: Number && !(eltype(prod_val) <: Num)
                  return LinearAlgebra.tr(prod_val)
            end
            if prod_val isa AbstractMatrix && eltype(prod_val) <: Num
                  all_num = true
                  for x in prod_val
                      if !is_number(x)
                          all_num = false; break
                      end
                  end
                  if all_num
                      return LinearAlgebra.tr(prod_val)
                  end
            end
        catch
        end
    end

    # Normalize trace by circular shifting and taking the smallest lexicographically.
    # We also consider the adjoint representation.
    function get_norm_string(fs)
        n = length(fs)
        s_list = [string(f) for f in fs]
        min_s = join(s_list, "*")
        for i = 1:n-1
            p = vcat(s_list[i+1:n], s_list[1:i])
            curr_s = join(p, "*")
            if curr_s < min_s
                min_s = curr_s
            end
        end
        return min_s
    end

    s1 = get_norm_string(factors)
    adj_factors = reverse([adjoint(f) for f in factors])
    s2 = get_norm_string(adj_factors)
    
    name = "tr(" * (s1 < s2 ? s1 : s2) * ")"
    # Use T=Number to preserve symbolic structure.
    return Num(Symbolics.variable(Symbol(name); T = Number))
end

function is_number(x)
    x = Symbolics.unwrap(x)
    return x isa Number
end
