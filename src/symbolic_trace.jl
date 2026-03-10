# Symbolics' Num is real-valued; we use AbstractMatrix{Any} to support complex symbols.
"""
    SymbolicMatrix(name::Symbol)
    SymbolicMatrix(name::Symbol, special_type::Symbol)
    SymbolicMatrix(name::Symbol, special_type::Symbol, dim)
    SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol, dim)

A wrapper associated with a symbolic name to represent a matrix in a coordinate-free way.
Used in the symbolic trace logic (via `tr_lazy`) and for metadata-driven element-wise integration.
"""
struct SymbolicMatrix <: AbstractMatrix{Any}
    name::Symbol
    is_adj::Bool
    is_trans::Bool
    special_type::Symbol
    dim::Union{Nothing,Integer,Num,Tuple{Union{Integer,Num},Union{Integer,Num}},Any}

    function SymbolicMatrix(
        name::Symbol,
        is_adj::Bool,
        is_trans::Bool,
        special_type::Symbol,
        dim::Any,
    )
        new(name, is_adj, is_trans, special_type, dim)
    end
end

"""
    IntegrationError(msg::String)

Custom exception thrown when a symbolic integration step fails (e.g., non-integer powers of traces).
"""
struct IntegrationError <: Exception
    msg::String
end

Base.showerror(io::IO, e::IntegrationError) = print(io, "IntegrationError: ", e.msg)

struct MatrixMetadata end

SymbolicMatrix(name::Symbol) = SymbolicMatrix(name, false, false, :Constant, nothing)
SymbolicMatrix(name::Symbol, special_type::Symbol) =
    SymbolicMatrix(name, false, false, special_type, nothing)
SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol) =
    SymbolicMatrix(name, is_adj, false, special_type, nothing)
SymbolicMatrix(name::Symbol, special_type::Symbol, dim) =
    SymbolicMatrix(name, false, false, special_type, dim)
SymbolicMatrix(name::Symbol, is_adj::Bool, special_type::Symbol, dim) =
    SymbolicMatrix(name, is_adj, false, special_type, dim)

import Base: *, adjoint, transpose, show, ^, size, getindex
import LinearAlgebra: tr

function Base.size(A::SymbolicMatrix)
    d = A.dim === nothing ? (nothing, nothing) : A.dim
    if d isa Tuple
        rows, cols = d
        return A.is_adj ? (cols, rows) : (rows, cols)
    else
        return (d, d)
    end
end

function Base.axes(A::SymbolicMatrix)
    sz = size(A)
    return map(s -> s isa Integer ? (Base.OneTo(s)) : (1:1), sz)
end

function _getindex_scalar(A::SymbolicMatrix, i, j)
    # Bounds checking
    if A.dim !== nothing
        rows, cols = size(A)
        # Check rows
        if i isa Integer
            if i < 1 || (rows isa Integer && i > rows)
                throw(BoundsError(A, (i, j)))
            end
        else
            # Symbolic check for i < 1 or i > rows
            # We use subtraction because simplify(d+1 > d) doesn't always return true
            si = Symbolics.simplify(i)
            if isequal(si, 0) ||
               isequal(si, -1) ||
               isequal(Symbolics.simplify(si < 1), true)
                throw(BoundsError(A, (i, j)))
            end

            # Check if i - rows is a positive constant
            sir = Symbolics.simplify(i - rows)
            if isequal(sir, 1) ||
               isequal(sir, 2) ||
               (sir isa Real && !(sir isa Num) && sir > 0)
                throw(BoundsError(A, (i, j)))
            end
        end
        # Check cols
        if j isa Integer
            if j < 1 || (cols isa Integer && j > cols)
                throw(BoundsError(A, (i, j)))
            end
        else
            # Symbolic check for j < 1 or j > cols
            sj = Symbolics.simplify(j)
            if isequal(sj, 0) ||
               isequal(sj, -1) ||
               isequal(Symbolics.simplify(sj < 1), true)
                throw(BoundsError(A, (i, j)))
            end

            sjr = Symbolics.simplify(j - cols)
            if isequal(sjr, 1) ||
               isequal(sjr, 2) ||
               (sjr isa Real && !(sjr isa Num) && sjr > 0)
                throw(BoundsError(A, (i, j)))
            end
        end
    end

    # If transposed, the row index of A^T is the col index of A
    actual_i = A.is_trans ? j : i
    actual_j = A.is_trans ? i : j

    s_name = Symbol(A.name, :_, actual_i, :_, actual_j)

    meta = Dict(
        :name => A.name,
        :type => A.special_type,
        :indices => (actual_i, actual_j),
        :is_adj => A.is_adj,
    )

    # Use T=Number to ensure Symbolics/SymbolicUtils does not incorrectly simplify conj(v)
    v = Symbolics.variable(s_name, T = Number)
    v_un = Symbolics.unwrap(v)
    v_meta_un = SymbolicUtils.setmetadata(v_un, MatrixMetadata, meta)
    v_meta = Symbolics.wrap(v_meta_un)

    return A.is_adj ? conj(v_meta) : v_meta
end

Base.getindex(A::SymbolicMatrix, i::Integer, j::Integer) = _getindex_scalar(A, i, j)
Base.getindex(A::SymbolicMatrix, i::Num, j::Num) = _getindex_scalar(A, i, j)
Base.getindex(A::SymbolicMatrix, i::Integer, j::Num) = _getindex_scalar(A, i, j)
Base.getindex(A::SymbolicMatrix, i::Num, j::Integer) = _getindex_scalar(A, i, j)

function Base.getindex(
    A::SymbolicMatrix,
    i::Union{Integer,AbstractVector,Colon},
    j::Union{Integer,AbstractVector,Colon},
)
    rows = (i isa Colon) ? (1:size(A, 1)) : i
    cols = (j isa Colon) ? (1:size(A, 2)) : j

    # Handle single element access by dispatching back to (Integer, Integer)
    if rows isa Integer && cols isa Integer
        return invoke(getindex, Tuple{SymbolicMatrix,Integer,Integer}, A, rows, cols)
    end

    res = Any[A[r, c] for r in rows, c in cols]
    # Return vector if single column/row was requested (standard Julia behavior)
    if length(cols) == 1 && (j isa Integer)
        return res[:, 1]
    elseif length(rows) == 1 && (i isa Integer)
        return res[1, :]
    end
    return SymbolicMatrixProduct([res])
end

function Base.adjoint(A::SymbolicMatrix)
    return conj(transpose(A))
end

function Base.transpose(A::SymbolicMatrix)
    # Pure transpose
    return SymbolicMatrix(A.name, A.is_adj, !A.is_trans, A.special_type, A.dim)
end

function Base.conj(A::SymbolicMatrix)
    # real special types: O, GOE, GinOE, Perm, CPerm
    if A.special_type in (:O, :GOE, :GinOE, :Perm, :CPerm)
        return A
    end
    # Pure conjugation
    return SymbolicMatrix(A.name, !A.is_adj, A.is_trans, A.special_type, A.dim)
end

# Factory functions for symbolic matrices
symbolic_unitary(name, d) = SymbolicMatrix(name, false, :U, d)
symbolic_orthogonal(name, d) = SymbolicMatrix(name, false, :O, d)
symbolic_symplectic(name, d) = SymbolicMatrix(name, false, :Sp, d)
symbolic_pure_state(name, d) = SymbolicMatrix(name, false, :psi, d)
symbolic_permutation(name, d) = SymbolicMatrix(name, false, :Perm, d)

function Base.show(io::IO, A::SymbolicMatrix)
    print(io, A.name)
    if A.is_adj && A.is_trans
        print(io, "'")
    elseif A.is_trans
        print(io, ".'")
    elseif A.is_adj
        print(io, "ᴴ") # Or conj... maybe just denote it somehow.
    end
end

function Base.show(io::IO, ::MIME"text/plain", A::SymbolicMatrix)
    print(io, A.name)
    if A.is_adj && A.is_trans
        print(io, "'")
    elseif A.is_trans
        print(io, ".'")
    elseif A.is_adj
        print(io, "ᴴ")
    end
    print(io, " (Symbolic Matrix")
    if A.dim !== nothing
        print(io, ", dim=", A.dim)
    else
        print(io, ", unspecified dimension")
    end
    if A.special_type !== :Constant
        print(io, ", type=", A.special_type)
    end
    print(io, ")")
end

"""
    LazyTrace(cycles::Vector{Vector{AbstractMatrix}}, prefactor::Union{Num, Number})

A lazy representation of a trace (or product of traces) of matrix products.
Used to represent expressions like `tr(A*B) * tr(C)` symbolically before integration.
"""
struct LazyTrace
    cycles::Vector{Vector{AbstractMatrix}}
    prefactor::Union{Num,Number}
end

function Base.getproperty(t::LazyTrace, s::Symbol)
    if s == :factors
        return vcat(t.cycles...)
    end
    return getfield(t, s)
end

"""
    LazySum(terms::Vector{LazyTrace})

A lazy representation of a sum of `LazyTrace` objects.
Enables symbolic integration of expressions like `tr(A*B) + tr(C*D)`.
"""
struct LazySum
    terms::Vector{LazyTrace}
end

"""
    LazyPower(base, exponent)

A lazy representation of a power of a `LazyTrace` or `LazySum`.
Used to represent expressions like `abs(tr(U))` as `(tr(U)*tr(U'))^0.5`.
"""
struct LazyPower
    base::Union{LazyTrace,LazySum}
    exponent::Any
end

struct SymbolicMatrixProduct <: AbstractMatrix{Any}
    factors::Vector{AbstractMatrix}
end

struct SymbolicKron <: AbstractMatrix{Any}
    A::AbstractMatrix
    B::AbstractMatrix
end

function Base.show(io::IO, P::SymbolicMatrixProduct)
    for (i, f) in enumerate(P.factors)
        show(io, f)
        if i < length(P.factors)
            print(io, " * ")
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", P::SymbolicMatrixProduct)
    show(io, P)
    print(io, " (Symbolic Matrix Product")
    sz = size(P)
    print(io, ", size=$(sz[1])x$(sz[2])")
    print(io, ", factors=$(length(P.factors))")
    print(io, ")")
end

const SymbolicAny = Union{SymbolicMatrix,SymbolicMatrixProduct,SymbolicKron}

function Base.size(P::SymbolicMatrixProduct)
    if isempty(P.factors)
        return (0, 0)
    end
    return (size(P.factors[1], 1), size(P.factors[end], 2))
end

function Base.size(P::SymbolicMatrixProduct, i::Integer)
    sz = size(P)
    return sz[i]
end

function Base.getindex(P::SymbolicMatrixProduct, i::Integer, j::Integer)
    factors = P.factors
    if length(factors) == 1
        return factors[1][i, j]
    end
    A = factors[1]
    B = length(factors) == 2 ? factors[2] : SymbolicMatrixProduct(factors[2:end])

    # Resolve inner dimension
    dimA = size(A, 2)
    dimB = size(B, 1)

    dim = nothing
    if dimA isa Integer && dimB isa Integer
        if dimA != dimB
            throw(
                DimensionMismatch(
                    "matrix A has dimensions $(size(A)), matrix B has dimensions $(size(B))",
                ),
            )
        end
        dim = dimA
    elseif dimA isa Integer
        dim = dimA
    elseif dimB isa Integer
        dim = dimB
    end

    if dim isa Integer
        return Symbolics.wrap(sum(A[i, k] * B[k, j] for k = 1:dim))
    else
        # Fallback to a symbolic sum variable if dimension is unknown
        return Num(Symbolics.variable(Symbol("sum_$(A)_$(B)_$(i)_$(j)"); T = Number))
    end
end

function Base.size(K::SymbolicKron)
    szA = size(K.A)
    szB = size(K.B)

    # helper to multiply sizes part by part
    mul_dim(a, b) = (a === nothing || b === nothing) ? nothing : a * b

    return (mul_dim(szA[1], szB[1]), mul_dim(szA[2], szB[2]))
end

function Base.getindex(K::SymbolicKron, i::Integer, j::Integer)
    szB = size(K.B)
    rowsB = szB[1]
    colsB = szB[2]

    if rowsB === nothing || colsB === nothing
        # Throw a helpful error instead of MethodError in divrem
        throw(
            ArgumentError(
                "Cannot index into SymbolicKron with unknown dimensions in factor $(K.B). Specify dimensions or use tr() for scalar results.",
            ),
        )
    end

    iA, iB = divrem(i - 1, rowsB) .+ 1
    jA, jB = divrem(j - 1, colsB) .+ 1
    return K.A[iA, jA] * K.B[iB, jB]
end

function Base.adjoint(K::SymbolicKron)
    return SymbolicKron(adjoint(K.A), adjoint(K.B))
end

function Base.transpose(K::SymbolicKron)
    return SymbolicKron(transpose(K.A), transpose(K.B))
end

function Base.show(io::IO, K::SymbolicKron)
    print(io, "(")
    show(io, K.A)
    print(io, " ⊗ ")
    show(io, K.B)
    print(io, ")")
end

# Multiplication logic
_factors(A::SymbolicMatrix) = AbstractMatrix[A]
_factors(P::SymbolicMatrixProduct) = P.factors
_factors(A::SymbolicKron) = AbstractMatrix[A]
_factors(A::AbstractMatrix) = AbstractMatrix[A]

function *(A::SymbolicAny, B::SymbolicAny)
    if A isa SymbolicKron && B isa SymbolicKron
        return SymbolicKron(A.A * B.A, A.B * B.B)
    end
    return SymbolicMatrixProduct(vcat(_factors(A), _factors(B)))
end
function *(A::SymbolicAny, B::AbstractMatrix)
    return SymbolicMatrixProduct(vcat(_factors(A), AbstractMatrix[B]))
end
function *(A::AbstractMatrix, B::SymbolicAny)
    return SymbolicMatrixProduct(vcat(AbstractMatrix[A], _factors(B)))
end

function *(A::SymbolicAny, B::SymbolicAny, Cs::SymbolicAny...)
    res_factors = vcat(_factors(A), _factors(B))
    for C in Cs
        append!(res_factors, _factors(C))
    end
    # Check if we can collapse all as kron
    if all(f -> f isa SymbolicKron, res_factors)
        return reduce((a, b) -> SymbolicKron(a.A * b.A, a.B * b.B), res_factors)
    end
    return SymbolicMatrixProduct(res_factors)
end

# kron overloads
function LinearAlgebra.kron(A::SymbolicAny, B::SymbolicAny)
    return SymbolicKron(A, B)
end
function LinearAlgebra.kron(A::SymbolicAny, B::AbstractMatrix)
    return SymbolicKron(A, B)
end
function LinearAlgebra.kron(A::AbstractMatrix, B::SymbolicAny)
    return SymbolicKron(A, B)
end

# Resolve ambiguities with LinearAlgebra and Symbolics
using LinearAlgebra
for T in [Adjoint, Transpose]
    @eval function *(A::$T{<:Any,<:AbstractVector}, B::SymbolicAny)
        return SymbolicMatrixProduct(vcat(Any[A], _factors(B)))
    end
    @eval function *(A::SymbolicAny, B::$T{<:Any,<:AbstractVector})
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
    @eval function *(A::$T{<:Any,<:AbstractVector}, B::SymbolicAny, C::SymbolicAny)
        return (A * B) * C
    end
    @eval function *(A::$T{<:Any,<:AbstractVector}, B::AbstractMatrix, C::SymbolicAny)
        return (A * B) * C
    end
    @eval function *(A::$T{<:Any,<:AbstractVector}, B::SymbolicAny, C::AbstractMatrix)
        return (A * B) * C
    end

    # 4-arg (to prevent LinearAlgebra from catching it)
    @eval function *(
        A::$T{<:Any,<:AbstractVector},
        B::AbstractMatrix,
        C::AbstractMatrix,
        D::SymbolicAny,
    )
        return (A * B) * (C * D)
    end
    @eval function *(
        A::$T{<:Any,<:AbstractVector},
        B::AbstractMatrix,
        C::SymbolicAny,
        D::AbstractMatrix,
    )
        return (A * B) * (C * D)
    end

    # Disambiguate with internal overlaps
    @eval function *(
        A::$T{<:Any,<:AbstractVector},
        B::AbstractMatrix,
        C::SymbolicAny,
        D::SymbolicAny,
    )
        return (A * B) * (C * D)
    end

    # Disambiguate with Symbolics.Arr (3-arg)
    @eval function *(A::$T{W,<:AbstractVector}, B::Symbolics.Arr, C::SymbolicAny) where {W}
        return (A * B) * C
    end
    @eval function *(
        A::$T{W,<:AbstractVector},
        B::Symbolics.Arr,
        C::SymbolicAny,
    ) where {W<:Number}
        return (A * B) * C
    end

    # Disambiguate with Symbolics.Arr (4-arg)
    for W_type in [Any, Number]
        @eval function *(
            A::$T{W,<:AbstractVector},
            B::Symbolics.Arr,
            C::AbstractMatrix,
            D::SymbolicAny,
        ) where {W<:$W_type}
            return (A * B) * (C * D)
        end
        @eval function *(
            A::$T{W,<:AbstractVector},
            B::Symbolics.Arr,
            C::SymbolicAny,
            D::AbstractMatrix,
        ) where {W<:$W_type}
            return (A * B) * (C * D)
        end
        @eval function *(
            A::$T{W,<:AbstractVector},
            B::Symbolics.Arr,
            C::SymbolicAny,
            D::SymbolicAny,
        ) where {W<:$W_type}
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

function _is_identity(A)
    if A isa AbstractMatrix && !(A isa IntU.SymbolicAny)
        return A == I || (size(A, 1) == size(A, 2) && A == I(size(A, 1)))
    end
    return false
end

function _are_inverses(A, B)
    if A isa IntU.SymbolicMatrix && B isa IntU.SymbolicMatrix
        if A.name === B.name &&
           A.special_type === B.special_type &&
           isequal(A.dim, B.dim) &&
           A.is_adj != B.is_adj
            if A.special_type in
               (:U, :O, :Sp, :CUE, :COE, :CSE, :Perm, :CPerm, :DiagUnitary)
                return true
            end
        end
    end

    if A isa IntU.SymbolicKron && B isa IntU.SymbolicKron
        cancel1 = _are_inverses(A.A, B.A) || (_is_identity(A.A) && _is_identity(B.A))
        cancel2 = _are_inverses(A.B, B.B) || (_is_identity(A.B) && _is_identity(B.B))
        return cancel1 && cancel2
    end

    return false
end

function _simplify_cycle(factors::AbstractVector)
    if isempty(factors)
        return factors
    end

    changed = true
    current_factors = copy(factors)

    while changed && length(current_factors) >= 2
        changed = false
        n = length(current_factors)

        for i = 1:n
            j = (i % n) + 1
            f1 = current_factors[i]
            f2 = current_factors[j]

            cancels = _are_inverses(f1, f2) || _are_inverses(f2, f1)

            if cancels
                if i < j
                    deleteat!(current_factors, j)
                    deleteat!(current_factors, i)
                else
                    deleteat!(current_factors, i)
                    deleteat!(current_factors, j)
                end
                changed = true
                break
            end
        end
    end
    return current_factors
end

"""
    tr(A::SymbolicMatrix)
    tr(A::SymbolicMatrixProduct)

Symbolic trace of a coordinate-free matrix expression. 
Returns a `LazyTrace` object that can be integrated.
"""
function tr(A::SymbolicMatrix)
    return tr_lazy(A)
end
function tr(A::SymbolicMatrixProduct)
    # Check if any factor is a non-symbolic matrix (e.g. Matrix{Num} from kron)
    # Matrix{Num} is usually the result of kron(U, U) or similar.
    # We must expand such products to allow element-wise integration to find the unitaries.
    is_dirty = any(f -> !(f isa SymbolicAny) || f isa SymbolicKron, A.factors)

    if is_dirty
        rows, cols = size(A)
        # We can only expand if dimensions are concrete integers
        if rows isa Integer && cols isa Integer && rows == cols
            # Return a Num expression (sum of diagonals)
            return sum(i -> A[i, i], 1:rows)
        end
        # Fallback to lazy if we can't expand, though it might fail integration
    end

    return tr_lazy(A.factors)
end

function tr(K::SymbolicKron)
    return tr(K.A) * tr(K.B)
end

"""
    tr_lazy(product)

Creates a `LazyTrace` representing the symbolic trace of a matrix product.
The product can be a `SymbolicMatrix`, `SymbolicMatrixProduct`, or a vector of matrices.
"""
function tr_lazy(product::AbstractVector)
    simplified_factors = _simplify_cycle(product)
    return LazyTrace([collect(Any, simplified_factors)], 1)
end
function tr_lazy(product::SymbolicMatrix)
    return LazyTrace([[product]], 1)
end
function tr_lazy(product::SymbolicMatrixProduct)
    simplified_factors = _simplify_cycle(product.factors)
    return LazyTrace([simplified_factors], 1)
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
    new_cycles = [reverse([adjoint(f) for f in cycle]) for cycle in t.cycles]
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
    if n == 0
        return LazyTrace([], 1)
    end
    if n == 1
        return a
    end
    return LazyTrace(repeat(a.cycles, n), a.prefactor^n)
end

function Base.:^(a::Union{LazyTrace,LazySum}, n::Any)
    return LazyPower(a, n)
end

function Base.:^(a::LazyPower, n::Any)
    return LazyPower(a.base, a.exponent * n)
end

function Base.abs(t::Union{LazyTrace,LazySum})
    return (t * conj(t))^0.5
end

function Base.sqrt(t::Union{LazyTrace,LazySum})
    return t^0.5
end

function show(io::IO, t::LazyTrace)
    if t.prefactor != 1
        print(io, t.prefactor, "*")
    end
    if isempty(t.cycles)
        print(io, "1");
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

function show(io::IO, lp::LazyPower)
    print(io, "(")
    show(io, lp.base)
    print(io, ")^", lp.exponent)
end

"""
    tr_val(factors::AbstractVector)

Evaluates the trace of a product of matrices. If all factors are concrete, returns 
the numeric trace. If any factor is symbolic, returns a symbolic representation 
normalized by circular shifts and adjoints to ensure unique naming.
"""
function tr_val(factors::AbstractVector)
    if isempty(factors)
        return 1
    end

    # Try to evaluate if all factors are concrete matrices
    # Check if any factor is symbolic
    is_symbolic = any(
        f ->
            f isa SymbolicMatrix ||
            f isa SymbolicMatrixProduct ||
            f isa LazyTrace ||
            f isa SymbolicKron,
        factors,
    )

    if !is_symbolic
        prod_val = prod(factors)
        if prod_val isa AbstractMatrix &&
           eltype(prod_val) <: Number &&
           !(eltype(prod_val) <: Num)
            return LinearAlgebra.tr(prod_val)
        end
        if prod_val isa AbstractMatrix && eltype(prod_val) <: Num
            all_num = true
            for x in prod_val
                if !is_number(x)
                    all_num = false;
                    break
                end
            end
            if all_num
                return LinearAlgebra.tr(prod_val)
            end
        end
    end

    # Normalize trace string for unique symbolic representation
    function get_norm_string(fs)
        n = length(fs)
        s_list = [string(f) for f in fs]
        min_s = join(s_list, "*")
        for i = 1:(n-1)
            p = vcat(s_list[(i+1):n], s_list[1:i])
            curr_s = join(p, "*")
            if curr_s < min_s
                min_s = curr_s
            end
        end
        return min_s
    end

    s1 = get_norm_string(factors)
    trans_factors = reverse([transpose(f) for f in factors])
    s2 = get_norm_string(trans_factors)

    name = "tr(" * (s1 < s2 ? s1 : s2) * ")"
    # Use T=Number to preserve symbolic structure.
    return Num(Symbolics.variable(Symbol(name); T = Number))
end

function is_number(x)
    x = Symbolics.unwrap(x)
    if x isa Number
        return true
    end
    return hasproperty(x, :val) ? x.val isa Number : false
end
