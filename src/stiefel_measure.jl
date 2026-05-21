
@doc raw"""
    dStiefel(dim, k)

Defines the measure for integration over the Stiefel manifold $V_k(\mathbb{C}^d)$.
This manifold represents the set of $d \times k$ matrices with orthonormal columns.

The integration is performed by mapping $V$ to the first $k$ columns of a Haar-random unitary matrix $U(d)$.
If called with `dim`, it integrates entries tagged with `:V` via `SymbolicMatrix(:V, :V)`.

Reference:
- Edelman, A., Arias, T. A., & Smith, S. T. (1998). The geometry of algorithms with orthogonality constraints.
"""
function dStiefel(dim, k)
    _assert_no_float_param(dim, "dim", "dStiefel")
    _assert_no_float_param(k, "k", "dStiefel")
    return StiefelMeasure(dim, k)
end


"""
    StiefelMeasure(dim, k)

Internal type representing the measure on the Stiefel manifold. 
Users should use `dStiefel` constructors.
"""
struct StiefelMeasure{D,K,M} <: AbstractMeasure
    dim::D
    k::K
    matcher::M
    function StiefelMeasure(dim::D, k::K, matcher::M) where {D,K,M}
        d_int = _try_extract_int(dim)
        k_int = _try_extract_int(k)
        if d_int !== nothing && k_int !== nothing && k_int > d_int
            throw(
                ArgumentError(
                    "Stiefel manifold V_k(C^d) requires k <= d, got k=$k, d=$dim",
                ),
            )
        end
        new{D,K,M}(dim, k, matcher)
    end
end
StiefelMeasure(dim, k) = StiefelMeasure(dim, k, nothing)

IntU._measure_tag(::StiefelMeasure) = :V

function IntU.measure_info(measure::StiefelMeasure)
    subs_dict = Dict{Any,Any}()
    tag = IntU._measure_tag(measure)
    matcher = measure.matcher === nothing ? MetadataMatcher(tag) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    dim = _assert_no_float_param(dim, "dim", "StiefelMeasure")
    k = _assert_no_float_param(measure.k, "k", "StiefelMeasure")
    return (subs_dict, matcher, dim, (tag, k))
end

function integrate(P::SymbolicMatrixProduct, measure::StiefelMeasure)
    _validate_measure_discrete_params(measure)

    if isempty(P.factors)
        return Num(1)
    end

    dim_d = measure.dim
    dim_k = measure.k

    rows_sym = size(P.factors[1], 1)
    cols_sym = size(P.factors[end], 2)

    function resolve_dim(d_sym)
        d_un = Symbolics.unwrap(d_sym)
        if d_un isa Integer && d_un != typemax(Int)
            return d_un
        end
        return dim_d # preferring measure dim d for expansion if matrix dim is unknown/symbolic
    end


    if dim_d isa Integer && dim_k isa Integer

        function get_size(F)
            wrapper_adj = false
            if F isa Adjoint || F isa Transpose
                wrapper_adj = true
                F = parent(F)
            end

            if F isa SymbolicMatrix

                effective_adj = F.is_adj
                if wrapper_adj
                    effective_adj = !effective_adj
                end

                if effective_adj
                    return (dim_k, dim_d)
                else
                    return (dim_d, dim_k)
                end
            end

            sz = size(F)
            if wrapper_adj
                return (sz[2], sz[1])
            end
            return sz
        end

        numeric_factors = Vector{Matrix{Any}}(undef, length(P.factors))

        for (idx, f) in enumerate(P.factors)
            rows, cols = get_size(f)
            mat = Matrix{Any}(undef, rows, cols)
            for i = 1:rows
                for j = 1:cols
                    mat[i, j] = f[i, j]
                end
            end
            numeric_factors[idx] = mat
        end

        total_prod = numeric_factors[1]
        for idx = 2:length(numeric_factors)
            total_prod = total_prod * numeric_factors[idx]
        end

        return map(x -> integrate(x, measure), total_prod)
    end

    throw(
        ArgumentError(
            "Direct integration of SymbolicMatrixProduct for Stiefel requires concrete dimensions. Try integrating individual elements instead.",
        ),
    )
end

IntU._reconstruct_symbolic(m::StiefelMeasure, d_asymp) = dStiefel(d_asymp, m.k)
