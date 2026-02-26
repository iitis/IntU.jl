# Stiefel Manifold integration

@doc raw"""
    dStiefel(dim, k)
    dStiefel(V::SymbolicMatrix, k)

Defines the measure for integration over the Stiefel manifold $V_k(\mathbb{C}^d)$.
This manifold represents the set of $d \times k$ matrices with orthonormal columns.

The integration is performed by mapping $V$ to the first $k$ columns of a Haar-random unitary matrix $U(d)$.
If called with `dim`, it integrates entries tagged with `:U` via `SymbolicMatrix(:V, :U)`.

Reference:
- Edelman, A., Arias, T. A., & Smith, S. T. (1998). The geometry of algorithms with orthogonality constraints.
"""
dStiefel(dim, k) = StiefelMeasure(dim, k)


"""
    StiefelMeasure(dim, k)

Internal type representing the measure on the Stiefel manifold. 
Users should use `dStiefel` constructors.
"""
struct StiefelMeasure{D,K,M} <: AbstractMeasure
    dim::D
    k::K
    matcher::M
end
StiefelMeasure(dim, k) = StiefelMeasure(dim, k, nothing)

function IntU.measure_info(measure::StiefelMeasure)
    subs_dict = Dict{Any,Any}()
    # Default to matching :U, as per docs and examples
    matcher = measure.matcher === nothing ? MetadataMatcher(:U) : measure.matcher
    dim = measure.dim
    if dim isa SymbolicMatrix
        dim = dim.dim
    end
    # The integration rule uses :U logic (mapping V -> U P, where P is projection)
    return (subs_dict, matcher, dim, :U)
end

function integrate(P::SymbolicMatrixProduct, measure::StiefelMeasure)
    if isempty(P.factors)
        return Num(1)
    end

    dim_d = measure.dim
    dim_k = measure.k

    rows_sym = size(P.factors[1], 1)
    cols_sym = size(P.factors[end], 2)

    # Helper to resolve dimension symbols to measure dimensions
    function resolve_dim(d_sym)
        d_un = Symbolics.unwrap(d_sym)
        if d_un isa Integer && d_un != typemax(Int)
            return d_un
        end
        return dim_d # preferring measure dim d for expansion if matrix dim is unknown/symbolic
    end


    if dim_d isa Integer && dim_k isa Integer

        # We assume standard SymbolicMatrix/Adjoint wrappers
        function get_size(F)
            # Unwrap Adjoint/Transpose
            wrapper_adj = false
            if F isa Adjoint || F isa Transpose
                wrapper_adj = true
                F = parent(F)
            end

            if F isa SymbolicMatrix

                # Combine wrapper adjacency with internal adjacency
                effective_adj = F.is_adj
                if wrapper_adj
                    effective_adj = !effective_adj
                end

                # We treat it as d x k.
                if effective_adj
                    return (dim_k, dim_d)
                else
                    return (dim_d, dim_k)
                end
            end

            # Fallback for known matrices
            sz = size(F)
            if wrapper_adj
                return (sz[2], sz[1])
            end
            return sz
        end

        # 1. Expand all factors into dense Matrix{Any} of symbols
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

        # 2. Multiply them using standard linear algebra
        # This expands the sums symbolically
        total_prod = numeric_factors[1]
        for idx = 2:length(numeric_factors)
            total_prod = total_prod * numeric_factors[idx]
        end

        # 3. Integrate the result element-wise
        return map(x -> integrate(x, measure), total_prod)
    end

    # Fallback to standard symbolic integration if dimensions are not concrete integers
    throw(ArgumentError("Direct integration of SymbolicMatrixProduct for Stiefel requires concrete dimensions. Try integrating individual elements instead."))
end

"""
    asymptotic(expr, measure::StiefelMeasure, order=1)

Returns the series expansion of the integral in powers of `1/d`.
"""
function asymptotic(expr, measure::StiefelMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dStiefel(d_asymp, measure.k)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end
