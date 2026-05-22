
struct OrthogonalMeasure{D} <: AbstractMeasure
    dim::D
end

struct SymplecticMeasure{D} <: AbstractMeasure
    dim::D
    function SymplecticMeasure(dim::D) where {D}
        d_int = _try_extract_int(dim)
        if d_int !== nothing && isodd(d_int)
            throw(
                ArgumentError(
                    "Dimension dim must be even for SymplecticMeasure, got $dim.",
                ),
            )
        end
        new{D}(dim)
    end
end

@doc raw"""
    dO(dim)

Defines the Haar measure for the real Orthogonal group $O(d)$ with dimension `dim`.
Integration engine identifies variables via metadata tag `:O`.
"""
function dO(dim)
    _assert_no_float_param(dim, "dim", "dO")
    return OrthogonalMeasure(dim)
end

@doc raw"""
    dSp(dim)

Defines the Haar measure for the Symplectic group $Sp(d)$. 
The dimension `dim` must be even.
Integration engine identifies variables via metadata tag `:Sp`.
"""
function dSp(dim)
    _assert_no_float_param(dim, "dim", "dSp")
    return SymplecticMeasure(dim)
end


IntegrateUnitary._measure_tag(::OrthogonalMeasure) = :O
IntegrateUnitary._measure_tag(::SymplecticMeasure) = :Sp

function _j_pair_sign(idx, n)
    if idx <= n
        return idx + n, 1
    else
        return idx - n, -1
    end
end

IntegrateUnitary._reconstruct_symbolic(::OrthogonalMeasure, d_asymp) = dO(d_asymp)
IntegrateUnitary._reconstruct_symbolic(::SymplecticMeasure, d_asymp) = dSp(d_asymp)
