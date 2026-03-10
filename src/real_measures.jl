# Real and Symplectic measures

# Dummy types to represent the measures
struct OrthogonalMeasure{D} <: AbstractMeasure
    dim::D
end

struct SymplecticMeasure{D} <: AbstractMeasure
    dim::D
    function SymplecticMeasure(dim::D) where {D}
        if dim isa Integer && isodd(dim)
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
dO(dim) = OrthogonalMeasure(dim)

@doc raw"""
    dSp(dim)

Defines the Haar measure for the Symplectic group $Sp(d)$. 
The dimension `dim` must be even.
Integration engine identifies variables via metadata tag `:Sp`.
"""
dSp(dim) = SymplecticMeasure(dim)


IntU._measure_tag(::OrthogonalMeasure) = :O
IntU._measure_tag(::SymplecticMeasure) = :Sp

function _j_pair_sign(idx, n)
    if idx <= n
        return idx + n, 1
    else
        return idx - n, -1
    end
end

IntU._reconstruct_symbolic(::OrthogonalMeasure, d_asymp) = dO(d_asymp)
IntU._reconstruct_symbolic(::SymplecticMeasure, d_asymp) = dSp(d_asymp)
