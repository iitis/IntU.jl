# Real and Symplectic measures

# Dummy types to represent the measures
struct OrthogonalMeasure{T, N, D}
    O::AbstractArray{T, N}
    dim::D
end

struct SymplecticMeasure{T, N, D}
    S::AbstractArray{T, N}
    dim::D
end

"""
    dO(O, dim)

Defines the Haar measure for the real Orthogonal group O(d).

The integration of a monomial of entries is given by:
```math
\\int_{O(d)} O_{i_1 j_1} \\dots O_{i_{2n} j_{2n}} dO = \\sum_{\\pi, \\sigma \\in M_{2n}} \\delta_{\\pi}(i) \\delta_{\\sigma}(j) \\text{Wg}^O(\\pi, \\sigma, d)
```
where M_{2n} is the set of pair partitions.

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
dO(O::AbstractArray{T,N}, dim) where {T,N} = OrthogonalMeasure{T,N,typeof(dim)}(O, dim)

"""
    dSp(S, dim)

Defines the Haar measure for the Symplectic group Sp(d). 
The dimension `dim` must be even.

The integration formula uses the symplectic metric J and pair partitions:
```math
\\int_{Sp(d)} S_{i_1 j_1} \\dots S_{i_{2n} j_{2n}} dS = \\sum_{\\pi, \\sigma \\in M_{2n}} J_{\\pi}(i) J_{\\sigma}(j) \\text{Wg}^{Sp}(\\pi, \\sigma, d)
```

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
dSp(S::AbstractArray{T,N}, dim) where {T,N} = SymplecticMeasure{T,N,typeof(dim)}(S, dim)


"""
    integrate(expr, measure::OrthogonalMeasure)
"""
function integrate(expr::AbstractArray, measure::OrthogonalMeasure)
    return map(e -> integrate(e, measure), expr)
end

function integrate(expr::AbstractArray, measure::SymplecticMeasure)
    return map(e -> integrate(e, measure), expr)
end

function fallback_integrate(expr, measure::OrthogonalMeasure)
    O_sym = measure.O
    dim = measure.dim
    
    subs_dict = Dict{Any, Any}()
    O_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    
    if O_sym isa AbstractArray
        for i in 1:size(O_sym, 1)
            for j in 1:size(O_sym, 2)
                o_ij_num = _safe_Num(O_sym[i,j])
                o_ij_un = Symbolics.unwrap(o_ij_num)
                o_atomic = Symbolics.variable(:O_atomic, i, j)
                
                O_atomic_lookup[Symbolics.unwrap(o_atomic)] = (i, j)
                
                # O is real, so conj(O) = O
                subs_dict[o_ij_un] = o_atomic
                
                # Handle conjugates by mapping them to the same atomic variable
                subs_dict[Symbolics.unwrap(conj(o_ij_un))] = o_atomic
                subs_dict[Symbolics.unwrap(Base.conj(o_ij_un))] = o_atomic
            end
        end
    end

    # Pass empty U_bar_lookup because O is orthogonal (real for our integration purposes effectively)
    # We will reuse the core logic but identifying it as orthogonal measure
    matcher = LookupMatcher(O_atomic_lookup, Dict{Any, Tuple{Int, Int}}())
    return _robust_real_num(_integrate_core(expr, dim, subs_dict, matcher, :O))
end

function fallback_integrate(expr, measure::SymplecticMeasure)
    S_sym = measure.S
    dim = measure.dim
    
    subs_dict = Dict{Any, Any}()
    S_atomic_lookup = Dict{Any, Tuple{Int, Int}}()
    
    # Check dimensions
    N = size(S_sym, 1)
    if isodd(N)
        error("Symplectic matrix dimension must be even, got \$N")
    end
    n_half = N ÷ 2
    
    # Pre-create all atomic variables first
    atomics = Matrix{Any}(undef, N, N)
    
    if S_sym isa AbstractArray
        for i in 1:N
            for j in 1:N
                s_atomic = Symbolics.variable(:S_atomic, i, j)
                atomics[i, j] = s_atomic
                
                # Store lookup
                S_atomic_lookup[Symbolics.unwrap(s_atomic)] = (i, j)
            end
        end
        
        # Now define substitutions
        for i in 1:N
            for j in 1:N
                s_ij_num = _safe_Num(S_sym[i,j])
                s_ij_un = Symbolics.unwrap(s_ij_num)
                
                s_atomic = atomics[i, j]
                subs_dict[s_ij_un] = s_atomic
                
                # Handle conj(S_ij) -> S_pq * coeff
                # \bar{S}_{ij} = J_{ip} S_{pq} (J^T)_{qj} = J_{ip} S_{pq} (-J_{qj})
                # p = pair(i), q = pair(j)
                
                # Helper for J-pair and sign
                function j_pair_sign(idx, n)
                    if idx <= n
                        return idx + n, 1
                    else
                        return idx - n, -1
                    end
                end
                
                p, sign_i = j_pair_sign(i, n_half)
                q, sign_j = j_pair_sign(j, n_half) # Note: this calculates sign for J_{j, q} NOT J_{q, j}
                # We need -J_{q, j}.
                # J_{pair(j), j}.
                # if j <= n, q=j+n. J_{j+n, j} = -1. -J = 1.
                # if j > n, q=j-n. J_{j-n, j} = 1. -J = -1.
                # My logic in thought process: "sign(l, n)" was "if l <= n then 1 else -1".
                # Let's re-verify.
                # Term is -J_{q, j}. q = pair(j).
                # If j <= n, q = j+n. J_{q, j} = J_{j+n, j} = -1. Term = 1.
                # If j > n, q = j-n. J_{q, j} = J_{j-n, j} = 1. Term = -1.
                # So coeff_j = (j <= n ? 1 : -1).
                # My j_pair_sign returns 1 if idx<=n. Correct.
                
                coeff = sign_i * sign_j
                
                s_mapped = atomics[p, q]
                
                subs_dict[Symbolics.unwrap(conj(s_ij_un))] = coeff * s_mapped
                subs_dict[Symbolics.unwrap(Base.conj(s_ij_un))] = coeff * s_mapped
            end
        end
    end

    matcher = LookupMatcher(S_atomic_lookup, Dict{Any, Tuple{Int, Int}}())
    return _robust_real_num(_integrate_core(expr, dim, subs_dict, matcher, :Sp))
end

"""
    asymptotic(expr, measure::Union{OrthogonalMeasure, SymplecticMeasure}, order=1)
"""
function asymptotic(expr, measure::Union{OrthogonalMeasure, SymplecticMeasure}, order=1)
    d = measure.dim
    # If d is symbolic or not an integer, we can proceed directly
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end
    
    # If d is integer, checking asymptotic might require symbolic d
    d_asymp = Symbolics.variable(:d_asymp)
    # Reconstruct measure with symbolic dim
    m_sym = if measure isa OrthogonalMeasure
        dO(measure.O, d_asymp)
    else
        dSp(measure.S, d_asymp)
    end
    
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end
