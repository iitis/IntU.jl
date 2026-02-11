# Haar measure integration

"""
    LazySymbolicMatrix{D, T}

A lazy, "infinite" representation of a matrix with entries of type `T`.
It does not store entries in memory. Instead, accessing `M[i,j]` generates
a symbolic variable on the fly using `func`. This allows for operations
on matrices of symbolic or very large dimensions without memory overhead.
"""
struct LazySymbolicMatrix{D,T} <: AbstractMatrix{T}
    name::Symbol
    func::Function # Generates the variable
    dim::D
end

# Infinite size for lazy matrix to allow M[100, 100] etc
Base.size(::LazySymbolicMatrix) = (typemax(Int), typemax(Int))
Base.getindex(M::LazySymbolicMatrix, i::Integer, j::Integer) = M.func(i, j)

# Custom show to avoid printing infinite matrix
function Base.show(io::IO, M::LazySymbolicMatrix{D,T}) where {D,T}
    print(io, "LazySymbolicMatrix{$(T)}($(M.name))")
end
function Base.show(io::IO, ::MIME"text/plain", M::LazySymbolicMatrix{D,T}) where {D,T}
    print(io, "LazySymbolicMatrix{$(T)}($(M.name), dim=$(M.dim))")
end

"""
    symbolic_dimension_matrix(dim; name=:M, T=Complex{Num})

Creates a lazy symbolic matrix `M` with symbolic dimension `dim` and entry type `T`.
`M` allows accessing indices of any size (e.g. `M[i,j]`).
"""
function symbolic_dimension_matrix(dim; name = :M, T = Complex{Num})
    # Create variables on demand with consistent naming
    gen = (i, j) -> Symbolics.variable(Symbol(name, :_, i, :_, j), T = T)
    return LazySymbolicMatrix{typeof(dim),T}(name, gen, dim)
end

"""
    @symbolic_dimension M[1:d, 1:d]
    @symbolic_dimension M[1:d, 1:d] T=Real

Macro to define a `LazySymbolicMatrix` matrix `M` with symbolic dimension `d`.
Optionally specify the element type `T` (e.g. `Real` or `Complex`). Defaults to `Complex{Num}`.

Example:
```julia
@variables d
@symbolic_dimension U[1:d, 1:d]
measure = dU(U)

@symbolic_dimension O[1:d, 1:d] T=Real
measure = dO(O, d)
```
"""
macro symbolic_dimension(expr, options...)
    if !Meta.isexpr(expr, :ref) || length(expr.args) != 3
        error("Usage: @symbolic_dimension M[1:d, 1:d] [T=Type]")
    end

    name = expr.args[1]

    # Parse dimensions from ranges
    # Expecting 1:d
    range1 = expr.args[2]

    # Simple check for 1:d format
    if !Meta.isexpr(range1, :call) || range1.args[1] != :(:) || range1.args[2] != 1
        error("Dimensions must be in format 1:d")
    end

    dim_sym = range1.args[3]
    
    # Parse options
    T_val = :(Complex{Num})
    for opt in options
        if Meta.isexpr(opt, :(=)) && opt.args[1] == :T
             # Check if user passed Real, convert to valid symbolic type if needed
             # or just pass it through. Symbolics usually wants `Real` or `Complex{Num}`.
             # If user passes `Real`, `Symbolics.variable(..., T=Real)` works.
             T_val = opt.args[2]
        end
    end

    q = quote
        $(esc(name)) =
            symbolic_dimension_matrix($(esc(dim_sym)), name = $(QuoteNode(name)), T = $(esc(T_val)))
    end
    return q
end

# Matcher for LazySymbolicMatrix
struct SymbolicMatcher <: AbstractIndexMatcher
    regex::Regex
end

function match_index(m::SymbolicMatcher, t)
    # Unwrap to Sym
    s = Symbolics.unwrap(t)
    # We check string representation of the symbol
    s_str = string(s)

    mat = match(m.regex, s_str)
    if mat !== nothing
        try
            i = parse(Int, mat[1])
            j = parse(Int, mat[2])
            return (:U, i, j)
        catch
        end
    end
    return nothing
end


# Dummy type to represent the measure
struct HaarMeasure{M,D}
    U::M
    dim::D
end
"""
    dU(U, dim)
    dU(U::LazySymbolicMatrix)

Defines the Haar measure for the Unitary group U(d).

The integration of a monomial of entries is given by:
```math
\\int_{U(d)} U_{i_1 j_1} \\dots U_{i_n j_n} \\bar{U}_{k_1 l_1} \\dots \\bar{U}_{k_n l_n} dU = \\sum_{\\sigma, \\tau \\in S_n} \\delta_{i, k_\\sigma} \\delta_{j, l_\\tau} \\text{Wg}(\\sigma \\tau^{-1}, d)
```

Reference:
- Collins, B., & Śniady, P. (2006). Integration with respect to the Haar measure on unitary, orthogonal and symplectic groups.
"""
dU(U, dim) = HaarMeasure(U, dim)
dU(dim) = HaarMeasure(nothing, dim)
dU(S::LazySymbolicMatrix) = HaarMeasure(S, S.dim)


"""
    measure_info(measure)

Returns a tuple `(matcher, dim, type)` for the given measure. 
Internal function used for dispatching integration logic.
"""
function IntU.measure_info(measure::HaarMeasure)
    U_input = measure.U
    dim = measure.dim

    # Check if U is our new LazySymbolicMatrix or legacy AbstractArray
    if U_input isa LazySymbolicMatrix
        subs_dict = Dict{Any,Any}()
        matcher = SymbolicMatcher(Regex("^$(U_input.name)_(\\d+)_(\\d+)\$"))

        return (subs_dict, matcher, dim, :U)
    else
        U_sym = U_input
        subs_dict = Dict{Any,Any}()
        U_atomic_lookup = Dict{Any,Tuple}()
        U_bar_lookup = Dict{Any,Tuple}()

        if U_sym isa AbstractArray
            for i = 1:size(U_sym, 1)
                for j = 1:size(U_sym, 2)
                    u_ij_num = _safe_Num(U_sym[i, j])
                    u_ij_un = Symbolics.unwrap(u_ij_num)
                    u_atomic = Symbolics.variable(:U_atomic, i, j)
                    u_bar_atomic = Symbolics.variable(:U_bar_atomic, i, j)

                    U_atomic_lookup[Symbolics.unwrap(u_atomic)] = (i, j)
                    U_bar_lookup[Symbolics.unwrap(u_bar_atomic)] = (i, j)

                    subs_dict[u_ij_un] = u_atomic

                    c_ij_un = Symbolics.unwrap(conj(u_ij_num))
                    subs_dict[c_ij_un] = u_bar_atomic

                    bc_ij_un = Symbolics.unwrap(Base.conj(u_ij_num))
                    subs_dict[bc_ij_un] = u_bar_atomic
                end
            end
        end

        matcher = LookupMatcher(U_atomic_lookup, U_bar_lookup)
        return (subs_dict, matcher, dim, :U)
    end
end

function _manual_fallback(expr, measure::HaarMeasure)
    # LazyTrace expressions are handled by fallback_integrate dispatch, not here.
    error("HaarMeasure integration failed for: $(typeof(expr))")
end

"""
    asymptotic(expr, measure::HaarMeasure, order=1)

Returns the series expansion of the integral in powers of `1/d`.
"""
function asymptotic(expr, measure::HaarMeasure, order = 1)
    d = measure.dim
    if d isa Symbolics.Num || !(d isa Integer)
        exact_res = integrate(expr, measure)
        return _expand_asymptotic(exact_res, d, order)
    end

    d_asymp = Symbolics.variable(:d_asymp)
    m_sym = dU(measure.U, d_asymp)
    exact_res = integrate(expr, m_sym)
    return _expand_asymptotic(exact_res, d_asymp, order)
end

"""
    integrate(t::LazyTrace, measure::HaarMeasure)

Integrate a product of traces of matrices over the Haar measure.
Uses the graphical Weingarten calculus.
"""
function fallback_integrate(t::LazyTrace, measure::HaarMeasure)
    dim = measure.dim

    if isempty(t.cycles)
        return t.prefactor
    end

    # Identify U and U_dag instances across ALL cycles
    U_indices = Int[]
    U_bar_indices = Int[]

    total_factors = 0
    cycle_ranges = UnitRange{Int}[]
    all_factors = SymbolicMatrix[]

    for cycle in t.cycles
        start_idx = total_factors + 1
        append!(all_factors, cycle)
        total_factors += length(cycle)
        push!(cycle_ranges, start_idx:total_factors)
    end

    for (i, f) in enumerate(all_factors)
        if f.special_type == :U
            push!(U_indices, i)
        elseif f.special_type == :U_dag
            push!(U_bar_indices, i)
        end
    end

    n_U = length(U_indices)
    n_U_bar = length(U_bar_indices)

    if n_U != n_U_bar
        return 0
    end

    # Handle case with no U/U_bar (all constant traces)
    all_slots = sort([U_indices; U_bar_indices])
    if isempty(all_slots)
        return _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    end

    # Build Wires
    wires = _build_wires(U_indices, U_bar_indices, cycle_ranges, all_factors)

    # Calculate constant part (cycles with no U/U_bar)
    constant_part = _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)

    # Weingarten Sum
    u_map = Dict{Int,Int}(idx => m for (m, idx) in enumerate(U_indices))
    ub_map = Dict{Int,Int}(idx => m for (m, idx) in enumerate(U_bar_indices))

    permutations = collect(Combinatorics.permutations(1:n_U))
    total_val = Num(0)

    for sigma in permutations
        for tau in permutations

            inv_tau = invperm(tau)
            P = [sigma[inv_tau[i]] for i = 1:n_U]
            cycle_type = get_cycle_type(P)
            wg_val = weingarten(cycle_type, dim)

            if _symbolic_isequal(wg_val, 0)
                continue
            end

            visited_sigma_U = falses(n_U)
            visited_sigma_Ub = falses(n_U)
            visited_tau_U = falses(n_U)
            visited_tau_Ub = falses(n_U)
            current_term_traces = Num[]

            # Check U cycles (sigma-based)
            for start_m = 1:n_U
                if !visited_sigma_U[start_m]
                    val = _traverse_trace_cycle(
                        start_m,
                        1,
                        sigma,
                        inv_tau,
                        wires,
                        u_map,
                        ub_map,
                        visited_sigma_U,
                        visited_sigma_Ub,
                        U_indices,
                        U_bar_indices,
                        dim,
                    )
                    push!(current_term_traces, val)
                end
            end

            # Check Ub cycles (tau-based)
            for start_m = 1:n_U_bar
                if !visited_tau_Ub[start_m]
                    val = _traverse_trace_cycle(
                        start_m,
                        2,
                        sigma,
                        inv_tau,
                        wires,
                        u_map,
                        ub_map,
                        visited_tau_U,
                        visited_tau_Ub,
                        U_indices,
                        U_bar_indices,
                        dim,
                    )
                    push!(current_term_traces, val)
                end
            end

            term_prod = isempty(current_term_traces) ? Num(1) : prod(current_term_traces)
            total_val += term_prod * wg_val
        end
    end

    return constant_part * total_val
end

function _build_wires(U_indices, U_bar_indices, cycle_ranges, all_factors)
    wires = Dict{Int,Any}()
    all_slots = sort([U_indices; U_bar_indices])
    n_slots = length(all_slots)

    for k = 1:n_slots
        start_idx = all_slots[k]

        # Determine which cycle this slot belongs to
        cycle_id = 0
        for (cid, rng) in enumerate(cycle_ranges)
            if start_idx in rng
                cycle_id = cid
                break
            end
        end
        cycle_range = cycle_ranges[cycle_id]

        # Traverse forward from start_idx until we hit another U-slot OR wrap around
        consts = SymbolicMatrix[]
        curr = start_idx

        while true
            # Move to next in cycle
            if curr == last(cycle_range)
                curr = first(cycle_range)
            else
                curr += 1
            end

            # Check if we hit a U/U_bar
            if curr in all_slots
                # Destination found
                end_idx = curr
                wires[start_idx] = (end_idx, isempty(consts) ? nothing : consts)
                break
            end

            # Otherwise it's a constant
            push!(consts, all_factors[curr])

            # Safety break if cycle is full of constants (shouldn't happen as we started from a U slot)
            if curr == start_idx
                error("Cycle should contain at least one U/U_bar")
            end
        end
    end
    return wires
end

function _evaluate_constant_cycles(t, cycle_ranges, all_slots, dim)
    constant_part = t.prefactor
    for (cid, rng) in enumerate(cycle_ranges)
        # Check if any slot in rng is in all_slots
        has_U = false
        for idx in rng
            if idx in all_slots
                has_U = true;
                break
            end
        end

        if !has_U
            cycle = t.cycles[cid]
            if isempty(cycle)
                constant_part *= dim
            else
                constant_part *= tr_val(cycle)
            end
        end
    end
    return constant_part
end

function _traverse_trace_cycle(
    start_m,
    start_type,
    sigma,
    inv_tau,
    wires,
    u_map,
    ub_map,
    visited_U,
    visited_Ub,
    U_indices,
    U_bar_indices,
    dim,
)
    curr_trace_factors = SymbolicMatrix[]
    curr_type = start_type
    curr_idx = start_m

    while true
        if curr_type == 1
            if visited_U[curr_idx]
                break
            end
            visited_U[curr_idx] = true
            next_ub_m = sigma[curr_idx]
            visited_Ub[next_ub_m] = true # Match U_i to Ub_k
            start_factor_idx = U_bar_indices[next_ub_m]
            dest_factor_idx, mat_segment = wires[start_factor_idx]
            if mat_segment !== nothing
                append!(curr_trace_factors, mat_segment)
            end
            if haskey(u_map, dest_factor_idx)
                curr_type = 1;
                curr_idx = u_map[dest_factor_idx]
            else
                curr_type = 2;
                curr_idx = ub_map[dest_factor_idx]
            end
        else
            if visited_Ub[curr_idx]
                break
            end
            visited_Ub[curr_idx] = true
            next_u_m = inv_tau[curr_idx]
            visited_U[next_u_m] = true # Match Ub_l to U_j
            start_factor_idx = U_indices[next_u_m]
            dest_factor_idx, mat_segment = wires[start_factor_idx]
            if mat_segment !== nothing
                append!(curr_trace_factors, mat_segment)
            end
            if haskey(u_map, dest_factor_idx)
                curr_type = 1;
                curr_idx = u_map[dest_factor_idx]
            else
                curr_type = 2;
                curr_idx = ub_map[dest_factor_idx]
            end
        end
    end

    if isempty(curr_trace_factors)
        return dim
    else
        return tr_val(curr_trace_factors)
    end
end
