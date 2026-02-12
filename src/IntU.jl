module IntU

using Symbolics
using Combinatorics
using LinearAlgebra
using SymbolicUtils
using Memoization
using DataStructures

# Fix for ambiguity in Num(::Complex)
function Symbolics.Num(z::Complex{<:Real})
    return Symbolics.Term(complex, [z.re, z.im])
end

function Symbolics.Num(z::Complex{Num})
    return z.re + im*z.im
end

# Core Weingarten logic
include("weingarten.jl")

# Symbolic Trace Logic
include("symbolic_trace.jl")

# Shared integration engine
include("integration_core.jl")

# Specific measures
include("haar_measure.jl")
include("pure_states.jl")
include("real_measures.jl")
include("gaussian_measures.jl")
include("unitary_designs.jl")
include("circular_measures.jl")
include("permutation_measures.jl")
include("itensors_integration.jl")
include("hciz.jl")
include("su_measure.jl")
include("diagonal_unitary.jl")
include("stiefel_measure.jl")
include("macros.jl")


# Quantum Information helpers
include("qi.jl")

# Pre-computed library
include("library.jl")

import LinearAlgebra: det
import Symbolics: Num
export integrate,
    asymptotic,
    dU,
    hciz,
    vandermonde_det,
    dSU,
    SpecialUnitary,
    dPsi,
    dO,
    dSp,
    dGUE,
    dGOE,
    dGSE,
    dGinUE,
    dGinOE,
    dGinSE,
    dDiagUnitary,
    dStiefel,
    integrate_indices,
    tr,
    det,
    purity,
    average_purity,
    fidelity,
    average_fidelity,
    partial_trace,
    dDesign,
    @symbolic_dimension,
    @integrate
export dPerm, dCPerm
export dCOE, dCSE, dCUE
export weingarten, weingarten_orthogonal_val, weingarten_symplectic_val
export get_pair_partitions,
    canonicalize_pair_partition,
    conjugate_partition,
    count_loops,
    murnaghan_nakayama,
    character_at_id,
    irrep_dimension,
    get_weingarten_orthogonal_data
export AbstractIndexMatcher, LookupMatcher, _integrate_core, process_term

export SymbolicMatrix, tr_lazy, LazyTrace, LazySum
export integrate_graphical, GraphicalUnitary, ITensorUnitary

end # module
