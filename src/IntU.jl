module IntU

using Symbolics
using Combinatorics
using LinearAlgebra
using SymbolicUtils
using Memoization
using DataStructures
using MacroTools
using ProgressMeter


# Core Weingarten logic
include("weingarten.jl")

# Symbolic Trace Logic
include("symbolic_trace.jl")
include("centered_perm_impl.jl")

# Shared trace-graph helpers (used by Haar, Gaussian, Ginibre measures)
include("trace_helpers.jl")

# Symbolic utilities and index matching
include("symbolic_utils.jl")

# Core type definitions and integrate dispatches
include("integration_core.jl")

# Integration engine (rewrite rules, _integrate_core, process_term)
include("integration_engine.jl")

# Integration rule registrations and low-level integrate_indices_* functions
include("integration_rules.jl")

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
include("integrate_macro.jl")

# Public API: asymptotic, evaluate (after all measures for _reconstruct_symbolic)
include("integration_api.jl")


# Quantum Information helpers
include("qi.jl")

# Pre-computed library
include("library.jl")

import LinearAlgebra: det
import Symbolics: Num

# --- Exports ---

# Core integration API
export integrate, evaluate, asymptotic, @integrate

# Measures
export dU, dSU, SpecialUnitary
export dPsi
export dO, dSp
export dGUE, dGOE, dGSE
export dGinUE, dGinOE, dGinSE
export dDiagUnitary, dStiefel
export dDesign
export dPerm, dCPerm
export dCOE, dCSE, dCUE

# HCIZ
export hciz, vandermonde_det

# Weingarten calculus
export integrate_indices
export weingarten, weingarten_orthogonal_val, weingarten_symplectic_val
export get_pair_partitions,
    canonicalize_pair_partition,
    conjugate_partition,
    count_loops,
    murnaghan_nakayama,
    character_at_id,
    irrep_dimension,
    get_weingarten_orthogonal_data

# Quantum information helpers
export tr, det, partial_trace

# Internal API (exported for advanced use)
export AbstractIndexMatcher, MetadataMatcher, _integrate_core, process_term

# Symbolic matrix types
export SymbolicMatrix, tr_lazy, LazyTrace, LazySum
export symbolic_unitary,
    symbolic_orthogonal, symbolic_symplectic, symbolic_pure_state, symbolic_permutation
export integrate_graphical, GraphicalUnitary, ITensorUnitary

end # module
