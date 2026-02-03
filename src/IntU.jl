module IntU

using Symbolics
using Combinatorics
using LinearAlgebra
using SymbolicUtils
using Memoization

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
include("itensors_integration.jl")


# Quantum Information helpers
include("qi.jl")

# Pre-computed library
include("library.jl")

import LinearAlgebra: det
import Symbolics: Num
export integrate, asymptotic, dU, dPsi, dO, dSp, dGUE, dGOE, dGSE, integrate_indices, tr, det, purity, average_purity, fidelity, average_fidelity, partial_trace, dDesign, symbolic_dimension_unitary, @symbolic_dimension
export weingarten, weingarten_orthogonal_val, weingarten_symplectic_val
export get_pair_partitions, canonicalize_pair_partition, conjugate_partition, count_loops, murnaghan_nakayama, character_at_id, irrep_dimension, get_weingarten_orthogonal_data
export AbstractIndexMatcher, LookupMatcher, _integrate_core, SymbolicUnitary, process_term

export SymbolicMatrix, tr_lazy, LazyTrace, LazySum
export integrate_graphical, GraphicalUnitary, ITensorUnitary

end # module
