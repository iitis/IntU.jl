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

# Quantum Information helpers
include("qi.jl")

# Pre-computed library
include("library.jl")

import LinearAlgebra: det
import Symbolics: Num
export integrate, asymptotic, dU, dPsi, dO, dSp, dGUE, dGOE, dGSE, integrate_indices, tr, det, purity, average_purity, fidelity, average_fidelity, partial_trace
export SymbolicMatrix, tr_lazy, LazyTrace

end # module
