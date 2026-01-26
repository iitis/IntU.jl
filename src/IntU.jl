module IntU

using Symbolics
using Combinatorics
using LinearAlgebra
using SymbolicUtils
using Memoization

# Core Weingarten logic
include("weingarten.jl")

# Shared integration engine
include("integration_core.jl")

# Specific measures
# Specific measures
include("haar_measure.jl")
include("pure_states.jl")
include("real_measures.jl")
include("gaussian_measures.jl")

# Symbolic Trace Logic
include("symbolic_trace.jl")

# Quantum Information helpers
include("qi.jl")

import LinearAlgebra: det
import Symbolics: Num
export integrate, asymptotic, dU, dPsi, dO, dSp, dGUE, integrate_indices, tr, det, purity, average_purity, fidelity, average_fidelity, partial_trace
export SymbolicMatrix, tr_lazy, LazyTrace

end # module
