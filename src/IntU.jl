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

# Symbolic Trace Logic
include("symbolic_trace.jl")

# Quantum Information helpers
include("qi.jl")

import LinearAlgebra: tr, det
import Symbolics: Num
export integrate, asymptotic, dU, dPsi, integrate_indices, tr, det, purity, average_purity, fidelity, average_fidelity, partial_trace
export SymbolicMatrix, tr_lazy, LazyTrace

end # module
