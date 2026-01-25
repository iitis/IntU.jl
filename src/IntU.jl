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
include("haar_measure.jl")
include("pure_states.jl")

# Quantum Information helpers
include("qi.jl")

import LinearAlgebra: tr, det
import Symbolics: Num
export integrate, asymptotic, dU, dPsi, integrate_indices, tr, det, purity, average_purity, fidelity, average_fidelity, partial_trace

end # module
