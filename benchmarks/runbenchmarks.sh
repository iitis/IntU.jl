#!/bin/bash
set -e
# Navigate to project root
cd "$(dirname "$0")/.."

echo "========================================"
echo "    Running IntU.jl Benchmarks          "
echo "========================================"

echo ""
echo ">>> Running minors.jl"
julia --project=benchmarks benchmarks/minors.jl

echo ""
echo ">>> Running pure_states.jl"
julia --project=benchmarks benchmarks/pure_states.jl

echo ""
echo ">>> Running trace_moments.jl"
julia --project=benchmarks benchmarks/trace_moments.jl

echo ""
echo ">>> Running symbolic_trace.jl"
julia --project=benchmarks benchmarks/symbolic_trace.jl

echo ""
echo "All benchmarks completed successfully."
