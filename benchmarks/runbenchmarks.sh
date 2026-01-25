#!/bin/bash
set -e
# Navigate to project root
cd "$(dirname "$0")/.."

echo "========================================"
echo "    Running IntU.jl Benchmarks          "
echo "========================================"

echo ""
echo ">>> Running minors.jl"
julia --project=. benchmarks/minors.jl

echo ""
echo ">>> Running pure_states.jl"
julia --project=. benchmarks/pure_states.jl

echo ""
echo ">>> Running trace_moments.jl"
julia --project=. benchmarks/trace_moments.jl

echo ""
echo "All benchmarks completed successfully."
