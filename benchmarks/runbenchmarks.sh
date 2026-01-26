#!/bin/bash
set -e
# Navigate to project root
cd "$(dirname "$0")/.."

echo "========================================"
    echo "    Running IntU.jl Benchmarks          "
echo "========================================"

# Instantiate benchmarks environment to ensure dependencies are available
julia --project=benchmarks -e 'import Pkg; Pkg.instantiate()'

for f in benchmarks/[0-9]*.jl; do
    echo ""
    echo ">>> Running $f"
    julia --project=benchmarks "$f"
done

echo ""
echo "All benchmarks completed successfully."
