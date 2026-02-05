#!/bin/bash
set -e
# Navigate to project root
cd "$(dirname "$0")/.."

echo "========================================"
    echo "    Running IntU.jl Benchmarks          "
echo "========================================"

# Instantiate benchmarks environment to ensure dependencies are available
julia --project=benchmarks -e 'import Pkg; Pkg.instantiate()'

# Run specific scripts if provided, otherwise run all benchmarks
if [ $# -gt 0 ]; then
    for arg in "$@"; do
        found=0
        # Check if arg is a file, or if benchmarks/arg*.jl exists
        for f in "$arg" benchmarks/"$arg"*.jl; do
            if [ -f "$f" ]; then
                echo ""
                echo ">>> Running $f"
                julia --project=benchmarks "$f"
                found=1
                break
            fi
        done
        if [ $found -eq 0 ]; then
            echo "Warning: No benchmark matching '$arg' found."
        fi
    done
else
    for f in benchmarks/[0-9]*.jl; do
        echo ""
        echo ">>> Running $f"
        julia --project=benchmarks "$f"
    done
fi

echo ""
echo "All benchmarks completed successfully."
