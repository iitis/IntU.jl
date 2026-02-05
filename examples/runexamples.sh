#!/bin/bash
set -e
# Navigate to project root
cd "$(dirname "$0")/.."

echo "========================================"
echo "    Running IntU.jl Examples            "
echo "========================================"

# Run specific scripts if provided, otherwise run all examples
if [ $# -gt 0 ]; then
    for arg in "$@"; do
        found=0
        # Check direct file or examples/ arg*.jl
        for f in "$arg" examples/"$arg"*.jl; do
            if [ -f "$f" ]; then
                echo ""
                echo ">>> Running $f"
                julia --project="examples/" "$f"
                found=1
                break
            fi
        done
        if [ $found -eq 0 ]; then
            echo "Warning: No example matching '$arg' found."
        fi
    done
else
    for f in examples/*.jl; do
        echo ""
        echo ">>> Running $f"
        julia --project="examples/" "$f"
    done
fi

echo ""
echo "All examples completed successfully."
