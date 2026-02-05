#!/bin/bash
set -e
# Navigate to project root
cd "$(dirname "$0")/.."

echo "========================================"
echo "    Running IntU.jl Examples            "
echo "========================================"

for f in examples/*.jl; do
    echo ""
    echo ">>> Running $f"
    julia --project="examples/" "$f"
done

echo ""
echo "All examples completed successfully."
