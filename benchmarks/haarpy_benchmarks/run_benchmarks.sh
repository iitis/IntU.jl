#!/bin/bash
# Run both benchmark suites and produce comparison table.
#
# Prerequisites:
#   conda create -n haarpy_bench python=3.11 && conda activate haarpy_bench && pip install haarpy
#   Julia with IntU.jl and JSON3 available
#   (Julia dependencies are pinned by ../Manifest.toml and instantiated below)
#
# Usage:
#   cd benchmarks/haarpy_benchmarks
#   bash run_benchmarks.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "=== Running IntU.jl benchmarks (Julia) ==="
cd "$SCRIPT_DIR"
julia --project="$BENCH_ROOT" -e 'using Pkg; Pkg.instantiate()' 2>/dev/null
julia --project="$BENCH_ROOT" bench_intu.jl

echo ""
echo "=== Running Haarpy benchmarks (Python) ==="
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate haarpy_bench
cd "$SCRIPT_DIR"
python bench_haarpy.py

echo ""
echo "=== Comparison ==="
python compare_results.py
