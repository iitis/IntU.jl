#!/bin/bash
# Run both benchmark suites and produce comparison table.
#
# Prerequisites:
#   conda env create -f environment.haarpy_bench.yml
#   Julia with IntU.jl and JSON3 available
#   (Julia dependencies are pinned by ../Manifest.toml and instantiated below)
#
# Usage:
#   cd benchmarks/haarpy_benchmarks
#   bash run_benchmarks.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
HAARPY_VERSION="${HAARPY_VERSION:-0.0.6}"
HAARPY_CONDA_ENV="${HAARPY_CONDA_ENV:-haarpy_bench}"
HAARPY_ENV_YML="${HAARPY_ENV_YML:-$SCRIPT_DIR/environment.haarpy_bench.yml}"

echo "=== Running IntU.jl benchmarks (Julia) ==="
cd "$SCRIPT_DIR"
julia --project="$BENCH_ROOT" -e 'using Pkg; Pkg.instantiate()' 2>/dev/null
julia --project="$BENCH_ROOT" bench_intu.jl

echo ""
echo "=== Running Haarpy benchmarks (Python) ==="
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda is required for Haarpy benchmarks."
    echo "Install conda, then create env with:"
    echo "  conda env create -f \"$HAARPY_ENV_YML\""
    exit 1
fi
eval "$(conda shell.bash hook 2>/dev/null)"
if ! conda env list | awk '{print $1}' | grep -Fxq "$HAARPY_CONDA_ENV"; then
    echo "ERROR: conda environment '$HAARPY_CONDA_ENV' not found."
    echo "Create it with:"
    echo "  conda env create -f \"$HAARPY_ENV_YML\""
    exit 1
fi
conda activate "$HAARPY_CONDA_ENV"
cd "$SCRIPT_DIR"
python - <<'PY'
import os, sys
from importlib import metadata

expected = os.environ.get("HAARPY_VERSION", "0.0.6")
try:
    actual = metadata.version("haarpy")
except Exception:
    print(f"ERROR: haarpy is not installed in the active environment (expected {expected}).")
    print(f"Install with: pip install 'haarpy=={expected}'")
    sys.exit(1)

if actual != expected:
    print(f"ERROR: haarpy version mismatch: found {actual}, expected {expected}.")
    print(f"Install with: pip install 'haarpy=={expected}'")
    sys.exit(1)

print(f"Using haarpy=={actual}")
PY
python bench_haarpy.py

echo ""
echo "=== Comparison ==="
python compare_results.py
