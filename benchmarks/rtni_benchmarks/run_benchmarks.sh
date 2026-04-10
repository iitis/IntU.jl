#!/bin/bash
# Run IntU.jl and RTNI benchmarks and produce comparison table.
#
# Prerequisites:
#   - Julia with IntU.jl available
#     (Julia dependencies are pinned by ../Manifest.toml and instantiated below)
#   - Wolfram Mathematica with RTNI package installed
#     (RTNI.wl and precomputedWG/ must be in this directory or on $Path)
#   - For reproducibility, pin RTNI.wl to a specific source revision
#     (recommended: verify SHA-256 via RTNI_EXPECTED_SHA256 environment variable)
#
# Usage:
#   cd benchmarks/rtni_benchmarks
#   bash run_benchmarks.sh
#
# If Mathematica is on a different machine:
#   1. Run bench_rtni.wl on the Mathematica machine:
#        math -script bench_rtni.wl
#   2. Copy results_rtni.json to this directory
#   3. Run bench_intu.jl on the Julia machine:
#        julia --project=.. bench_intu.jl
#   4. Compare:
#        python compare_results.py

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BENCH_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "=== Running IntU.jl benchmarks (Julia) ==="
cd "$SCRIPT_DIR"
julia --project="$BENCH_ROOT" -e 'using Pkg; Pkg.instantiate()' 2>/dev/null
julia --project="$BENCH_ROOT" bench_intu.jl

echo ""
echo "=== Running RTNI benchmarks (Mathematica) ==="
cd "$SCRIPT_DIR"
if [ -n "${RTNI_EXPECTED_SHA256:-}" ] && [ -f "RTNI.wl" ]; then
    actual_sha="$(sha256sum RTNI.wl | awk '{print $1}')"
    if [ "$actual_sha" != "$RTNI_EXPECTED_SHA256" ]; then
        echo "ERROR: RTNI.wl SHA-256 mismatch."
        echo "  expected: $RTNI_EXPECTED_SHA256"
        echo "  actual:   $actual_sha"
        exit 1
    fi
fi
if command -v math &> /dev/null; then
    math -script bench_rtni.wl
else
    echo "Mathematica not found. Run bench_rtni.wl on a machine with Mathematica"
    echo "and copy results_rtni.json to this directory."
    if [ ! -f results_rtni.json ]; then
        echo "ERROR: results_rtni.json not found. Cannot produce comparison."
        exit 1
    fi
fi

echo ""
echo "=== Comparison ==="
python3 compare_results.py
