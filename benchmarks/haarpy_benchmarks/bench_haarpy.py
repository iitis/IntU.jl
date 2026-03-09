"""
Performance comparison benchmarks for Haarpy.
Computes the same integrals as bench_intu.jl for a head-to-head comparison.

Usage:
    conda activate haarpy_bench
    python bench_haarpy.py
"""

import time
import statistics
import json
import sys
from sympy import Symbol

import haarpy

d = Symbol("d")

N_WARMUP = 2
N_SAMPLES = 30
# If a single call takes longer than this, reduce to fewer samples
SLOW_THRESHOLD_S = 1.0
SLOW_SAMPLES = 5


def benchmark(func, n_warmup=N_WARMUP, n_samples=N_SAMPLES):
    """Benchmark a function, returning median time in seconds."""
    # Warmup and detect slow benchmarks
    for _ in range(n_warmup):
        func()

    # Probe: time one call to decide sample count
    start = time.perf_counter()
    func()
    probe_time = time.perf_counter() - start

    if probe_time > SLOW_THRESHOLD_S:
        n_samples = SLOW_SAMPLES

    times = []
    for _ in range(n_samples):
        # Clear lru_cache between samples for fair comparison
        # (IntU.jl benchmarks also re-run from scratch each sample)
        haarpy.haar_integral_unitary.cache_clear()
        haarpy.haar_integral_orthogonal.cache_clear()
        haarpy.haar_integral_circular_orthogonal.cache_clear()
        haarpy.haar_integral_permutation.cache_clear()

        start = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        "median_s": statistics.median(times),
        "min_s": min(times),
        "max_s": max(times),
        "samples": n_samples,
    }


# ============================================================================
# Define benchmark cases
# ============================================================================

results = {}


def run_and_report(name, func):
    print(f"  Running: {name} ...", end=" ", flush=True)
    try:
        # Verify it produces a result first
        res = func()
        stats = benchmark(func)
        ms = stats["median_s"] * 1000
        n = stats["samples"]
        print(f"{ms:.2f} ms  (N={n}, result: {res})")
        results[name] = {**stats, "median_ms": ms, "result": str(res)}
    except Exception as e:
        print(f"FAILED: {e}")
        results[name] = {"error": str(e)}


# --- Section 1: Unitary |U_11|^{2k}, symbolic d ---
# This is the primary comparison requested by the reviewer (k=3,4,5)

print("\n=== Unitary: |U_11|^{2k}, symbolic d ===")

run_and_report(
    "U_|U11|^6_sym",
    lambda: haarpy.haar_integral_unitary(
        ((1, 1, 1), (1, 1, 1), (1, 1, 1), (1, 1, 1)), d
    ),
)

run_and_report(
    "U_|U11|^8_sym",
    lambda: haarpy.haar_integral_unitary(
        ((1, 1, 1, 1), (1, 1, 1, 1), (1, 1, 1, 1), (1, 1, 1, 1)), d
    ),
)

run_and_report(
    "U_|U11|^10_sym",
    lambda: haarpy.haar_integral_unitary(
        (
            (1, 1, 1, 1, 1),
            (1, 1, 1, 1, 1),
            (1, 1, 1, 1, 1),
            (1, 1, 1, 1, 1),
        ),
        d,
    ),
)


# --- Section 2: Unitary |U_11|^{2k}, numeric d ---
print("\n=== Unitary: |U_11|^{2k}, numeric d ===")

for d_val in [10, 50]:
    run_and_report(
        f"U_|U11|^10_d={d_val}",
        lambda dv=d_val: haarpy.haar_integral_unitary(
            (
                (1, 1, 1, 1, 1),
                (1, 1, 1, 1, 1),
                (1, 1, 1, 1, 1),
                (1, 1, 1, 1, 1),
            ),
            dv,
        ),
    )


# --- Section 3: Orthogonal O_11^k, symbolic d ---
print("\n=== Orthogonal: O_11^k, symbolic d ===")

run_and_report(
    "O_O11^2_sym",
    lambda: haarpy.haar_integral_orthogonal(((1, 1), (1, 1)), d),
)

run_and_report(
    "O_O11^4_sym",
    lambda: haarpy.haar_integral_orthogonal(((1, 1, 1, 1), (1, 1, 1, 1)), d),
)


# --- Section 4: Orthogonal O_11^k, numeric d ---
print("\n=== Orthogonal: O_11^k, numeric d ===")

run_and_report(
    "O_O11^6_d=10",
    lambda: haarpy.haar_integral_orthogonal(
        ((1, 1, 1, 1, 1, 1), (1, 1, 1, 1, 1, 1)), 10
    ),
)

run_and_report(
    "O_O11^8_d=20",
    lambda: haarpy.haar_integral_orthogonal(
        ((1, 1, 1, 1, 1, 1, 1, 1), (1, 1, 1, 1, 1, 1, 1, 1)), 20
    ),
)

run_and_report(
    "O_O11^10_d=20",
    lambda: haarpy.haar_integral_orthogonal(
        (
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        ),
        20,
    ),
)

run_and_report(
    "O_O11^10_d=50",
    lambda: haarpy.haar_integral_orthogonal(
        (
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        ),
        50,
    ),
)


# --- Section 5: COE |S_11|^{2k}, symbolic d ---
print("\n=== Circular Orthogonal (COE): |S_11|^{2k}, symbolic d ===")

run_and_report(
    "COE_|S11|^2_sym",
    lambda: haarpy.haar_integral_circular_orthogonal(((1, 1), (1, 1)), d),
)

run_and_report(
    "COE_|S11|^4_sym",
    lambda: haarpy.haar_integral_circular_orthogonal(
        ((1, 1, 1, 1), (1, 1, 1, 1)), d
    ),
)

run_and_report(
    "COE_|S11|^6_sym",
    lambda: haarpy.haar_integral_circular_orthogonal(
        ((1, 1, 1, 1, 1, 1), (1, 1, 1, 1, 1, 1)), d
    ),
)


# --- Section 6: Permutation P_11^k ---
print("\n=== Permutation: P_11^k ===")

run_and_report(
    "Perm_P11^10_d=100",
    lambda: haarpy.haar_integral_permutation(
        (1, 1, 1, 1, 1, 1, 1, 1, 1, 1), (1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 100
    ),
)


# ============================================================================
# Save results
# ============================================================================
print("\n" + "=" * 72)
print("Summary (median times in ms)")
print("=" * 72)
print(f"{'Benchmark':<30s} {'Median (ms)':>12s}")
print("-" * 42)
for name, data in results.items():
    if "error" in data:
        print(f"{name:<30s} {'FAILED':>12s}")
    else:
        print(f"{name:<30s} {data['median_ms']:12.2f}")

with open("results_haarpy.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to results_haarpy.json")
