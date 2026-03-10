"""
Compare IntU.jl and RTNI benchmark results side-by-side.

Usage:
    python compare_results.py

Reads results_intu.json and results_rtni.json produced by the benchmark scripts.
"""

import json
import sys

# Ordered list of (key, display_name) for the comparison table
BENCHMARKS = [
    # Easy: diagonal symbolic
    ("U_|U11|^2_sym",    "U(d), |U₁₁|², symbolic d"),
    ("U_|U11|^4_sym",    "U(d), |U₁₁|⁴, symbolic d"),
    ("U_|U11|^6_sym",    "U(d), |U₁₁|⁶, symbolic d"),
    # Harder: diagonal symbolic
    ("U_|U11|^8_sym",    "U(d), |U₁₁|⁸, symbolic d"),
    ("U_|U11|^10_sym",   "U(d), |U₁₁|¹⁰, symbolic d"),
    # Harder: numeric
    ("U_|U11|^10_d=10",  "U(d), |U₁₁|¹⁰, d=10"),
    ("U_|U11|^10_d=50",  "U(d), |U₁₁|¹⁰, d=50"),
    # Trace moments
    ("U_|trU|^4_sym",    "U(d), |tr(U)|⁴, symbolic d"),
    ("U_|trU|^6_sym",    "U(d), |tr(U)|⁶, symbolic d"),
    ("U_|trU|^8_sym",    "U(d), |tr(U)|⁸, symbolic d"),
    # Trace polynomials
    ("U_trUAUdB_sym",    "U(d), tr(UAU*B), symbolic d"),
    ("U_tr(UAUdB)^2_sym","U(d), tr((UAU*B)²), symbolic d"),
]


def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: {path} not found. Run the benchmark script first.")
        sys.exit(1)


def main():
    intu = load_json("results_intu.json")
    rtni = load_json("results_rtni.json")

    header = f"{'Integral':<35s} {'IntU.jl (ms)':>14s} {'RTNI (ms)':>14s} {'Speedup':>10s}"
    sep = "-" * len(header)

    print("\n" + "=" * len(header))
    print("IntU.jl vs RTNI (Mathematica) — Performance Comparison")
    print("=" * len(header))
    print(header)
    print(sep)

    for key, label in BENCHMARKS:
        i_data = intu.get(key, {})
        r_data = rtni.get(key, {})

        i_ms = i_data.get("median_ms")
        r_ms = r_data.get("median_ms")

        i_str = f"{i_ms:.2f}" if i_ms is not None else "N/A"
        r_str = f"{r_ms:.2f}" if r_ms is not None else "N/A"

        if i_ms is not None and r_ms is not None and i_ms > 0:
            speedup = r_ms / i_ms
            sp_str = f"{speedup:.1f}x"
        else:
            sp_str = "—"

        print(f"{label:<35s} {i_str:>14s} {r_str:>14s} {sp_str:>10s}")

    print(sep)
    print("Speedup = RTNI / IntU.jl (higher = IntU.jl is faster).")

    # Also produce a LaTeX-ready table
    print("\n\n% LaTeX table (paste into manuscript)")
    print(r"\begin{table}")
    print(r"  \centering")
    print(r"  \caption{Performance comparison: \texttt{IntU.jl} vs.\ \texttt{RTNI}")
    print(r"  (Mathematica). Median of cold-cache runs.}")
    print(r"  \label{tab:rtni_comparison}")
    print(r"  \begin{tabular}{llrrr}")
    print(r"    \hline")
    print(r"    Group & Integrand & IntU.jl (ms) & RTNI (ms) & Speedup \\")
    print(r"    \hline")

    for key, label in BENCHMARKS:
        i_data = intu.get(key, {})
        r_data = rtni.get(key, {})
        i_ms = i_data.get("median_ms")
        r_ms = r_data.get("median_ms")

        i_str = f"{i_ms:.2f}" if i_ms is not None else "---"
        r_str = f"{r_ms:.2f}" if r_ms is not None else "---"

        if i_ms is not None and r_ms is not None and i_ms > 0:
            speedup = r_ms / i_ms
            sp_str = f"{speedup:.1f}$\\times$"
        else:
            sp_str = "---"

        # Split label into group and integrand for LaTeX columns
        parts = label.split(", ", 1)
        group = parts[0]
        integrand = parts[1] if len(parts) > 1 else ""

        print(f"    {group} & {integrand} & {i_str} & {r_str} & {sp_str} \\\\")

    print(r"    \hline")
    print(r"  \end{tabular}")
    print(r"\end{table}")


if __name__ == "__main__":
    main()
