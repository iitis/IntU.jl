"""
Compare IntU.jl and Haarpy benchmark results side-by-side.

Usage:
    python compare_results.py

Reads results_intu.json and results_haarpy.json produced by the benchmark scripts.
"""

import json
import sys

BENCHMARKS = [
    # Unitary symbolic
    ("U_|U11|^6_sym", "U(d), |U₁₁|⁶, symbolic d"),
    ("U_|U11|^8_sym", "U(d), |U₁₁|⁸, symbolic d"),
    ("U_|U11|^10_sym", "U(d), |U₁₁|¹⁰, symbolic d"),
    # Unitary numeric
    ("U_|U11|^10_d=10", "U(d), |U₁₁|¹⁰, d=10"),
    ("U_|U11|^10_d=50", "U(d), |U₁₁|¹⁰, d=50"),
    # Orthogonal symbolic
    ("O_O11^2_sym", "O(d), O₁₁², symbolic d"),
    ("O_O11^4_sym", "O(d), O₁₁⁴, symbolic d"),
    # Orthogonal numeric
    ("O_O11^6_d=10", "O(d), O₁₁⁶, d=10"),
    ("O_O11^8_d=20", "O(d), O₁₁⁸, d=20"),
    ("O_O11^10_d=20", "O(d), O₁₁¹⁰, d=20"),
    ("O_O11^10_d=50", "O(d), O₁₁¹⁰, d=50"),
    # COE
    ("COE_|S11|^2_sym", "COE, |S₁₁|², symbolic d"),
    ("COE_|S11|^4_sym", "COE, |S₁₁|⁴, symbolic d"),
    ("COE_|S11|^6_sym", "COE, |S₁₁|⁶, symbolic d"),
    # Off-diagonal: Unitary
    ("U_offdiag_4_sym", "U(d), |U₁₁|²|U₁₂|², symbolic d"),
    ("U_offdiag_8_sym", "U(d), |U₁₁|⁴|U₁₂|⁴, symbolic d"),
    ("U_cross_4_sym", "U(d), |U₁₁|²|U₂₂|², symbolic d"),
    # Off-diagonal: Orthogonal
    ("O_offdiag_4_sym", "O(d), O₁₁²O₁₂², symbolic d"),
    ("O_cross_4_sym", "O(d), O₁₁O₁₂O₂₁O₂₂, symbolic d"),
    # Off-diagonal: COE
    ("COE_offdiag_2_sym", "COE, |S₁₂|², symbolic d"),
    ("COE_offdiag_4_sym", "COE, |S₁₂|⁴, symbolic d"),
    ("COE_mixed_4_sym", "COE, |S₁₁|²|S₁₂|², symbolic d"),
    # Permutation
    ("Perm_P11^10_d=100", "Perm, P₁₁¹⁰, d=100"),
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
    haarpy = load_json("results_haarpy.json")

    header = f"{'Integral':<35s} {'IntU.jl (ms)':>14s} {'Haarpy (ms)':>14s} {'Speedup':>10s}"
    sep = "-" * len(header)

    print("\n" + "=" * len(header))
    print("IntU.jl vs Haarpy — Performance Comparison")
    print("=" * len(header))
    print(header)
    print(sep)

    for key, label in BENCHMARKS:
        i_data = intu.get(key, {})
        h_data = haarpy.get(key, {})

        i_ms = i_data.get("median_ms")
        h_ms = h_data.get("median_ms")

        i_str = f"{i_ms:.2f}" if i_ms is not None else "N/A"
        h_str = f"{h_ms:.2f}" if h_ms is not None else "N/A"

        if i_ms is not None and h_ms is not None and i_ms > 0:
            speedup = h_ms / i_ms
            sp_str = f"{speedup:.1f}x"
        else:
            sp_str = "—"

        print(f"{label:<35s} {i_str:>14s} {h_str:>14s} {sp_str:>10s}")

    print(sep)
    print(f"N = {N_SAMPLES} samples, median reported. Speedup = Haarpy / IntU.jl.")

    print("\n\n% LaTeX table (paste into manuscript)")
    print(r"\begin{table}")
    print(r"  \centering")
    print(r"  \caption{Performance comparison: IntU.jl vs.\ Haarpy (median of 30 runs).}")
    print(r"  \label{tab:haarpy_comparison}")
    print(r"  \begin{tabular}{llrrr}")
    print(r"    \hline")
    print(r"    Group & Integrand & IntU.jl (ms) & Haarpy (ms) & Speedup \\")
    print(r"    \hline")

    for key, label in BENCHMARKS:
        i_data = intu.get(key, {})
        h_data = haarpy.get(key, {})
        i_ms = i_data.get("median_ms")
        h_ms = h_data.get("median_ms")

        i_str = f"{i_ms:.2f}" if i_ms is not None else "---"
        h_str = f"{h_ms:.2f}" if h_ms is not None else "---"

        if i_ms is not None and h_ms is not None and i_ms > 0:
            speedup = h_ms / i_ms
            sp_str = f"{speedup:.1f}$\\times$"
        else:
            sp_str = "---"

        parts = label.split(", ", 1)
        group = parts[0]
        integrand = parts[1] if len(parts) > 1 else ""

        print(f"    {group} & {integrand} & {i_str} & {h_str} & {sp_str} \\\\")

    print(r"    \hline")
    print(r"  \end{tabular}")
    print(r"\end{table}")


N_SAMPLES = 30

if __name__ == "__main__":
    main()
