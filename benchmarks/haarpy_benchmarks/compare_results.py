"""
Compare IntegrateUnitary.jl and Haarpy benchmark results side-by-side.

Usage:
    python compare_results.py

Reads results_intu.json and results_haarpy.json produced by the benchmark scripts.
"""

import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent


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
    path_obj = SCRIPT_DIR / path
    try:
        with path_obj.open(encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: {path_obj} not found. Run the benchmark script first.")
        sys.exit(1)


def get_meta(blob):
    raw = blob.get("_meta")
    return raw if isinstance(raw, dict) else None


def nested_get(dct, keys):
    cur = dct
    for key in keys:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    return cur


def print_meta(label, meta):
    print(f"\n{label} metadata")
    if not meta:
        print("  warning: missing _meta block")
        return

    ts = meta.get("timestamp_utc", "unknown")
    runtime_name = nested_get(meta, ("runtime", "name")) or "unknown"
    runtime_ver = nested_get(meta, ("runtime", "version")) or "unknown"
    host = nested_get(meta, ("host", "hostname")) or "unknown"
    os_name = nested_get(meta, ("host", "os")) or "unknown"
    arch = nested_get(meta, ("host", "arch")) or "unknown"
    packages = meta.get("packages", {})
    script = meta.get("script", "unknown")

    print(f"  timestamp_utc: {ts}")
    print(f"  runtime: {runtime_name} {runtime_ver}")
    print(f"  host: {host} ({os_name}, {arch})")
    print(f"  packages: {packages if isinstance(packages, dict) else 'unknown'}")
    print(f"  script: {script}")


def warn_missing_meta_field(meta, path, label, side):
    if nested_get(meta, path) is None:
        print(f"Warning: metadata missing {label} in {side} result.")


def warn_mismatch_field(left_meta, right_meta, path, label, left_side, right_side):
    left_val = nested_get(left_meta, path)
    right_val = nested_get(right_meta, path)
    if left_val is None or right_val is None:
        return
    if left_val != right_val:
        print(
            f"Warning: metadata mismatch for {label}: "
            f"{left_side}={left_val!r}, {right_side}={right_val!r}."
        )


def warn_meta_consistency(intu_meta, haarpy_meta):
    if not intu_meta or not haarpy_meta:
        print("Warning: metadata consistency check skipped because one or both _meta blocks are missing.")
        return

    for path, label in [
        (("host", "hostname"), "host.hostname"),
        (("host", "os"), "host.os"),
        (("host", "arch"), "host.arch"),
    ]:
        warn_mismatch_field(intu_meta, haarpy_meta, path, label, "IntegrateUnitary", "Haarpy")

    for path, label, side in [
        (("runtime", "name"), "runtime.name", "IntegrateUnitary"),
        (("runtime", "version"), "runtime.version", "IntegrateUnitary"),
        (("packages", "IntegrateUnitary"), "packages.IntegrateUnitary", "IntegrateUnitary"),
    ]:
        warn_missing_meta_field(intu_meta, path, label, side)

    for path, label, side in [
        (("runtime", "name"), "runtime.name", "Haarpy"),
        (("runtime", "version"), "runtime.version", "Haarpy"),
        (("packages", "haarpy"), "packages.haarpy", "Haarpy"),
    ]:
        warn_missing_meta_field(haarpy_meta, path, label, side)

    intu_runtime = nested_get(intu_meta, ("runtime", "name"))
    if intu_runtime is not None and intu_runtime != "Julia":
        print(f"Warning: unexpected IntegrateUnitary runtime.name={intu_runtime!r} (expected 'Julia').")

    haarpy_runtime = nested_get(haarpy_meta, ("runtime", "name"))
    if haarpy_runtime is not None and haarpy_runtime != "Python":
        print(f"Warning: unexpected Haarpy runtime.name={haarpy_runtime!r} (expected 'Python').")


def main():
    intu = load_json("results_intu.json")
    haarpy = load_json("results_haarpy.json")
    intu_meta = get_meta(intu)
    haarpy_meta = get_meta(haarpy)

    print_meta("IntegrateUnitary.jl", intu_meta)
    print_meta("Haarpy", haarpy_meta)
    warn_meta_consistency(intu_meta, haarpy_meta)

    header = f"{'Integral':<35s} {'IntegrateUnitary.jl (ms)':>14s} {'Haarpy (ms)':>14s} {'Speedup':>10s}"
    sep = "-" * len(header)

    print("\n" + "=" * len(header))
    print("IntegrateUnitary.jl vs Haarpy — Performance Comparison")
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
    print("IntegrateUnitary.jl rows: fixed N=30 cold-cache samples (median reported).")
    print("Haarpy rows: default N=30, adaptively reduced to N=5 for slow cases (median reported).")
    print("Speedup = Haarpy / IntegrateUnitary.jl.")

    print("\n\n% LaTeX table (paste into manuscript)")
    print(r"\begin{table}")
    print(r"  \centering")
    print(
        r"  \caption{Performance comparison: IntegrateUnitary.jl vs.\ Haarpy (IntegrateUnitary.jl: fixed $N=30$ cold-cache samples; Haarpy: default $N=30$ adaptively reduced to $N=5$ for slow rows).}"
    )
    print(r"  \label{tab:haarpy_comparison}")
    print(r"  \begin{tabular}{llrrr}")
    print(r"    \hline")
    print(r"    Group & Integrand & IntegrateUnitary.jl (ms) & Haarpy (ms) & Speedup \\")
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


if __name__ == "__main__":
    main()
