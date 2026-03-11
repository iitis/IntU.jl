"""
Compare IntU.jl and RTNI benchmark results side-by-side.

Usage:
    python compare_results.py

Reads results_intu.json and results_rtni.json produced by the benchmark scripts.
"""

import json
import sys
from pathlib import Path

# Ordered list of benchmark rows.
# Each row can map an IntU key to a different RTNI key when comparing
# alternative IntU implementations of the same integral.
BENCHMARKS = [
    # Element API: diagonal moments
    {"intu_key": "U_|U11|^2_sym", "rtni_key": "U_|U11|^2_sym", "label": "U(d) (Element API), |U₁₁|², symbolic d"},
    {"intu_key": "U_|U11|^4_sym", "rtni_key": "U_|U11|^4_sym", "label": "U(d) (Element API), |U₁₁|⁴, symbolic d"},
    {"intu_key": "U_|U11|^6_sym", "rtni_key": "U_|U11|^6_sym", "label": "U(d) (Element API), |U₁₁|⁶, symbolic d"},
    {"intu_key": "U_|U11|^8_sym", "rtni_key": "U_|U11|^8_sym", "label": "U(d) (Element API), |U₁₁|⁸, symbolic d"},
    {"intu_key": "U_|U11|^10_sym", "rtni_key": "U_|U11|^10_sym", "label": "U(d) (Element API), |U₁₁|¹⁰, symbolic d"},
    {"intu_key": "U_|U11|^10_d=10", "rtni_key": "U_|U11|^10_d=10", "label": "U(d) (Element API), |U₁₁|¹⁰, d=10"},
    {"intu_key": "U_|U11|^10_d=50", "rtni_key": "U_|U11|^10_d=50", "label": "U(d) (Element API), |U₁₁|¹⁰, d=50"},
    # ITensors graphical engine: same integrals as above
    {"intu_key": "U_|U11|^2_sym_itensor", "rtni_key": "U_|U11|^2_sym", "label": "U(d) (ITensors), |U₁₁|², symbolic d"},
    {"intu_key": "U_|U11|^4_sym_itensor", "rtni_key": "U_|U11|^4_sym", "label": "U(d) (ITensors), |U₁₁|⁴, symbolic d"},
    {"intu_key": "U_|U11|^6_sym_itensor", "rtni_key": "U_|U11|^6_sym", "label": "U(d) (ITensors), |U₁₁|⁶, symbolic d"},
    {"intu_key": "U_|U11|^8_sym_itensor", "rtni_key": "U_|U11|^8_sym", "label": "U(d) (ITensors), |U₁₁|⁸, symbolic d"},
    {"intu_key": "U_|U11|^10_sym_itensor", "rtni_key": "U_|U11|^10_sym", "label": "U(d) (ITensors), |U₁₁|¹⁰, symbolic d"},
    {"intu_key": "U_|U11|^10_d=10_itensor", "rtni_key": "U_|U11|^10_d=10", "label": "U(d) (ITensors), |U₁₁|¹⁰, d=10"},
    {"intu_key": "U_|U11|^10_d=50_itensor", "rtni_key": "U_|U11|^10_d=50", "label": "U(d) (ITensors), |U₁₁|¹⁰, d=50"},
    # Trace moments
    {"intu_key": "U_|trU|^4_sym", "rtni_key": "U_|trU|^4_sym", "label": "U(d), |tr(U)|⁴, symbolic d"},
    {"intu_key": "U_|trU|^6_sym", "rtni_key": "U_|trU|^6_sym", "label": "U(d), |tr(U)|⁶, symbolic d"},
    {"intu_key": "U_|trU|^8_sym", "rtni_key": "U_|trU|^8_sym", "label": "U(d), |tr(U)|⁸, symbolic d"},
    {"intu_key": "U_|trU|^4_sym_itensor", "rtni_key": None, "label": "U(d) (ITensors), |tr(U)|⁴, symbolic d"},
    {"intu_key": "U_|trU|^6_sym_itensor", "rtni_key": None, "label": "U(d) (ITensors), |tr(U)|⁶, symbolic d"},
    {"intu_key": "U_|trU|^8_sym_itensor", "rtni_key": None, "label": "U(d) (ITensors), |tr(U)|⁸, symbolic d"},
    # Trace polynomials
    {"intu_key": "U_trUAUdB_sym", "rtni_key": "U_trUAUdB_sym", "label": "U(d), tr(UAU*B), symbolic d"},
    {"intu_key": "U_tr(UAUdB)^2_sym", "rtni_key": "U_tr(UAUdB)^2_sym", "label": "U(d), tr((UAU*B)²), symbolic d"},
    {"intu_key": "U_trUAUdB_sym_itensor", "rtni_key": None, "label": "U(d) (ITensors), tr(UAU*B), symbolic d"},
    {"intu_key": "U_tr(UAUdB)^2_sym_itensor", "rtni_key": None, "label": "U(d) (ITensors), tr((UAU*B)²), symbolic d"},
]

TRACE_SCALAR_ROWS = {
    "U_|trU|^4_sym",
    "U_|trU|^6_sym",
    "U_|trU|^8_sym",
}

LATEX_ROW_LABELS = {
    "U_|U11|^2_sym": ("U($d$) (Element API)", r"$|U_{11}|^2$, symbolic $d$"),
    "U_|U11|^4_sym": ("U($d$) (Element API)", r"$|U_{11}|^4$, symbolic $d$"),
    "U_|U11|^6_sym": ("U($d$) (Element API)", r"$|U_{11}|^6$, symbolic $d$"),
    "U_|U11|^8_sym": ("U($d$) (Element API)", r"$|U_{11}|^8$, symbolic $d$"),
    "U_|U11|^10_sym": ("U($d$) (Element API)", r"$|U_{11}|^{10}$, symbolic $d$"),
    "U_|U11|^10_d=10": ("U($d$) (Element API)", r"$|U_{11}|^{10}$, $d=10$"),
    "U_|U11|^10_d=50": ("U($d$) (Element API)", r"$|U_{11}|^{10}$, $d=50$"),
    "U_|U11|^2_sym_itensor": ("U($d$) (ITensors)", r"$|U_{11}|^2$, symbolic $d$"),
    "U_|U11|^4_sym_itensor": ("U($d$) (ITensors)", r"$|U_{11}|^4$, symbolic $d$"),
    "U_|U11|^6_sym_itensor": ("U($d$) (ITensors)", r"$|U_{11}|^6$, symbolic $d$"),
    "U_|U11|^8_sym_itensor": ("U($d$) (ITensors)", r"$|U_{11}|^8$, symbolic $d$"),
    "U_|U11|^10_sym_itensor": ("U($d$) (ITensors)", r"$|U_{11}|^{10}$, symbolic $d$"),
    "U_|U11|^10_d=10_itensor": ("U($d$) (ITensors)", r"$|U_{11}|^{10}$, $d=10$"),
    "U_|U11|^10_d=50_itensor": ("U($d$) (ITensors)", r"$|U_{11}|^{10}$, $d=50$"),
    "U_|trU|^4_sym": ("U($d$)", r"$|tr(U)|^4$, symbolic $d$"),
    "U_|trU|^6_sym": ("U($d$)", r"$|tr(U)|^6$, symbolic $d$"),
    "U_|trU|^8_sym": ("U($d$)", r"$|tr(U)|^8$, symbolic $d$"),
    "U_|trU|^4_sym_itensor": ("U($d$) (ITensors)", r"$|tr(U)|^4$, symbolic $d$"),
    "U_|trU|^6_sym_itensor": ("U($d$) (ITensors)", r"$|tr(U)|^6$, symbolic $d$"),
    "U_|trU|^8_sym_itensor": ("U($d$) (ITensors)", r"$|tr(U)|^8$, symbolic $d$"),
    "U_trUAUdB_sym": ("U($d$)", r"$\mathrm{tr}(UAU^\ast B)$, symbolic $d$"),
    "U_tr(UAUdB)^2_sym": ("U($d$)", r"$\mathrm{tr}((UAU^\ast B)^2)$, symbolic $d$"),
    "U_trUAUdB_sym_itensor": ("U($d$) (ITensors)", r"$\mathrm{tr}(UAU^\ast B)$, symbolic $d$"),
    "U_tr(UAUdB)^2_sym_itensor": ("U($d$) (ITensors)", r"$\mathrm{tr}((UAU^\ast B)^2)$, symbolic $d$"),
}
LATEX_ROWS_OUTPUT = Path("rtni_table_rows.tex")


def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: {path} not found. Run the benchmark script first.")
        sys.exit(1)

def rtni_row_is_scalarized(rtni_key, r_data):
    if not rtni_key:
        return False
    if rtni_key not in TRACE_SCALAR_ROWS:
        return True
    if not r_data:
        return False
    raw = str(r_data.get("result", "")).strip()
    # RTNI graph placeholder in our existing outputs for trace moments.
    if raw == "{}" or raw == "":
        return False
    return True


def row_metrics(row, intu, rtni):
    i_data = intu.get(row["intu_key"], {})
    rtni_key = row.get("rtni_key")
    r_data = rtni.get(rtni_key, {}) if rtni_key else {}
    i_ms = i_data.get("median_ms")
    r_ms = r_data.get("median_ms") if rtni_row_is_scalarized(rtni_key, r_data) else None
    return i_ms, r_ms


def speedup_string(i_ms, r_ms, latex=False):
    if i_ms is not None and r_ms is not None and i_ms > 0:
        speedup = r_ms / i_ms
        return f"{speedup:.1f}$\\times$" if latex else f"{speedup:.1f}x"
    return "---" if latex else "—"


def render_console_table(intu, rtni):
    header = f"{'Integral':<35s} {'IntU.jl (ms)':>14s} {'RTNI (ms)':>14s} {'Speedup':>10s}"
    sep = "-" * len(header)

    print("\n" + "=" * len(header))
    print("IntU.jl vs RTNI (Mathematica) — Performance Comparison")
    print("=" * len(header))
    print(header)
    print(sep)

    for row in BENCHMARKS:
        i_ms, r_ms = row_metrics(row, intu, rtni)
        i_str = f"{i_ms:.2f}" if i_ms is not None else "N/A"
        r_str = f"{r_ms:.2f}" if r_ms is not None else "N/A"
        sp_str = speedup_string(i_ms, r_ms, latex=False)
        print(f"{row['label']:<35s} {i_str:>14s} {r_str:>14s} {sp_str:>10s}")

    print(sep)
    print("Speedup = RTNI / IntU.jl (higher = IntU.jl is faster).")


def render_latex_table_stdout(intu, rtni):
    print("\n\n% LaTeX table (paste into manuscript)")
    print(r"\begin{table}")
    print(r"  \centering")
    print(r"  \caption{Performance comparison: \texttt{IntU.jl} vs.\ \texttt{RTNI}")
    print(r"  (Mathematica). Median runtime over repeated runs.}")
    print(r"  \label{tab:rtni_comparison}")
    print(r"  \begin{tabular}{llrrr}")
    print(r"    \hline")
    print(r"    Group & Integrand & IntU.jl (ms) & RTNI (ms) & Speedup \\")
    print(r"    \hline")

    for row in BENCHMARKS:
        i_ms, r_ms = row_metrics(row, intu, rtni)
        i_str = f"{i_ms:.2f}" if i_ms is not None else "---"
        r_str = f"{r_ms:.2f}" if r_ms is not None else "---"
        sp_str = speedup_string(i_ms, r_ms, latex=True)

        parts = row["label"].split(", ", 1)
        group = parts[0]
        integrand = parts[1] if len(parts) > 1 else ""
        print(f"    {group} & {integrand} & {i_str} & {r_str} & {sp_str} \\\\")

    print(r"    \hline")
    print(r"  \end{tabular}")
    print(r"\end{table}")


def write_latex_rows_file(intu, rtni, out_path):
    missing = [row["intu_key"] for row in BENCHMARKS if row["intu_key"] not in LATEX_ROW_LABELS]
    if missing:
        raise KeyError(f"Missing LaTeX labels for keys: {missing}")

    lines = [
        "% Auto-generated by benchmarks/rtni_benchmarks/compare_results.py",
        "% Do not edit manually.",
    ]
    for row in BENCHMARKS:
        key = row["intu_key"]
        i_ms, r_ms = row_metrics(row, intu, rtni)
        i_str = f"{i_ms:.2f}" if i_ms is not None else "---"
        r_str = f"{r_ms:.2f}" if r_ms is not None else "---"
        sp_str = speedup_string(i_ms, r_ms, latex=True)
        group, integrand = LATEX_ROW_LABELS[key]
        lines.append(f"{group} & {integrand} & {i_str} & {r_str} & {sp_str} \\\\")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    intu = load_json("results_intu.json")
    rtni = load_json("results_rtni.json")
    render_console_table(intu, rtni)
    render_latex_table_stdout(intu, rtni)
    write_latex_rows_file(intu, rtni, LATEX_ROWS_OUTPUT)
    print(f"\nLaTeX rows saved to {LATEX_ROWS_OUTPUT}")


if __name__ == "__main__":
    main()
