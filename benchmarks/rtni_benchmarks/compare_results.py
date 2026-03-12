"""
Compare IntU.jl and RTNI benchmark results side-by-side.

Usage:
    python compare_results.py

Reads results_intu.json and results_rtni.json produced by the benchmark scripts.
"""

import json
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def row(*, intu_key, rtni_key, label, group, integrand, trace_row=False):
    return {
        "intu_key": intu_key,
        "rtni_key": rtni_key,
        "label": label,
        "group": group,
        "integrand": integrand,
        "trace_row": trace_row,
    }


MAIN_COMPARABLE_ROWS = [
    row(
        intu_key="U_|U11|^2_sym",
        rtni_key="U_|U11|^2_sym",
        label="U(d), |U₁₁|², symbolic d",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^2$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^4_sym",
        rtni_key="U_|U11|^4_sym",
        label="U(d), |U₁₁|⁴, symbolic d",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^4$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^6_sym",
        rtni_key="U_|U11|^6_sym",
        label="U(d), |U₁₁|⁶, symbolic d",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^6$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^8_sym",
        rtni_key="U_|U11|^8_sym",
        label="U(d), |U₁₁|⁸, symbolic d",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^8$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^8_d=10",
        rtni_key="U_|U11|^8_d=10",
        label="U(d), |U₁₁|⁸, d=10",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^8$, $d=10$",
    ),
    row(
        intu_key="U_|U11|^2|U12|^2_sym",
        rtni_key="U_|U11|^2|U12|^2_sym",
        label="U(d), |U₁₁|²|U₁₂|², symbolic d",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^2|U_{12}|^2$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^2|U12|^4_sym",
        rtni_key="U_|U11|^2|U12|^4_sym",
        label="U(d), |U₁₁|²|U₁₂|⁴, symbolic d",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^2|U_{12}|^4$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^2|U22|^2_sym",
        rtni_key="U_|U11|^2|U22|^2_sym",
        label="U(d), |U₁₁|²|U₂₂|², symbolic d",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^2|U_{22}|^2$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^2|U12|^2_d=10",
        rtni_key="U_|U11|^2|U12|^2_d=10",
        label="U(d), |U₁₁|²|U₁₂|², d=10",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^2|U_{12}|^2$, $d=10$",
    ),
    row(
        intu_key="U_|U11|^2|U12|^4_d=10",
        rtni_key="U_|U11|^2|U12|^4_d=10",
        label="U(d), |U₁₁|²|U₁₂|⁴, d=10",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^2|U_{12}|^4$, $d=10$",
    ),
    row(
        intu_key="U_|U11|^2|U22|^2_d=10",
        rtni_key="U_|U11|^2|U22|^2_d=10",
        label="U(d), |U₁₁|²|U₂₂|², d=10",
        group="U($d$) (Element API)",
        integrand=r"$|U_{11}|^2|U_{22}|^2$, $d=10$",
    ),
    row(
        intu_key="U_|trU|^4_sym",
        rtni_key="U_|trU|^4_sym",
        label="U(d), |tr(U)|⁴, symbolic d",
        group="U($d$)",
        integrand=r"$|tr(U)|^4$, symbolic $d$",
        trace_row=True,
    ),
    row(
        intu_key="U_|trU|^6_sym",
        rtni_key="U_|trU|^6_sym",
        label="U(d), |tr(U)|⁶, symbolic d",
        group="U($d$)",
        integrand=r"$|tr(U)|^6$, symbolic $d$",
        trace_row=True,
    ),
    row(
        intu_key="U_|trU|^8_sym",
        rtni_key="U_|trU|^8_sym",
        label="U(d), |tr(U)|⁸, symbolic d",
        group="U($d$)",
        integrand=r"$|tr(U)|^8$, symbolic $d$",
        trace_row=True,
    ),
    row(
        intu_key="U_trUAUdB_sym",
        rtni_key="U_trUAUdB_sym",
        label="U(d), tr(UAU*B), symbolic d",
        group="U($d$)",
        integrand=r"$\mathrm{tr}(UAU^\ast B)$, symbolic $d$",
    ),
    row(
        intu_key="U_tr(UAUdB)^2_sym",
        rtni_key="U_tr(UAUdB)^2_sym",
        label="U(d), tr((UAU*B)²), symbolic d",
        group="U($d$)",
        integrand=r"$\mathrm{tr}((UAU^\ast B)^2)$, symbolic $d$",
    ),
]

SUPPLEMENTARY_ROWS = [
    row(
        intu_key="U_|U11|^2_sym_itensor",
        rtni_key=None,
        label="U(d) (ITensors), |U₁₁|², symbolic d",
        group="U($d$) (ITensors)",
        integrand=r"$|U_{11}|^2$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^4_sym_itensor",
        rtni_key=None,
        label="U(d) (ITensors), |U₁₁|⁴, symbolic d",
        group="U($d$) (ITensors)",
        integrand=r"$|U_{11}|^4$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^6_sym_itensor",
        rtni_key=None,
        label="U(d) (ITensors), |U₁₁|⁶, symbolic d",
        group="U($d$) (ITensors)",
        integrand=r"$|U_{11}|^6$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^8_sym_itensor",
        rtni_key=None,
        label="U(d) (ITensors), |U₁₁|⁸, symbolic d",
        group="U($d$) (ITensors)",
        integrand=r"$|U_{11}|^8$, symbolic $d$",
    ),
    row(
        intu_key="U_|U11|^8_d=10_itensor",
        rtni_key=None,
        label="U(d) (ITensors), |U₁₁|⁸, d=10",
        group="U($d$) (ITensors)",
        integrand=r"$|U_{11}|^8$, $d=10$",
    ),
    row(
        intu_key="U_|trU|^4_sym_itensor",
        rtni_key=None,
        label="U(d) (ITensors), |tr(U)|⁴, symbolic d",
        group="U($d$) (ITensors)",
        integrand=r"$|tr(U)|^4$, symbolic $d$",
        trace_row=True,
    ),
    row(
        intu_key="U_|trU|^6_sym_itensor",
        rtni_key=None,
        label="U(d) (ITensors), |tr(U)|⁶, symbolic d",
        group="U($d$) (ITensors)",
        integrand=r"$|tr(U)|^6$, symbolic $d$",
        trace_row=True,
    ),
    row(
        intu_key="U_|trU|^8_sym_itensor",
        rtni_key=None,
        label="U(d) (ITensors), |tr(U)|⁸, symbolic d",
        group="U($d$) (ITensors)",
        integrand=r"$|tr(U)|^8$, symbolic $d$",
        trace_row=True,
    ),
    row(
        intu_key="U_trUAUdB_sym_itensor",
        rtni_key=None,
        label="U(d) (ITensors), tr(UAU*B), symbolic d",
        group="U($d$) (ITensors)",
        integrand=r"$\mathrm{tr}(UAU^\ast B)$, symbolic $d$",
    ),
    row(
        intu_key="U_tr(UAUdB)^2_sym_itensor",
        rtni_key=None,
        label="U(d) (ITensors), tr((UAU*B)²), symbolic d",
        group="U($d$) (ITensors)",
        integrand=r"$\mathrm{tr}((UAU^\ast B)^2)$, symbolic $d$",
    ),
]

LATEX_ROWS_OUTPUT = SCRIPT_DIR / "rtni_table_rows.tex"


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
    script = meta.get("script", "unknown")
    sources = meta.get("sources", {})

    print(f"  timestamp_utc: {ts}")
    print(f"  runtime: {runtime_name} {runtime_ver}")
    print(f"  host: {host} ({os_name}, {arch})")
    if isinstance(sources, dict) and sources:
        print(f"  sources: {sources}")
    else:
        print("  sources: unknown")
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


def warn_meta_consistency(intu_meta, rtni_meta):
    if not intu_meta or not rtni_meta:
        print("Warning: metadata consistency check skipped because one or both _meta blocks are missing.")
        return

    for path, label in [
        (("host", "hostname"), "host.hostname"),
        (("host", "os"), "host.os"),
        (("host", "arch"), "host.arch"),
    ]:
        warn_mismatch_field(intu_meta, rtni_meta, path, label, "IntU", "RTNI")

    for path, label, side in [
        (("runtime", "name"), "runtime.name", "IntU"),
        (("runtime", "version"), "runtime.version", "IntU"),
        (("packages", "IntU"), "packages.IntU", "IntU"),
        (("packages", "ITensors"), "packages.ITensors", "IntU"),
    ]:
        warn_missing_meta_field(intu_meta, path, label, side)

    for path, label, side in [
        (("runtime", "name"), "runtime.name", "RTNI"),
        (("runtime", "version"), "runtime.version", "RTNI"),
        (("packages", "RTNI"), "packages.RTNI", "RTNI"),
    ]:
        warn_missing_meta_field(rtni_meta, path, label, side)

    intu_runtime = nested_get(intu_meta, ("runtime", "name"))
    if intu_runtime is not None and intu_runtime != "Julia":
        print(f"Warning: unexpected IntU runtime.name={intu_runtime!r} (expected 'Julia').")

    rtni_runtime = nested_get(rtni_meta, ("runtime", "name"))
    if rtni_runtime is not None and rtni_runtime != "Mathematica":
        print(
            f"Warning: unexpected RTNI runtime.name={rtni_runtime!r} "
            "(expected 'Mathematica')."
        )

    for path, label in [
        (("sources", "RTNI.wl", "sha256"), "sources.RTNI.wl.sha256"),
        (("sources", "RTNI.wl", "path"), "sources.RTNI.wl.path"),
    ]:
        warn_missing_meta_field(intu_meta, path, label, "IntU")
        warn_missing_meta_field(rtni_meta, path, label, "RTNI")
        warn_mismatch_field(intu_meta, rtni_meta, path, label, "IntU", "RTNI")


def get_entry(blob, key):
    if not key:
        return {}
    data = blob.get(key)
    return data if isinstance(data, dict) else {}


def get_status(entry):
    status = entry.get("status")
    if isinstance(status, str):
        return status
    return "ok" if isinstance(entry.get("median_ms"), (int, float)) else "missing"


def get_median_ms(entry):
    value = entry.get("median_ms")
    return float(value) if isinstance(value, (int, float)) else None


def is_non_scalar_rtni(entry):
    status = get_status(entry)
    raw = str(entry.get("result", "")).strip()
    if status.startswith("non_scalar"):
        return True
    return raw == "{}"


def speedup_str(i_ms, r_ms, *, latex=False, non_scalar=False):
    if i_ms is not None and r_ms is not None and i_ms > 0:
        value = r_ms / i_ms
        if latex:
            s = f"{value:.1f}$\\times$"
            return s + r" \textsuperscript{\ensuremath{\dagger}}" if non_scalar else s
        s = f"{value:.1f}x"
        return s + "†" if non_scalar else s
    return "---" if latex else "—"


def row_metrics(row, intu, rtni):
    i_entry = get_entry(intu, row["intu_key"])
    r_entry = get_entry(rtni, row["rtni_key"])
    i_ms = get_median_ms(i_entry)
    r_ms = get_median_ms(r_entry)
    return {
        "i_ms": i_ms,
        "r_ms": r_ms,
        "i_status": get_status(i_entry),
        "r_status": get_status(r_entry),
        "rtni_non_scalar": bool(row.get("trace_row")) and is_non_scalar_rtni(r_entry),
    }


def render_console_table(rows, intu, rtni, *, title):
    header = f"{'Integral':<40s} {'IntU.jl (ms)':>14s} {'RTNI (ms)':>14s} {'Speedup':>13s}"
    sep = "-" * len(header)

    print("\n" + "=" * len(header))
    print(title)
    print("=" * len(header))
    print(header)
    print(sep)

    for r in rows:
        m = row_metrics(r, intu, rtni)
        i_str = f"{m['i_ms']:.2f}" if m["i_ms"] is not None else "N/A"
        r_str = f"{m['r_ms']:.2f}" if m["r_ms"] is not None else "N/A"
        sp = speedup_str(
            m["i_ms"],
            m["r_ms"],
            latex=False,
            non_scalar=m["rtni_non_scalar"],
        )
        print(f"{r['label']:<40s} {i_str:>14s} {r_str:>14s} {sp:>13s}")

    print(sep)
    print("Speedup = RTNI / IntU.jl (higher = IntU.jl is faster).")
    print("† RTNI row timed successfully, but result remained non-scalar (graph placeholder).")


def render_supplementary_console(rows, intu):
    header = f"{'Supplementary IntU-only row':<45s} {'IntU.jl (ms)':>14s} {'Status':>12s}"
    sep = "-" * len(header)

    print("\n" + "=" * len(header))
    print("Supplementary IntU Rows (Not in Main RTNI Table)")
    print("=" * len(header))
    print(header)
    print(sep)
    for r in rows:
        entry = get_entry(intu, r["intu_key"])
        ms = get_median_ms(entry)
        status = get_status(entry)
        ms_str = f"{ms:.2f}" if ms is not None else "N/A"
        print(f"{r['label']:<45s} {ms_str:>14s} {status:>12s}")
    print(sep)


def write_latex_rows_file(rows, intu, rtni, out_path):
    lines = [
        "% Auto-generated by benchmarks/rtni_benchmarks/compare_results.py",
        "% Main comparable Element-API rows only.",
        r"\providecommand{\RTNIComparisonRows}{}",
        r"\renewcommand{\RTNIComparisonRows}{%",
    ]

    for r in rows:
        m = row_metrics(r, intu, rtni)
        i_str = f"{m['i_ms']:.2f}" if m["i_ms"] is not None else "---"
        r_str = f"{m['r_ms']:.2f}" if m["r_ms"] is not None else "---"
        sp_str = speedup_str(
            m["i_ms"],
            m["r_ms"],
            latex=True,
            non_scalar=m["rtni_non_scalar"],
        )
        lines.append(
            f"{r['group']} & {r['integrand']} & {i_str} & {r_str} & {sp_str} \\\\%"
        )

    lines.append("}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    intu = load_json("results_intu.json")
    rtni = load_json("results_rtni.json")
    intu_meta = get_meta(intu)
    rtni_meta = get_meta(rtni)

    print_meta("IntU.jl", intu_meta)
    print_meta("RTNI", rtni_meta)
    warn_meta_consistency(intu_meta, rtni_meta)

    render_console_table(
        MAIN_COMPARABLE_ROWS,
        intu,
        rtni,
        title="IntU.jl vs RTNI (Mathematica) — Main Comparable Rows",
    )
    render_supplementary_console(SUPPLEMENTARY_ROWS, intu)
    write_latex_rows_file(MAIN_COMPARABLE_ROWS, intu, rtni, LATEX_ROWS_OUTPUT)
    print(f"\nLaTeX rows saved to {LATEX_ROWS_OUTPUT}")


if __name__ == "__main__":
    main()
