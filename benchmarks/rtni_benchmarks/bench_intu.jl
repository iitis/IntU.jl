"""
Performance comparison benchmarks for IntU.jl.
Computes the same unitary integrals as bench_rtni.wl for a head-to-head
comparison with RTNI (Mathematica).

Usage:
    julia --project=/path/to/IntU.jl/benchmarks bench_intu.jl
"""

using IntU
using Symbolics
using BenchmarkTools
using Memoization
using Printf
using ITensors
using Dates
using SHA

@variables d

# Match RTNI benchmark settings
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10
const FAST_SAMPLES = 10
const SLOW_SAMPLES = 3
const SLOW_THRESHOLD_S = 2.0
const RATIONAL_TOL_ABS = BigFloat("1e-70")
const RATIONAL_TOL_REL = BigFloat("1e-60")
const RATIONAL_DEN_MAX = big(10)^12

# High-order symbolic ITensor contractions can trigger large-order warnings;
# disable warning spam so benchmark output remains readable.
ITensors.disable_warn_order()

function median_ms(b)
    m = median(b)
    return m.time / 1e6  # ns -> ms
end

results = Dict{String,Any}()

function module_version(mod)
    try
        v = Base.pkgversion(mod)
        return v === nothing ? "unknown" : string(v)
    catch
        return "unknown"
    end
end

function file_sha256(path::AbstractString)
    if !isfile(path)
        return "missing"
    end
    open(path, "r") do io
        return bytes2hex(SHA.sha256(read(io)))
    end
end

function benchmark_meta()
    hostname = try
        gethostname()
    catch
        "unknown"
    end
    rtni_path = abspath(joinpath(@__DIR__, "RTNI.wl"))
    return Dict(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "host" => Dict(
            "hostname" => hostname,
            "os" => string(Sys.KERNEL),
            "arch" => string(Sys.ARCH),
            "machine" => string(Sys.MACHINE),
        ),
        "runtime" => Dict("name" => "Julia", "version" => string(VERSION)),
        "packages" => Dict(
            "IntU" => module_version(IntU),
            "ITensors" => module_version(ITensors),
        ),
        "sources" => Dict(
            "RTNI.wl" => Dict("path" => rtni_path, "sha256" => file_sha256(rtni_path)),
        ),
        "script" => abspath(@__FILE__),
    )
end

function json_escape(s::AbstractString)
    return replace(
        s,
        "\\" => "\\\\",
        "\"" => "\\\"",
        "\n" => "\\n",
        "\r" => "\\r",
        "\t" => "\\t",
    )
end

function to_json(x)
    if x isa Dict
        parts = String[]
        for (k, v) in sort(collect(x), by = p -> string(p[1]))
            push!(parts, "\"" * json_escape(string(k)) * "\": " * to_json(v))
        end
        return "{" * join(parts, ", ") * "}"
    elseif x isa AbstractVector
        return "[" * join((to_json(v) for v in x), ", ") * "]"
    elseif x isa AbstractString
        return "\"" * json_escape(x) * "\""
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x === nothing
        return "null"
    else
        return string(x)
    end
end

function canonicalize_result(res)
    # ITensors graphical integration returns rank-0 ITensor for scalar diagrams.
    if res isa ITensor
        if length(inds(res)) == 0
            return simplify(ITensors.scalar(res))
        end
        return res
    end
    return simplify(res)
end

function maybe_rationalize_numeric(x::BigFloat)
    tol = max(RATIONAL_TOL_ABS, abs(x) * RATIONAL_TOL_REL)
    q = rationalize(BigInt, x; tol = tol)
    if denominator(q) > RATIONAL_DEN_MAX
        return nothing
    end
    if abs(x - BigFloat(q)) <= tol
        return q
    end
    return nothing
end

function numeric_result_string(x::BigFloat, fallback::String)
    q = maybe_rationalize_numeric(x)
    if q === nothing
        return fallback
    end
    if denominator(q) == 1
        return string(numerator(q))
    end
    return string(q)
end

function format_result_string(res)
    if res isa Union{Integer,Rational}
        return string(res)
    end

    if res isa AbstractFloat
        s = string(res)
        return numeric_result_string(BigFloat(res), s)
    end

    # Keep symbolic expressions untouched; only parse plain numeric strings.
    s = string(res)
    x = try
        parse(BigFloat, s)
    catch
        nothing
    end
    if x === nothing
        return s
    end

    return numeric_result_string(x, s)
end

function projector_e11(idx_left::Index, idx_right::Index)
    t = ITensor(idx_left, idx_right)
    t[idx_left => 1, idx_right => 1] = 1
    return t
end

function u11_moment_itensor_network(k::Int, idx_dim::Int)
    tensors = Any[]
    for n in 1:k
        out_idx = Index(idx_dim, "u_out,$n")
        in_idx = Index(idx_dim, "u_in,$n")
        out_idx_adj = Index(idx_dim, "uadj_out,$n")
        in_idx_adj = Index(idx_dim, "uadj_in,$n")

        U = ITensorUnitary(out_indices = [out_idx], in_indices = [in_idx], is_adj = false)
        U_dag = ITensorUnitary(
            out_indices = [out_idx_adj],
            in_indices = [in_idx_adj],
            is_adj = true,
        )

        # |U_11|^(2k): each U/U† pair is pinned by rank-1 projectors E11 on input and output legs.
        push!(tensors, U, projector_e11(in_idx, out_idx_adj), U_dag, projector_e11(in_idx_adj, out_idx))
    end
    return tensors
end

function integrate_u11_itensor(k::Int, measure::IntU.AbstractMeasure; idx_dim::Int)
    tensors = u11_moment_itensor_network(k, idx_dim)
    return integrate(tensors, measure)
end

function matrix_constant_tensor(vals::AbstractMatrix{<:Real}, idx_left::Index, idx_right::Index)
    t = ITensor(idx_left, idx_right)
    for a in 1:size(vals, 1), b in 1:size(vals, 2)
        t[idx_left => a, idx_right => b] = vals[a, b]
    end
    return t
end

function trace_moment_itensor_network(k::Int, idx_dim::Int)
    tensors = Any[]
    for n in 1:k
        out_idx = Index(idx_dim, "trU_out,$n")
        in_idx = Index(idx_dim, "trU_in,$n")
        out_idx_adj = Index(idx_dim, "trUb_out,$n")
        in_idx_adj = Index(idx_dim, "trUb_in,$n")

        U = ITensorUnitary(out_indices = [out_idx], in_indices = [in_idx], is_adj = false)
        U_dag = ITensorUnitary(
            out_indices = [out_idx_adj],
            in_indices = [in_idx_adj],
            is_adj = true,
        )

        # Close each U and U† into a trace loop.
        push!(tensors, U, delta(out_idx, in_idx), U_dag, delta(out_idx_adj, in_idx_adj))
    end
    return tensors
end

function integrate_trace_moment_itensor(k::Int, measure::IntU.AbstractMeasure; idx_dim::Int)
    tensors = trace_moment_itensor_network(k, idx_dim)
    return integrate(tensors, measure)
end

function tr_uau_db_itensor_network(idx_dim::Int, avals::AbstractMatrix{<:Real}, bvals::AbstractMatrix{<:Real})
    i = Index(idx_dim, "u_out,1")
    j = Index(idx_dim, "u_in,1")
    i2 = Index(idx_dim, "ud_in,1")
    j2 = Index(idx_dim, "ud_out,1")
    U = ITensorUnitary(out_indices = [i], in_indices = [j], is_adj = false)
    U_dag = ITensorUnitary(out_indices = [j2], in_indices = [i2], is_adj = true)
    A = matrix_constant_tensor(avals, j, j2)
    B = matrix_constant_tensor(bvals, i2, i)
    return Any[U, A, U_dag, B]
end

function tr_uau_db_sq_itensor_network(
    idx_dim::Int,
    avals::AbstractMatrix{<:Real},
    bvals::AbstractMatrix{<:Real},
)
    i1 = Index(idx_dim, "u1_out")
    j1 = Index(idx_dim, "u1_in")
    j2 = Index(idx_dim, "ud1_out")
    i2 = Index(idx_dim, "ud1_in")
    i3 = Index(idx_dim, "u2_out")
    j3 = Index(idx_dim, "u2_in")
    j4 = Index(idx_dim, "ud2_out")
    i4 = Index(idx_dim, "ud2_in")

    U1 = ITensorUnitary(out_indices = [i1], in_indices = [j1], is_adj = false)
    U1_dag = ITensorUnitary(out_indices = [j2], in_indices = [i2], is_adj = true)
    U2 = ITensorUnitary(out_indices = [i3], in_indices = [j3], is_adj = false)
    U2_dag = ITensorUnitary(out_indices = [j4], in_indices = [i4], is_adj = true)

    A1 = matrix_constant_tensor(avals, j1, j2)
    B1 = matrix_constant_tensor(bvals, i2, i3)
    A2 = matrix_constant_tensor(avals, j3, j4)
    B2 = matrix_constant_tensor(bvals, i4, i1)

    return Any[U1, A1, U1_dag, B1, U2, A2, U2_dag, B2]
end

function integrate_tr_uau_db_itensor(
    measure::IntU.AbstractMeasure;
    idx_dim::Int,
    avals::AbstractMatrix{<:Real},
    bvals::AbstractMatrix{<:Real},
)
    return integrate(tr_uau_db_itensor_network(idx_dim, avals, bvals), measure)
end

function integrate_tr_uau_db_sq_itensor(
    measure::IntU.AbstractMeasure;
    idx_dim::Int,
    avals::AbstractMatrix{<:Real},
    bvals::AbstractMatrix{<:Real},
)
    return integrate(tr_uau_db_sq_itensor_network(idx_dim, avals, bvals), measure)
end

function choose_samples(f)
    # Probe once with cache reset to mimic cold-ish run behavior and
    # dynamically reduce samples for very slow benchmarks.
    Memoization.empty_all_caches!()
    probe_s = @elapsed f()
    return probe_s > SLOW_THRESHOLD_S ? SLOW_SAMPLES : FAST_SAMPLES
end

function run_and_report(name, f)
    print("  Running: $name ...")
    # Verify it works
    res = f()
    res = canonicalize_result(res)
    res_str = format_result_string(res)
    nsamp = choose_samples(f)
    b = @benchmark $f() evals=1 samples=nsamp setup=(Memoization.empty_all_caches!())
    ms = median_ms(b)
    @printf(" %.2f ms  (N=%d, result: %s)\n", ms, length(b.times), res_str)
    results[name] = Dict(
        "median_ms" => ms,
        "result" => res_str,
        "samples" => length(b.times),
        "status" => "ok",
    )
    return ms
end

# ============================================================================
# Section 1: Easy - Unitary |U_11|^{2k}, symbolic d (k=1,2,3)
# ============================================================================
println("\n=== Easy: Unitary |U_11|^{2k}, symbolic d ===")

U = SymbolicMatrix(:U, :U)
measure_sym = dU(d)

run_and_report("U_|U11|^2_sym", () -> integrate(abs(U[1, 1])^2, measure_sym))
run_and_report("U_|U11|^4_sym", () -> integrate(abs(U[1, 1])^4, measure_sym))
run_and_report("U_|U11|^6_sym", () -> integrate(abs(U[1, 1])^6, measure_sym))

# ============================================================================
# Section 2: Harder - Unitary |U_11|^{2k}, symbolic d (k=4)
# ============================================================================
println("\n=== Harder: Unitary |U_11|^{2k}, symbolic d ===")

run_and_report("U_|U11|^8_sym", () -> integrate(abs(U[1,1])^8, measure_sym))

# ============================================================================
# Section 2b: Same moments via ITensors graphical engine
# ============================================================================
println("\n=== Harder: Unitary |U_11|^{2k}, symbolic d (ITensors engine) ===")
println("    (all benchmarked k; uses rank-1 projector network)")

const ITENSOR_U11_IDX_DIM = 1

run_and_report("U_|U11|^2_sym_itensor", () -> integrate_u11_itensor(1, measure_sym; idx_dim = ITENSOR_U11_IDX_DIM))
run_and_report("U_|U11|^4_sym_itensor", () -> integrate_u11_itensor(2, measure_sym; idx_dim = ITENSOR_U11_IDX_DIM))
run_and_report("U_|U11|^6_sym_itensor", () -> integrate_u11_itensor(3, measure_sym; idx_dim = ITENSOR_U11_IDX_DIM))
run_and_report("U_|U11|^8_sym_itensor", () -> integrate_u11_itensor(4, measure_sym; idx_dim = ITENSOR_U11_IDX_DIM))

# ============================================================================
# Section 3: Practical numeric coverage (no 10th-power rows)
# ============================================================================
println("\n=== Harder: Unitary |U_11|^{2k}, numeric d ===")

for d_val in [10]
    U_n = SymbolicMatrix(:U, :U, d_val)
    m_n = dU(d_val)
    run_and_report("U_|U11|^8_d=$d_val", () -> integrate(abs(U_n[1,1])^8, m_n))
    run_and_report("U_|U11|^8_d=$(d_val)_itensor", () -> integrate_u11_itensor(4, m_n; idx_dim = ITENSOR_U11_IDX_DIM))
end

# ============================================================================
# Section 4: Mixed element moments
# ============================================================================
println("\n=== Mixed element moments, symbolic d ===")

run_and_report(
    "U_|U11|^2|U12|^2_sym",
    () -> integrate(abs(U[1, 1])^2 * abs(U[1, 2])^2, measure_sym),
)
run_and_report(
    "U_|U11|^2|U12|^4_sym",
    () -> integrate(abs(U[1, 1])^2 * abs(U[1, 2])^4, measure_sym),
)
run_and_report(
    "U_|U11|^2|U22|^2_sym",
    () -> integrate(abs(U[1, 1])^2 * abs(U[2, 2])^2, measure_sym),
)

println("\n=== Mixed element moments, numeric d ===")
for d_val in [10]
    U_n = SymbolicMatrix(:U, :U, d_val)
    m_n = dU(d_val)
    run_and_report(
        "U_|U11|^2|U12|^2_d=$d_val",
        () -> integrate(abs(U_n[1, 1])^2 * abs(U_n[1, 2])^2, m_n),
    )
    run_and_report(
        "U_|U11|^2|U12|^4_d=$d_val",
        () -> integrate(abs(U_n[1, 1])^2 * abs(U_n[1, 2])^4, m_n),
    )
    run_and_report(
        "U_|U11|^2|U22|^2_d=$d_val",
        () -> integrate(abs(U_n[1, 1])^2 * abs(U_n[2, 2])^2, m_n),
    )
end

# ============================================================================
# Section 5: Trace moments (graph-based in RTNI)
# ============================================================================
println("\n=== Trace moments: |tr(U)|^{2k}, symbolic d ===")

run_and_report("U_|trU|^4_sym", () -> integrate(abs(tr(U))^4, measure_sym))
run_and_report("U_|trU|^6_sym", () -> integrate(abs(tr(U))^6, measure_sym))
run_and_report("U_|trU|^8_sym", () -> integrate(abs(tr(U))^8, measure_sym))

println("\n=== Trace moments: |tr(U)|^{2k}, symbolic d (ITensors engine) ===")
println("    (fixed concrete trace-loop dimension for ITensor constants)")

const ITENSOR_TRACE_IDX_DIM = 2
run_and_report("U_|trU|^4_sym_itensor", () -> integrate_trace_moment_itensor(2, measure_sym; idx_dim = ITENSOR_TRACE_IDX_DIM))
run_and_report("U_|trU|^6_sym_itensor", () -> integrate_trace_moment_itensor(3, measure_sym; idx_dim = ITENSOR_TRACE_IDX_DIM))
run_and_report("U_|trU|^8_sym_itensor", () -> integrate_trace_moment_itensor(4, measure_sym; idx_dim = ITENSOR_TRACE_IDX_DIM))

# ============================================================================
# Section 6: Trace polynomial - tr(U A U^* B)
# ============================================================================
println("\n=== Trace polynomials: symbolic d ===")

A = SymbolicMatrix(:A, :Constant)
B = SymbolicMatrix(:B, :Constant)

run_and_report("U_trUAUdB_sym", () -> integrate(tr(U * A * U' * B), measure_sym))

run_and_report("U_tr(UAUdB)^2_sym",
    () -> integrate(tr(U * A * U' * B * U * A * U' * B), measure_sym))

println("\n=== Trace polynomials: symbolic d (ITensors engine, concrete A/B) ===")

const ITENSOR_POLY_IDX_DIM = 2
const A_IT = [1 2; 3 4]
const B_IT = [2 1; 0 1]

run_and_report("U_trUAUdB_sym_itensor",
    () -> integrate_tr_uau_db_itensor(measure_sym; idx_dim = ITENSOR_POLY_IDX_DIM, avals = A_IT, bvals = B_IT))

run_and_report("U_tr(UAUdB)^2_sym_itensor",
    () -> integrate_tr_uau_db_sq_itensor(measure_sym; idx_dim = ITENSOR_POLY_IDX_DIM, avals = A_IT, bvals = B_IT))

# ============================================================================
# Save results
# ============================================================================
results["_meta"] = benchmark_meta()
output_path = joinpath(@__DIR__, "results_intu.json")

println("\n" * "="^72)
println("Summary (median times in ms)")
println("="^72)
@printf("%-30s %12s\n", "Benchmark", "Median (ms)")
println("-"^42)
for (name, data) in sort(collect(results), by = x->x[1])
    if name == "_meta"
        continue
    end
    @printf("%-30s %12.2f\n", name, data["median_ms"])
end

# Write JSON manually (avoids JSON3 dependency)
open(output_path, "w") do io
    println(io, "{")
    entries = sort(collect(results), by = x->x[1])
    for (idx, (name, data)) in enumerate(entries)
        comma = idx < length(entries) ? "," : ""
        if name == "_meta"
            println(io, "  \"_meta\": $(to_json(data))$comma")
        else
            ms = data["median_ms"]
            res = data["result"]
            n = data["samples"]
            status = get(data, "status", "ok")
            println(
                io,
                "  \"$name\": {\"median_ms\": $ms, \"result\": \"$(escape_string(string(res)))\", \"samples\": $n, \"status\": \"$(escape_string(string(status)))\"}$comma",
            )
        end
    end
    println(io, "}")
end
println("\nResults saved to $(output_path)")
