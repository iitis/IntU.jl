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

@variables d

# Match RTNI benchmark settings
BenchmarkTools.DEFAULT_PARAMETERS.samples = 30

function median_ms(b)
    m = median(b)
    return m.time / 1e6  # ns -> ms
end

results = Dict{String,Any}()

function run_and_report(name, f)
    print("  Running: $name ...")
    # Verify it works
    res = f()
    res = simplify(res)
    b = @benchmark $f() evals=1 samples=30 setup=(Memoization.empty_all_caches!())
    ms = median_ms(b)
    @printf(" %.2f ms  (result: %s)\n", ms, string(res))
    results[name] = Dict("median_ms" => ms, "result" => string(res), "samples" => length(b.times))
    return ms
end

# ============================================================================
# Section 1: Easy - Unitary |U_11|^{2k}, symbolic d (k=1,2,3)
# ============================================================================
println("\n=== Easy: Unitary |U_11|^{2k}, symbolic d ===")

U = SymbolicMatrix(:U, :U)
measure_sym = dU(d)

run_and_report("U_|U11|^2_sym", () -> integrate(abs(U[1,1])^2, measure_sym))
run_and_report("U_|U11|^4_sym", () -> integrate(abs(U[1,1])^4, measure_sym))
run_and_report("U_|U11|^6_sym", () -> integrate(abs(U[1,1])^6, measure_sym))

# ============================================================================
# Section 2: Harder - Unitary |U_11|^{2k}, symbolic d (k=4,5)
# ============================================================================
println("\n=== Harder: Unitary |U_11|^{2k}, symbolic d ===")

run_and_report("U_|U11|^8_sym", () -> integrate(abs(U[1,1])^8, measure_sym))
run_and_report("U_|U11|^10_sym", () -> integrate(abs(U[1,1])^10, measure_sym))

# ============================================================================
# Section 3: Harder - Unitary |U_11|^{10}, numeric d
# ============================================================================
println("\n=== Harder: Unitary |U_11|^{2k}, numeric d ===")

for d_val in [10, 50]
    U_n = SymbolicMatrix(:U, :U, d_val)
    m_n = dU(d_val)
    run_and_report("U_|U11|^10_d=$d_val", () -> integrate(abs(U_n[1,1])^10, m_n))
end

# ============================================================================
# Section 4: Trace moments (graph-based in RTNI)
# ============================================================================
println("\n=== Trace moments: |tr(U)|^{2k}, symbolic d ===")

run_and_report("U_|trU|^4_sym", () -> integrate(abs(tr(U))^4, measure_sym))
run_and_report("U_|trU|^6_sym", () -> integrate(abs(tr(U))^6, measure_sym))
run_and_report("U_|trU|^8_sym", () -> integrate(abs(tr(U))^8, measure_sym))

# ============================================================================
# Section 5: Trace polynomial - tr(U A U^* B)
# ============================================================================
println("\n=== Trace polynomials: symbolic d ===")

A = SymbolicMatrix(:A, :Constant)
B = SymbolicMatrix(:B, :Constant)

run_and_report("U_trUAUdB_sym",
    () -> integrate(tr(U * A * U' * B), measure_sym))

run_and_report("U_tr(UAUdB)^2_sym",
    () -> integrate(tr(U * A * U' * B * U * A * U' * B), measure_sym))

# ============================================================================
# Save results
# ============================================================================
println("\n" * "="^72)
println("Summary (median times in ms)")
println("="^72)
@printf("%-30s %12s\n", "Benchmark", "Median (ms)")
println("-"^42)
for (name, data) in sort(collect(results), by=x->x[1])
    @printf("%-30s %12.2f\n", name, data["median_ms"])
end

# Write JSON manually (avoids JSON3 dependency)
open("results_intu.json", "w") do io
    println(io, "{")
    entries = sort(collect(results), by=x->x[1])
    for (idx, (name, data)) in enumerate(entries)
        ms = data["median_ms"]
        res = data["result"]
        n = data["samples"]
        comma = idx < length(entries) ? "," : ""
        println(io, "  \"$name\": {\"median_ms\": $ms, \"result\": \"$(escape_string(string(res)))\", \"samples\": $n}$comma")
    end
    println(io, "}")
end
println("\nResults saved to results_intu.json")
