"""
Performance comparison benchmarks for IntU.jl.
Computes the same integrals as bench_haarpy.py for a head-to-head comparison.

Usage:
    julia --project=/path/to/IntU.jl bench_intu.jl
"""

using IntU
using Symbolics
using BenchmarkTools
using Printf

@variables d

# Match Haarpy benchmark settings
BenchmarkTools.DEFAULT_PARAMETERS.samples = 30
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 120

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
    b = @benchmark $f()
    ms = median_ms(b)
    @printf(" %.2f ms  (result: %s)\n", ms, string(res))
    results[name] = Dict("median_ms" => ms, "result" => string(res), "samples" => length(b.times))
    return ms
end

# ============================================================================
# Section 1: Unitary |U_11|^{2k}, symbolic d
# ============================================================================
println("\n=== Unitary: |U_11|^{2k}, symbolic d ===")

U = SymbolicMatrix(:U, :U)
measure_sym = dU(d)

run_and_report("U_|U11|^6_sym", () -> integrate(abs(U[1,1])^6, measure_sym))
run_and_report("U_|U11|^8_sym", () -> integrate(abs(U[1,1])^8, measure_sym))
run_and_report("U_|U11|^10_sym", () -> integrate(abs(U[1,1])^10, measure_sym))

# ============================================================================
# Section 2: Unitary |U_11|^{10}, numeric d
# ============================================================================
println("\n=== Unitary: |U_11|^{2k}, numeric d ===")

for d_val in [10, 50]
    U_n = SymbolicMatrix(:U, :U, d_val)
    m_n = dU(d_val)
    run_and_report("U_|U11|^10_d=$d_val", () -> integrate(abs(U_n[1,1])^10, m_n))
end

# ============================================================================
# Section 3: Orthogonal O_11^k, symbolic d
# ============================================================================
println("\n=== Orthogonal: O_11^k, symbolic d ===")

O = SymbolicMatrix(:O, :O)
mO_sym = dO(d)

run_and_report("O_O11^2_sym", () -> integrate(O[1,1]^2, mO_sym))
run_and_report("O_O11^4_sym", () -> integrate(O[1,1]^4, mO_sym))

# ============================================================================
# Section 4: Orthogonal O_11^k, numeric d
# ============================================================================
println("\n=== Orthogonal: O_11^k, numeric d ===")

O10 = SymbolicMatrix(:O, :O, BigInt(10))
mO10 = dO(BigInt(10))
run_and_report("O_O11^6_d=10", () -> integrate(O10[1,1]^6, mO10))

O20 = SymbolicMatrix(:O, :O, BigInt(20))
mO20 = dO(BigInt(20))
run_and_report("O_O11^8_d=20", () -> integrate(O20[1,1]^8, mO20))
run_and_report("O_O11^10_d=20", () -> integrate(O20[1,1]^10, mO20))

O50 = SymbolicMatrix(:O, :O, BigInt(50))
mO50 = dO(BigInt(50))
run_and_report("O_O11^10_d=50", () -> integrate(O50[1,1]^10, mO50))

# ============================================================================
# Section 5: COE |S_11|^{2k}, symbolic d
# ============================================================================
println("\n=== Circular Orthogonal (COE): |S_11|^{2k}, symbolic d ===")

S_coe = SymbolicMatrix(:S, :COE)
mCOE = dCOE(d)

run_and_report("COE_|S11|^2_sym", () -> integrate(abs(S_coe[1,1])^2, mCOE))
run_and_report("COE_|S11|^4_sym", () -> integrate(abs(S_coe[1,1])^4, mCOE))
run_and_report("COE_|S11|^6_sym", () -> integrate(abs(S_coe[1,1])^6, mCOE))

# ============================================================================
# Section 6: Permutation P_11^k
# ============================================================================
println("\n=== Permutation: P_11^k ===")

P100 = SymbolicMatrix(:P, :P, 100)
mP100 = dPerm(100)
run_and_report("Perm_P11^10_d=100", () -> integrate(P100[1,1]^10, mP100))

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
