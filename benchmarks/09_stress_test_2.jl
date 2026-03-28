#!/usr/bin/env julia
# Usage:
#   julia benchmarks/09_stress_test_2.jl --out results.json --samples 50
#
# What it does:
#   - Runs a battery of higher-degree Haar integrals (8th/10th moments, mixed moments)
#     over U(d), O(d), and Sp(d) using IntU.jl.
#   - Compares results against closed-form expected expressions.
#   - Benchmarks runtime (median time, memory, allocs) per test case.
#   - Writes a JSON report for regression testing.

using IntU
using Symbolics
using BenchmarkTools
using JSON3
using Dates

# --- CLI args ----------------------------------------------------------------
function get_arg(flag::String, default::String)
    for (i, a) in pairs(ARGS)
        if a == flag && i < length(ARGS)
            return ARGS[i+1]
        end
    end
    return default
end

function has_flag(flag::String)
    any(==(flag), ARGS)
end

outpath = get_arg("--out", "intu_stress_results.json")
samples = parse(Int, get_arg("--samples", "30"))
quick = has_flag("--quick")

doublefactorial_odd(n::Int) = prod(1:2:n)

function sym_is_zero(x)
    return IntU._iszero(x)
end

function safe_string(x)
    return string(x)
end

function eval_at(expr, dvar, val::Int)
    subbed = Symbolics.substitute(expr, Dict(dvar => val))
    return Symbolics.simplify(subbed)
end

function full_simplify(x)
    return Symbolics.simplify(Symbolics.expand(x))
end

# --- Expected-value constructors ----------------------------------------------
# U(d): if (|U_11|^2) = X, then X ~ Beta(1, d-1), so E[X^k] = k! / (d)(d+1)...(d+k-1)
unitary_abs2_moment(d, k::Int) = factorial(k) / prod(d + i for i = 0:(k-1))

# U(d): for distinct components in a row, (|U_1j|^2) are Dirichlet(1,...,1)
# E[∏_{j=1}^m |U_1j|^2] = 1 / (d)(d+1)...(d+m-1)
unitary_row_product(d, m::Int) = 1 / prod(d + i for i = 0:(m-1))

# U(d): E[(|U_11|^2)^2 (|U_12|^2)^2] under Dirichlet(1,...,1)
unitary_row_square_square(d) = (2*2) / prod(d + i for i = 0:3)  # 4 / d(d+1)(d+2)(d+3)

# O(d): if O_11^2 ~ Beta(1/2, (d-1)/2), then E[O_11^(2k)] = (2k-1)!! / d(d+2)...(d+2k-2)
orth_entry_even_moment(d, k::Int) = doublefactorial_odd(2k-1) / prod(d + 2i for i = 0:(k-1))

# O(d): squared coordinates Dirichlet(1/2,...,1/2)
# E[∏_{j=1}^m O_1j^2] = 1 / d(d+2)...(d+2m-2)
orth_row_product(d, m::Int) = 1 / prod(d + 2i for i = 0:(m-1))

# O(d): E[(O_11^2)^2 (O_12^2)^2] with Dirichlet(1/2,...,1/2) gives 9 / d(d+2)(d+4)(d+6)
orth_row_square_square(d) = 9 / prod(d + 2i for i = 0:3)

# Sp(d): we use the same |S_11|^2-moment formulas as U(d) (complex entry magnitude-squared moments),
symp_abs2_moment(d, k::Int) = unitary_abs2_moment(d, k)

function bench_integrate(expr, μ; samples::Int)
    r = integrate(expr, μ)
    t = @benchmark integrate($expr, $μ) evals=1 samples=samples
    med = median(t)
    return r,
    Dict(
        "median_ms" => med.time / 1e6,
        "memory_bytes" => med.memory,
        "allocs" => med.allocs,
        "samples" => length(t.times),
    )
end

results = Dict{String,Any}(
    "meta" => Dict(
        "timestamp" => string(Dates.now()),
        "julia_version" => string(VERSION),
        "samples" => samples,
        "quick" => quick,
    ),
    "cases" => Any[],
)

function push_case!(case_result)
    push!(results["cases"], case_result)
end

dvals_symbolic = quick ? [4, 8, 16] : [4, 8, 16, 32, 64]

println("=== IntU.jl stress/benchmark suite ===")
println("samples per benchmark: $samples")
println("JSON output: $outpath")
println()

# ---------------- U(d): symbolic-d stress tests -------------------------------
@variables d::Int
U = SymbolicMatrix(:U, :U, d)
μU = dU(d)

U_cases = [
    ("U_abs2_pow4__|U11|^8", abs2(U[1, 1])^4, unitary_abs2_moment(d, 4)),
    ("U_abs2_pow5__|U11|^10", abs2(U[1, 1])^5, unitary_abs2_moment(d, 5)),
    ("U_row_prod4__∏|U1j|^2", prod(abs2(U[1, j]) for j = 1:4), unitary_row_product(d, 4)),
    (
        "U_row_sq_sq__|U11|^4|U12|^4",
        (abs2(U[1, 1])^2) * (abs2(U[1, 2])^2),
        unitary_row_square_square(d),
    ),
]

function safe_eq(x, y)
    val = (x == y)
    return val isa Bool ? val : false
end

for (name, expr, expected) in U_cases
    got, bm = bench_integrate(expr, μU; samples = samples)
    diff0 = sym_is_zero(got - expected)

    numeric_checks = Dict{String,Any}()
    for dv in dvals_symbolic
        g = eval_at(got, d, dv)
        e = eval_at(expected, d, dv)
        numeric_checks[string(dv)] = Dict(
            "got" => safe_string(g),
            "expected" => safe_string(e),
            "ok" => safe_eq(g, e),
        )
    end

    println("[U(d)] $name")
    println("  integrand : $(safe_string(expr))")
    println("  expected  : $(safe_string(full_simplify(expected)))")
    println("  got       : $(safe_string(full_simplify(got)))")
    println("  symbolic ok: $diff0")
    println(
        "  median: $(bm["median_ms"]) ms, mem: $(bm["memory_bytes"]) B, allocs: $(bm["allocs"])",
    )
    println()

    push_case!(
        Dict(
            "group" => "U",
            "name" => name,
            "d" => "symbolic",
            "integrand" => safe_string(expr),
            "expected" => safe_string(expected),
            "got" => safe_string(got),
            "symbolic_ok" => diff0,
            "numeric_checks" => numeric_checks,
            "benchmark" => bm,
        ),
    )
end

# ------------- O(d): symbolic-d stress tests ----------------------------------
N_sym = 10
@variables d
O_sym = SymbolicMatrix(:O, :O, d)
O = O_sym[1:N_sym, 1:N_sym]
μO = dO(d)

O_cases = [
    # (name, expr, k, degree_for_expected)
    ("O_pow8__O11^8", O[1, 1]^8, 4),
    ("O_pow10__O11^10", O[1, 1]^10, 5),
    ("O_row_prod4__∏O1j^2", prod(O[1, j]^2 for j = 1:4), 4),
    ("O_row_sq_sq__O11^4 O12^4", (O[1, 1]^4) * (O[1, 2]^4), 4),
]

for (name, expr, k) in O_cases
    use_concrete_O = (k >= 4)
    local_μO = use_concrete_O ? dO(20) : μO
    local_d = use_concrete_O ? 20 : d

    local_expected = if name == "O_row_prod4__∏O1j^2"
        orth_row_product(local_d, 4)
    elseif name == "O_row_sq_sq__O11^4 O12^4"
        orth_row_square_square(local_d)
    else
        orth_entry_even_moment(local_d, k)
    end

    println("[O(d)] $name $(use_concrete_O ? "(concrete d=20)" : "(symbolic)")")
    got, bm = bench_integrate(expr, local_μO; samples = samples)

    ok = if use_concrete_O
        isapprox(Symbolics.value(got), Symbolics.value(local_expected); atol = 1e-13)
    else
        sym_is_zero(got - local_expected)
    end

    numeric_checks = Dict{String,Any}()
    if !use_concrete_O
        for dv in dvals_symbolic
            g = eval_at(got, d, dv)
            e = eval_at(local_expected, d, dv)
            numeric_checks[string(dv)] = Dict(
                "got" => safe_string(g),
                "expected" => safe_string(e),
                "ok" => safe_eq(g, e),
            )
        end
    end

    println("  integrand : $(safe_string(expr))")
    println("  expected  : $(safe_string(full_simplify(local_expected)))")
    println("  got       : $(safe_string(full_simplify(got)))")
    println("  numeric ok: $ok")
    println(
        "  median: $(bm["median_ms"]) ms, mem: $(bm["memory_bytes"]) B, allocs: $(bm["allocs"])",
    )
    println()

    push_case!(
        Dict(
            "group" => "O",
            "name" => name,
            "d" => use_concrete_O ? 20 : "symbolic",
            "integrand" => safe_string(expr),
            "expected" => safe_string(local_expected),
            "got" => safe_string(got),
            "symbolic_ok" => ok,
            "numeric_checks" => numeric_checks,
            "benchmark" => bm,
        ),
    )
end

# ------------- "Bigger d" sanity checks (explicit numeric matrices) -----------
dnums = quick ? [10, 20] : [10, 20, 50]

for dnum in dnums
    if quick && dnum > 10
        continue
    end
    Ubig_sym = SymbolicMatrix(:Ubig, :U, dnum)
    Ubig = Ubig_sym[1:dnum, 1:dnum]
    # Keep concrete dimensions integer-like for strict dimension validation.
    μUbig = dU(dnum)
    expr = abs2(Ubig[1, 1])^5
    expected = 120.0 / (dnum*(dnum+1)*(dnum+2)*(dnum+3)*(dnum+4))

    got, bm = bench_integrate(expr, μUbig; samples = samples)
    ok = isapprox(Symbolics.value(got), expected; atol = 1e-13)

    println("[U(d) numeric] d=$dnum  |U11|^10")
    println("  expected: $expected")
    println("  got     : $(safe_string(got))")
    println("  ok      : $ok")
    println(
        "  median: $(bm["median_ms"]) ms, mem: $(bm["memory_bytes"]) B, allocs: $(bm["allocs"])",
    )
    println()

    push_case!(
        Dict(
            "group" => "U",
            "name" => "U_numeric_d_$(dnum)__|U11|^10",
            "d" => dnum,
            "integrand" => safe_string(expr),
            "expected" => safe_string(expected),
            "got" => safe_string(got),
            "symbolic_ok" => ok,
            "benchmark" => bm,
        ),
    )
end

for dnum in dnums
    if quick && dnum > 10
        continue
    end
    Obig_sym = SymbolicMatrix(:Obig, :O, dnum)
    Obig = Obig_sym[1:dnum, 1:dnum]
    # Keep concrete dimensions integer-like for strict dimension validation.
    μObig = dO(dnum)
    expr = Obig[1, 1]^10
    expected =
        Float64((doublefactorial_odd(9)) // (dnum*(dnum+2)*(dnum+4)*(dnum+6)*(dnum+8)))  # 945/...

    got, bm = bench_integrate(expr, μObig; samples = samples)
    ok = isapprox(Symbolics.value(got), expected; atol = 1e-13)

    println("[O(d) numeric] d=$(dnum)  O11^10")
    println("  expected: $expected")
    println("  got     : $(safe_string(got))")
    println("  ok      : $ok")
    println(
        "  median: $(bm["median_ms"]) ms, mem: $(bm["memory_bytes"]) B, allocs: $(bm["allocs"])",
    )
    println()

    push_case!(
        Dict(
            "group" => "O",
            "name" => "O_numeric_d_$(dnum)__O11^10",
            "d" => dnum,
            "integrand" => safe_string(expr),
            "expected" => safe_string(expected),
            "got" => safe_string(got),
            "symbolic_ok" => ok,
            "benchmark" => bm,
        ),
    )
end

# ------------- Sp(d): stress tests -------------------------------------------
# (a) High moments of |S11| (handled by IntU.jl for Sp via dSp) :contentReference[oaicite:6]{index=6}
# (b) A “no conjugates” symplectic example with known closed form:
#     ∫ s_{1,1} s_{2,N+2} s_{N+1,2} s_{N+2,N+1} dSp
#     = 1 / (4 N (N-1) (2N+1)) for Sp(N) in 2N×2N complex form :contentReference[oaicite:7]{index=7}

Nvals = quick ? [3, 5] : [3, 5, 10]

for N in Nvals
    dSp_num = 2N
    S_sym = SymbolicMatrix(:S, :Sp, dSp_num)
    S = S_sym[1:dSp_num, 1:dSp_num]
    μSp = dSp(dSp_num)

    # (a) |S11|^8 and |S11|^10
    expr8 = abs2(S[1, 1])^4
    expr10 = abs2(S[1, 1])^5
    expected8 = 24 // (dSp_num*(dSp_num+1)*(dSp_num+2)*(dSp_num+3))
    expected10 = 120 // (dSp_num*(dSp_num+1)*(dSp_num+2)*(dSp_num+3)*(dSp_num+4))

    cases_Sp = Any[]
    # (a) |S11|^8 and |S11|^10
    # Weingarten formula for Sp(2N) is singular for degree 2k if N < k.
    if N >= 4
        push!(cases_Sp, ("Sp_numeric_d_$(dSp_num)__|S11|^8", expr8, expected8))
    end
    if N >= 5 && (!quick || N == 5)
        push!(cases_Sp, ("Sp_numeric_d_$(dSp_num)__|S11|^10", expr10, expected10))
    end

    for (nm, expr, expected) in cases_Sp
        # Keep concrete dimensions integer-like for strict dimension validation.
        local_μSp = dSp(dSp_num)
        local_expected = Float64(expected)

        got, bm = bench_integrate(expr, local_μSp; samples = samples)

        ok = isapprox(Symbolics.value(got), local_expected; atol = 1e-12)

        println("[Sp(d) numeric] N=$N (d=$dSp_num)  $nm")
        println("  expected: $local_expected")
        println("  got     : $(safe_string(got))")
        println("  ok      : $ok")
        println(
            "  median: $(bm["median_ms"]) ms, mem: $(bm["memory_bytes"]) B, allocs: $(bm["allocs"])",
        )
        println()

        push_case!(
            Dict(
                "group" => "Sp",
                "name" => nm,
                "d" => dSp_num,
                "integrand" => safe_string(expr),
                "expected" => safe_string(local_expected),
                "got" => safe_string(got),
                "symbolic_ok" => ok,
                "benchmark" => bm,
            ),
        )
    end

    # (b) Collins Example 4.7 style monomial (no conjugates)
    expr_collins = S[1, 1] * S[2, N+2] * S[N+1, 2] * S[N+2, N+1]
    expected_collins = Float64(1 // (4*N*(N-1)*(2N+1)))  # = 1 / (d(d-2)(d+1)) with d=2N

    μSp_f = dSp(dSp_num)
    gotC, bmC = bench_integrate(expr_collins, μSp_f; samples = samples)
    okC = isapprox(Symbolics.value(gotC), expected_collins; atol = 1e-14)

    println("[Sp(d) numeric] N=$N (d=$dSp_num)  Collins-example monomial")
    println("  integrand: $(safe_string(expr_collins))")
    println("  expected : $expected_collins")
    println("  got      : $(safe_string(gotC))")
    println("  ok       : $okC")
    println(
        "  median: $(bmC["median_ms"]) ms, mem: $(bmC["memory_bytes"]) B, allocs: $(bmC["allocs"])",
    )
    println()

    push_case!(
        Dict(
            "group" => "Sp",
            "name" => "Sp_numeric_N_$(N)__collins_example_degree4",
            "d" => dSp_num,
            "integrand" => safe_string(expr_collins),
            "expected" => safe_string(expected_collins),
            "got" => safe_string(gotC),
            "symbolic_ok" => okC,
            "benchmark" => bmC,
        ),
    )
end

open(outpath, "w") do io
    JSON3.write(io, results; indent = 2)
end

println("=== Done. Wrote JSON report to: $outpath ===")
