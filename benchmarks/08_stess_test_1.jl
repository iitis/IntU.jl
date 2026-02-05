# ============================================================
# IntU.jl Haar-integration test script (Unitary/Orthogonal/Symplectic)
#
# Run:
#   julia --project=benchmarks benchmarks/08_stess_test_1.jl
#
# Notes:
# - Keep integrands polynomial: use conj(z) explicitly instead of abs(z),
#   because abs(z) introduces sqrt(.) in many CASes.
# - The dimension 'd' is kept symbolic for the integration measure.
# - We iterate over concrete sizes N for the matrix variables to create valid Julia arrays.
# ============================================================
using Symbolics
using IntU
using LinearAlgebra
using BenchmarkTools
using Statistics

# --------------------------
# Helpers
# --------------------------
function equal_symbolic(got, expected; subs = Vector{Dict}())
    # Try direct symbolic simplification first
    try
        diff = Symbolics.simplify(got - expected)
        if diff == 0
            return true
        end
    catch
        # fall through to substitution checks
    end
    # Fallback: validate by substituting several integer values
    for s in subs
        try
            dg = Symbolics.simplify(Symbolics.substitute(got, s))
            de = Symbolics.simplify(Symbolics.substitute(expected, s))
            if Symbolics.simplify(dg - de) != 0
                return false
            end
        catch
            return false
        end
    end
    return !isempty(subs) # if we had substitutions and they all passed
end

function run_example(name, expr, μ, expected; subs = Vector{Dict}(), benchmark = false)
    # Warm-up/Correctness check
    got = IntU.integrate(expr, μ)

    # Simplify the result for display
    # We try simplify(expand(.)) as it is often stronger
    simplified_got = try
        Symbolics.simplify(Symbolics.expand(got))
    catch
        Symbolics.simplify(got)
    end

    println("== $name ==")
    # println("Integrand: $expr") 

    ok = equal_symbolic(got, expected; subs = subs)
    status = ok ? "PASS" : "FAIL"
    println("$status | Expect: $expected | Got: $simplified_got")

    if benchmark
        println("  -> Benchmarking...")
        # Benchmark with median time (excluding compilation)
        t = @benchmark IntU.integrate($expr, $μ) samples=30 evals=1
        med_time_ms = median(t).time / 1e6
        println("  -> Median Time: $(round(med_time_ms, digits=2)) ms")
    end
    println("")
    return got
end

# Define symbolic dimension globally
@variables d

# ============================================================
# 1) Unitary group U(d)
# ============================================================
function test_unitary()
    println("------------------------------------------------------------")
    println("Testing Unitary U(d) moments (degree <= 4) using Symbolic Dimension")
    println("------------------------------------------------------------")

    @variables d::Int
    @symbolic_dimension U[1:d, 1:d]
    μU = dU(U)
    # Using specific values for substitution checks
    subsU = [Dict(d => 3), Dict(d => 4), Dict(d => 7)]

    # Note: U is infinite lazy matrix, so we can access any index
    run_example("U1: ∫ 1 dU", 1, μU, 1; subs = subsU)
    run_example("U2: ∫ U₁₁ dU = 0", U[1, 1], μU, 0; subs = subsU)
    run_example("U3: ∫ |U₁₁|² dU = 1/d", U[1, 1]*conj(U[1, 1]), μU, 1/d; subs = subsU)

    run_example("U4: ∫ U₁₁ conj(U₁₂) dU = 0", U[1, 1]*conj(U[1, 2]), μU, 0; subs = subsU)
    run_example(
        "U5: ∫ |U₁₁|⁴ dU = 2/(d(d+1))",
        (U[1, 1]*conj(U[1, 1]))^2,
        μU,
        2/(d*(d+1));
        subs = subsU,
    )
    run_example(
        "U6: ∫ |U₁₁|²|U₁₂|² dU = 1/(d(d+1))",
        (U[1, 1]*conj(U[1, 1]))*(U[1, 2]*conj(U[1, 2])),
        μU,
        1/(d*(d+1));
        subs = subsU,
    )
    run_example(
        "U7: ∫ |U₁₁|²|U₂₁|² dU = 1/(d(d+1))",
        (U[1, 1]*conj(U[1, 1]))*(U[2, 1]*conj(U[2, 1])),
        μU,
        1/(d*(d+1));
        subs = subsU,
    )
    run_example(
        "U8: ∫ |U₁₁|²|U₂₂|² dU = 1/(d²-1)",
        (U[1, 1]*conj(U[1, 1]))*(U[2, 2]*conj(U[2, 2])),
        μU,
        1/(d^2 - 1);
        subs = subsU,
    )
    run_example(
        "U9: ∫ U₁₁U₂₂ conj(U₁₂)conj(U₂₁) dU = -1/(d(d²-1))",
        U[1, 1]*U[2, 2]*conj(U[1, 2])*conj(U[2, 1]),
        μU,
        -1/(d*(d^2 - 1));
        subs = subsU,
    )
end

# ============================================================
# 2) Orthogonal group O(d)
# ============================================================
function test_orthogonal(N::Int)
    println("------------------------------------------------------------")
    println("Testing Orthogonal O(d) moments (degree <= 4) using $N x $N matrix")
    println("------------------------------------------------------------")

    @variables O[1:N, 1:N]::Real
    μO = dO(O, d)
    subsO = [Dict(d => 3), Dict(d => 4), Dict(d => 8)]

    run_example("O1: ∫ O₁₁ dO = 0", O[1, 1], μO, 0; subs = subsO)
    run_example("O2: ∫ O₁₁² dO = 1/d", O[1, 1]^2, μO, 1/d; subs = subsO)
    run_example("O3: ∫ O₁₁⁴ dO = 3/(d(d+2))", O[1, 1]^4, μO, 3/(d*(d+2)); subs = subsO)

    if N >= 2
        run_example(
            "O4: ∫ O₁₁² O₁₂² dO = 1/(d(d+2))",
            O[1, 1]^2 * O[1, 2]^2,
            μO,
            1/(d*(d+2));
            subs = subsO,
        )
        run_example(
            "O5: ∫ O₁₁² O₂₂² dO = (d+1)/(d(d-1)(d+2))",
            O[1, 1]^2 * O[2, 2]^2,
            μO,
            (d+1)/(d*(d-1)*(d+2));
            subs = subsO,
        )
        run_example(
            "O6: ∫ O₁₁ O₁₂ O₂₁ O₂₂ dO = -1/(d(d-1)(d+2))",
            O[1, 1]*O[1, 2]*O[2, 1]*O[2, 2],
            μO,
            -1/(d*(d-1)*(d+2));
            subs = subsO,
        )
    end
end

# ============================================================
# 3) Symplectic group Sp(d)
# ============================================================
function test_symplectic(N::Int)
    if isodd(N) || N < 4
        # Skip N=2 (Sp(2)) because the standard Weingarten matrix for degree 6 
        # is singular at d=2 (Sp(2) uses O(-d) relation).
        # degree 4 (Sp2) should work, but for consistency in this example 
        # we start from N=4.
        return
    end
    println("------------------------------------------------------------")
    println("Testing Symplectic Sp(d) using $N x $N matrix (d=$N)")
    println("------------------------------------------------------------")

    @variables S[1:N, 1:N]::Complex
    μSp = dSp(S, N)
    subsSp = [Dict(d => N)]

    run_example("Sp1: ∫ |S₁₁|² dSp = 1/d", S[1, 1]*conj(S[1, 1]), μSp, 1/d; subs = subsSp)
    # For Sp(d), first row distribution matches U(d), so moments match.
    # E|S11|^4 = 2 / (d(d+1))
    run_example(
        "Sp2: ∫ |S₁₁|⁴ dSp = 2/(d(d+1))",
        (S[1, 1]*conj(S[1, 1]))^2,
        μSp,
        2/(d*(d+1));
        subs = subsSp,
    )
end

# ============================================================
# 4) High-degree Moments (Run once to save time)
# ============================================================
function test_high_moments()
    println("------------------------------------------------------------")
    println("Testing High-degree Moments (Degree 6) and Benchmarking")
    println("------------------------------------------------------------")

    # Unitary U(d) 6-th moment
    @variables d::Int
    @symbolic_dimension U[1:d, 1:d]
    μU = dU(U)
    run_example(
        "U10: ∫ |U₁₁|⁶ dU = 6/(d(d+1)(d+2)) [Symbolic d]",
        (U[1, 1]*conj(U[1, 1]))^3,
        μU,
        6/(d*(d+1)*(d+2));
        subs = [Dict(d => 3)],
        benchmark = true,
    )

    # Orthogonal O(d) 6-th moment - slow with symbolic inversion, using concrete dimension for demonstration
    @variables O[1:1, 1:1]::Real
    println("Note: Orthogonal 6th moment is computed with concrete d=10 for speed.")
    μO_concrete = dO(O, 10)
    run_example(
        "O7: ∫ O₁₁⁶ dO = 15/(d(d+2)(d+4)) [Concrete d=10]",
        O[1, 1]^6,
        μO_concrete,
        15//(10*(10+2)*(10+4)),
        benchmark = true,
    )
end

# ============================================================
# 5) Application-style example
# ============================================================
function test_application()
    println("------------------------------------------------------------")
    println("Application: bipartite random pure state (via Haar U(D))")
    println("------------------------------------------------------------")
    nA, nB = 2, 3
    D = nA*nB
    @variables Uψ[1:D, 1:D]::Complex
    μConcrete = dU(Uψ, D)
    psi(a, b) = Uψ[(a-1)*nB+b, 1]
    purity = zero(Num)
    for a = 1:nA, ap = 1:nA, b = 1:nB, bp = 1:nB
        purity += psi(a, b) * conj(psi(ap, b)) * psi(ap, bp) * conj(psi(a, bp))
    end
    expected_purity = (nA + nB) / (D + 1)
    run_example(
        "Purity (nA=2,nB=3,D=6): E[tr(ρ_A^2)]",
        purity,
        μConcrete,
        expected_purity,
        benchmark = true,
    )
end

# ============================================================
# Main Loop
# ============================================================
println("============================================================")
println("Running IntU examples with symbolic dimension 'd'")
println("============================================================\n")

# Run unitary once (no benchmark)
test_unitary()

# Run loops for low-degree moments (no benchmark)
for N = 2:4
    test_orthogonal(N)
    test_symplectic(N)
end

# Run high-degree moments ONCE (WITH BENCHMARK)
test_high_moments()

# Run application (WITH BENCHMARK)
test_application()

println("All examples completed.")
