using IntU
using Test
using Symbolics
using LinearAlgebra

@testset verbose=true "Pure States Integration" begin
    @variables d
    @variables d
    @variables psi[1:3]::Complex

    function is_really_zero(x)
        # 1. Structural simplified zero
        x_simp = Symbolics.simplify(x)
        if IntU._symbolic_isequal(x_simp, 0) || IntU._symbolic_isequal(x_simp, Num(0))
            return true
        end

        # 2. Numerical evaluation with random substitution
        # Replace all variables with random numbers
        vars = Symbolics.get_variables(x_simp)
        subs = Dict(v => rand() + 0.1 for v in vars) # Avoid 0
        try
            val = Symbolics.substitute(x_simp, subs)
            # Unwrap potential Num wrapper
            val_un = Symbolics.unwrap(val)
            if val_un isa Number
                return abs(val_un) < 1e-12
            end
            # If still symbolic, force eval?
            return abs(eval(Meta.parse(string(val)))) < 1e-12
        catch
            return false
        end
    end

    @testset verbose=true "Diagonal Term" begin
        # 1. <psi_i | psi_i> average
        expr1 = psi[1] * conj(psi[1])
        res1 = integrate(expr1, dPsi(psi, d))
        expected1 = 1 / d
        # Use robust symbolic check
        diff = res1 - expected1
        @test is_really_zero(diff)
    end

    @testset verbose=true "Off-Diagonal Term" begin
        # 2. <psi_i | psi_j> average (i != j)
        expr2 = psi[1] * conj(psi[2])
        res2 = integrate(expr2, dPsi(psi, d))
        @test Symbolics.iszero(Symbolics.simplify(res2))
    end

    @testset verbose=true "Fidelity Average" begin
        # 3. |<psi|phi>|^2 average where phi is fixed
        # This should be 1/d
        @variables phi[1:3]::Complex
        # Use explicit sum to ensure proper scalar expansion
        inner_prod = sum(conj(psi[i]) * phi[i] for i = 1:3)
        expr3 = inner_prod * conj(inner_prod)

        res3 = integrate(expr3, dPsi(psi, d))

        # Expected: sum_i |phi_i|^2 / d
        expected3 = sum(phi[i]*conj(phi[i]) for i = 1:3) / d

        # Ensure scalar shape
        res3_scalar = res3
        if res3 isa AbstractArray
            res3_scalar = res3[1]
        elseif !isa(res3, Number) && applicable(length, res3) && length(res3) == 1
            res3_scalar = first(res3)
        end

        diff = res3_scalar - expected3
        @test is_really_zero(diff)
    end
end
