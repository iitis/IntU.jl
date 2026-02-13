using IntU
using Test
using Symbolics
using LinearAlgebra

function is_really_zero(x)
    x = Symbolics.simplify(x)
    IntU._symbolic_isequal(x, 0) && return true
    x = Symbolics.expand(x)
    x = Symbolics.simplify(x)
    IntU._symbolic_isequal(x, 0) && return true
    
    vars = Symbolics.get_variables(x)
    if isempty(vars)
        v = Symbolics.value(x)
        return v isa Number && abs(v) < 1e-10
    end
    
    for i in 1:3
        subs = Dict(v => rand() + 0.1 for v in vars)
        try
            val_sub = Symbolics.substitute(x, subs)
            v = to_numeric(val_sub)
            if v isa Number
                if abs(v) < 1e-9 continue end
                println("DEBUG: is_really_zero sampling fail round $i: abs($v) = $(abs(v))")
                # println("  vars = ", vars)
                # println("  subs = ", subs)
                # println("  val_sub = ", val_sub)
                return false
            end
            println("DEBUG: is_really_zero sampling error round $i (not a number): $v")
            return false
        catch e
            println("DEBUG: is_really_zero substitution error round $i: $e")
            return false
        end
    end
    return true
end

@testset verbose=true "Pure States Integration" begin
    @variables d
    psi = SymbolicMatrix(:psi, :psi, d)

    @testset "Diagonal Term" begin
        res = integrate(psi[1, 1] * conj(psi[1, 1]), dPsi(d))
        @test is_really_zero(res - 1/d)
    end

    @testset "Off-Diagonal Term" begin
        res = integrate(psi[1, 1] * conj(psi[2, 1]), dPsi(d))
        @test is_really_zero(res)
    end

    @testset "Fidelity Average" begin
        # Use real components to avoid complex variable simplification issues
        @variables r[1:2] i[1:2]
        phi = [r[1] + im*i[1], r[2] + im*i[2]]
        
        inner_prod = conj(psi[1, 1])*phi[1] + conj(psi[2, 1])*phi[2]
        expr = inner_prod * conj(inner_prod)
        res = integrate(expr, dPsi(d))
        
        # Expected: sum_j |phi_j|^2 / d
        expected = (phi[1]*conj(phi[1]) + phi[2]*conj(phi[2])) / d
        @test is_really_zero(res - expected)
    end
end
