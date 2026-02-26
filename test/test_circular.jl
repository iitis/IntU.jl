using IntU
using Symbolics
using LinearAlgebra
using Test

@testset "Circular Ensembles" begin
    @variables d

    # helper for symbolic zero check
    function sym_iszero(x)
        # 1. Quick check
        if isequal(x, 0) || isequal(x, 0.0)
            return true
        end

        # 2. Simplify and check
        s = simplify(expand(x))
        if isequal(s, 0) || isequal(s, 0.0)
            return true
        end
        un = Symbolics.unwrap(s)
        if isequal(Num(un), 0)
            return true
        end

        # 3. Numerical fallback
        for d_val in [3.14, 1.23, 7.89]
            val = Symbolics.substitute(s, Dict(d => d_val))
            num_val = to_numeric(val)
            if num_val isa Number
                if abs(num_val) > 1e-12
                    return false
                end
            else
                return false
            end
        end
        return true
    end

    # 1. COE
    @testset "COE Properties" begin
        S = SymbolicMatrix(:S, :COE, d)
        m = dCOE(d)

        @test sym_iszero(integrate(S[1, 1], m))
        @test sym_iszero(integrate(S[1, 2], m))

        # E[S_11 S*_11] = 2/(d+1)
        val_11 = integrate(S[1, 1] * conj(S[1, 1]), m)
        @test sym_iszero(val_11 - 2 / (d + 1))

        # E[S_12 S*_12] = 1/(d+1)
        val_12 = integrate(S[1, 2] * conj(S[1, 2]), m)
        @test sym_iszero(val_12 - 1 / (d + 1))

        # Symmetry check
        val_sym = integrate(S[1, 2] * conj(S[2, 1]), m)
        @test sym_iszero(val_sym - val_12)
    end

    # 2. CSE
    @testset "CSE Properties" begin
        S_cse = SymbolicMatrix(:S, :CSE, d)
        mc = dCSE(d)

        @test sym_iszero(integrate(S_cse[1, 1], mc))

        # E[|S_11|^2] = 1/(d-1)
        val_11 = integrate(S_cse[1, 1] * conj(S_cse[1, 1]), mc)
        @test sym_iszero(val_11 - 1 / (d - 1))

        val_12 = integrate(S_cse[1, 2] * conj(S_cse[1, 2]), mc)
        @test sym_iszero(val_12)

        val_cross = integrate(S_cse[1, 2] * conj(S_cse[2, 1]), mc)
        @test sym_iszero(val_cross + val_12)
    end

    # 3. CUE
    @testset "CUE Properties" begin
        U = SymbolicMatrix(:U, :U, d)
        m = dCUE(d)

        @test sym_iszero(integrate(U[1, 1], m))
        val_11 = integrate(U[1, 1] * conj(U[1, 1]), m)
        @test sym_iszero(val_11 - 1/d)

        val_22 = integrate(abs(U[1, 1]*U[2, 2])^2, m)
        @test sym_iszero(val_22 - 1/(d^2-1))
    end

end
