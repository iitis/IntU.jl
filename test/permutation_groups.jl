using IntU
using Test
using Symbolics

@testset "Permutation Groups" begin
    @variables d
    P = [Symbolics.variable(:P, i, j) for i = 1:10, j = 1:10]
    measure = dPerm(P, d)

    function is_zero(x)
        # First try symbolic simplification
        simplified = Symbolics.simplify(Symbolics.expand(x))
        v = Symbolics.value(simplified)
        if v isa Number
            return iszero(v)
        end
        # For symbolic d, substitute with a large numeric value and check
        subs_val = Symbolics.substitute(simplified, d => 1000)
        num_val = Symbolics.value(subs_val)
        if num_val isa Number
            return abs(num_val) < 1e-10
        end
        return isequal(simplified, 0)
    end

    @testset "Basic integration" begin
        # E[P_11] = 1/d
        @test is_zero(integrate(P[1, 1], measure) - 1/d)

        # E[P_11 * P_22] = 1 / (d(d-1))
        @test is_zero(integrate(P[1, 1] * P[2, 2], measure) - 1/(d*(d-1)))

        # E[P_11 * P_12] = 0 (same row, different column)
        @test is_zero(integrate(P[1, 1] * P[1, 2], measure))

        # E[P_11 * P_21] = 0 (different row, same column)
        @test is_zero(integrate(P[1, 1] * P[2, 1], measure))

        # E[P_11^2] = E[P_11] = 1/d (since P_ij in {0,1})
        @test is_zero(integrate(P[1, 1]^2, measure) - 1/d)
    end

    @testset "Higher order" begin
        # E[P_11 * P_22 * P_33] = 1 / (d(d-1)(d-2))
        @test is_zero(integrate(P[1, 1] * P[2, 2] * P[3, 3], measure) - 1/(d*(d-1)*(d-2)))
    end

    @testset "Centered Permutations" begin
        Y = [Symbolics.variable(:Y, i, j) for i = 1:10, j = 1:10]
        m_centered = dCPerm(Y, d)

        # E[Y_11] = 0
        @test is_zero(integrate(Y[1, 1], m_centered))

        # E[Y_11 * Y_11] = (d-1)/d^2
        res2 = integrate(Y[1, 1]^2, m_centered)
        @test is_zero(res2 - (d-1)/d^2)

        # E[Y_11 * Y_12]
        res12 = integrate(Y[1, 1] * Y[1, 2], m_centered)
        @test is_zero(res12 - (-1/d^2))
    end

    @testset "Numeric dimension" begin
        P_num = [Symbolics.variable(:P, i, j) for i = 1:3, j = 1:3]
        m_num = dPerm(P_num, 3)

        @test integrate(P_num[1, 1], m_num) == 1//3
        @test integrate(P_num[1, 1] * P_num[2, 2], m_num) == 1//6
        @test integrate(P_num[1, 1] * P_num[1, 2], m_num) == 0
    end
end
