@testset "GOE Integration" begin
    # Test for N=3
    N = 3
    H = SymbolicMatrix(:H, :H, N)
    meas = dGOE(N)

    @testset "Tr(H)" begin
        expr = IntU.tr(H)
        res = integrate(expr, meas)
        @test to_numeric(res) == 0
    end

    @testset "Tr(H^2)" begin
        expr = IntU.tr(H^2)
        res = integrate(expr, meas)
        @test to_numeric(res) == N^2 + N
    end

    @testset "Tr(H^4)" begin
        expr = IntU.tr(H^4)
        res = integrate(expr, meas)
        expected = 2*N^3 + 5*N^2 + 5*N
        @test to_numeric(res) == expected
    end

    @testset "Component Moments" begin
        # < H_11^2 > = 2
        res1 = integrate(H[1, 1]^2, meas)
        @test to_numeric(res1) == 2

        # < H_12^2 > = 1
        res2 = integrate(H[1, 2]^2, meas)
        @test to_numeric(res2) == 1

        # < H_12 H_21 > = 1 (due to real symmetry)
        res3 = integrate(H[1, 2]*H[2, 1], meas)
        @test to_numeric(res3) == 1
    end
end
