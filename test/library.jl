using IntegrateUnitary
using Test
using Symbolics

@testset "Integral Library" begin
    @variables d
    sym_eq(a, b) = IntegrateUnitary._symbolic_isequal(Symbolics.simplify(Symbolics.expand(a - b)), 0)

    @testset "Haar Unitary Trace" begin
        U = SymbolicMatrix(:U, :U)
        A = SymbolicMatrix(:A)
        B = SymbolicMatrix(:B)

        expr = IntegrateUnitary.tr(U * A * U' * B)
        res = integrate(expr, dU(d))

        # Expected: (tr(A) * tr(B)) / d
        expected = (IntegrateUnitary.tr_val([A]) * IntegrateUnitary.tr_val([B])) / d
        @test isequal(res, expected)
    end

    @testset "GUE Moments" begin
        H = SymbolicMatrix(:H, :GUE)

        @test isequal(integrate(IntegrateUnitary.tr(H^2), dGUE(d)), d^2)
        @test isequal(integrate(IntegrateUnitary.tr(H^4), dGUE(d)), 2d^3 + d)
        @test isequal(integrate(IntegrateUnitary.tr(H^6), dGUE(d)), 5d^4 + 10d^2)
    end

    @testset "GOE Moments" begin
        H = SymbolicMatrix(:H, :GOE)
        @test isequal(integrate(IntegrateUnitary.tr(H^2), dGOE(d)), d^2 + d)
        @test isequal(integrate(IntegrateUnitary.tr(H^4), dGOE(d)), 2d^3 + 5d^2 + 5d)
    end

    @testset "GSE Moments" begin
        H = SymbolicMatrix(:H, :GSE)
        @test isequal(integrate(IntegrateUnitary.tr(H^2), dGSE(d)), d^2 - d)
        @test isequal(integrate(IntegrateUnitary.tr(H^4), dGSE(d)), 2d^3 - 5d^2 + 5d)
    end

    @testset "Fallback Check" begin
        # Something not in library
        U_sym = SymbolicMatrix(:U, :U)
        expr = U_sym[1, 1] * conj(U_sym[1, 1])
        res = integrate(expr, dU(d))
        @test isequal(res, 1/d)
    end

    @testset "Prefactor Handling" begin
        # Verify that check_gaussian_library correctly applies the prefactor
        H = SymbolicMatrix(:H, :GUE)
        # 3 * tr(H^2) should yield 3 * d^2, not d^2
        expr = IntegrateUnitary.LazyTrace(Vector{AbstractMatrix}[[H, H]], Num(3))
        res = IntegrateUnitary.check_gaussian_library(expr, dGUE(d), :GUE)
        @test res !== nothing
        @test isequal(res, 3 * d^2)

        # tr(H^2) with prefactor=1 should yield d^2
        expr1 = IntegrateUnitary.LazyTrace(Vector{AbstractMatrix}[[H, H]], Num(1))
        res1 = IntegrateUnitary.check_gaussian_library(expr1, dGUE(d), :GUE)
        @test isequal(res1, d^2)
    end

    @testset "Gaussian Element-Wise Second Moments" begin
        H_gue = SymbolicMatrix(:H, :GUE)
        H_goe = SymbolicMatrix(:H, :GOE)
        H_gse = SymbolicMatrix(:H, :GSE)

        cases = [
            (H_gue[1, 1]^2, dGUE(d), 1),
            (H_gue[1, 1] * conj(H_gue[1, 1]), dGUE(d), 1),
            (H_goe[1, 1]^2, dGOE(d), 2),
            (H_goe[1, 2]^2, dGOE(d), 1),
            (H_gse[1, 1]^2, dGSE(d), 1),
        ]

        for (expr, measure, expected) in cases
            lib = IntegrateUnitary.check_library(expr, measure)
            @test lib !== nothing
            @test sym_eq(lib, expected)
            @test sym_eq(integrate(expr, measure), expected)
            @test sym_eq(IntegrateUnitary.fallback_integrate(expr, measure), expected)
        end

        @test IntegrateUnitary.check_library(H_gue[1, 1]^4, dGUE(d)) === nothing
    end

    @testset "Ginibre Library Entries" begin
        G_ue = SymbolicMatrix(:G, :GinUE)
        G_oe = SymbolicMatrix(:G, :GinOE)
        G_se = SymbolicMatrix(:G, :GinSE)

        ginue_cases = [
            (IntegrateUnitary.tr(G_ue * G_ue'), dGinUE(d), d^2),
            (IntegrateUnitary.tr((G_ue * G_ue')^2), dGinUE(d), 2d^3),
            (IntegrateUnitary.tr(G_ue * G_ue')^2, dGinUE(d), d^4 + d^2),
        ]

        for (expr, measure, expected) in ginue_cases
            lib = IntegrateUnitary.check_library(expr, measure)
            @test lib !== nothing
            @test sym_eq(lib, expected)
            @test sym_eq(integrate(expr, measure), expected)
            @test sym_eq(IntegrateUnitary.fallback_integrate(expr, measure), expected)
        end

        expr_goe = IntegrateUnitary.tr(G_oe * transpose(G_oe))
        lib_goe = IntegrateUnitary.check_library(expr_goe, dGinOE(d))
        @test lib_goe !== nothing
        @test sym_eq(lib_goe, d^2)
        @test sym_eq(integrate(expr_goe, dGinOE(d)), d^2)

        expr_gse = IntegrateUnitary.tr(G_se * G_se')
        lib_gse = IntegrateUnitary.check_library(expr_gse, dGinSE(d))
        @test lib_gse !== nothing
        @test sym_eq(lib_gse, d^2)
        @test sym_eq(integrate(expr_gse, dGinSE(d)), d^2)

        @test IntegrateUnitary.check_library(IntegrateUnitary.tr(G_ue^2), dGinUE(d)) === nothing
    end

    @testset "Orthogonal/Symplectic/Circular Library Entries" begin
        O = SymbolicMatrix(:O, :O)
        Sp = SymbolicMatrix(:Sp, :Sp)
        S_coe = SymbolicMatrix(:S, :COE)
        S_cse = SymbolicMatrix(:S, :CSE)

        low_order_cases = [
            (O[1, 1]^2, dO(d), 1 / d),
            (O[1, 1]^4, dO(d), 3 / (d * (d + 2))),
            (abs(Sp[1, 1])^2, dSp(d), 1 / d),
            (abs(Sp[1, 1])^4, dSp(d), 2 / (d + d^2)),
            (abs(Sp[1, 1])^2 * abs(Sp[1, 2])^2, dSp(d), 1 / (d + d^2)),
            (abs(S_coe[1, 1])^2, dCOE(d), 2 / (d + 1)),
            (abs(S_coe[1, 2])^2, dCOE(d), 1 / (d + 1)),
            (abs(S_coe[1, 1])^4, dCOE(d), 8 / ((d + 1) * (d + 3))),
            (abs(S_coe[1, 2])^4, dCOE(d), 2 / (d * (d + 3))),
            (
                abs(S_coe[1, 1])^2 * abs(S_coe[1, 2])^2,
                dCOE(d),
                2 / ((d + 1) * (d + 3)),
            ),
            (abs(S_cse[1, 1])^2, dCSE(d), 1 / (d - 1)),
            (abs(S_cse[1, 1])^4, dCSE(d), 2 / (-d + d^2)),
        ]

        for (expr, measure, expected) in low_order_cases
            lib = IntegrateUnitary.check_library(expr, measure)
            @test lib !== nothing
            @test sym_eq(lib, expected)
            @test sym_eq(integrate(expr, measure), expected)
            @test sym_eq(IntegrateUnitary.fallback_integrate(expr, measure), expected)
        end

        high_order_cases = [
            O[1, 1]^2 * O[1, 2]^4 * O[1, 3]^6,
            O[1, 1]^2 * O[2, 2]^4 * O[1, 3]^6,
        ]
        for expr in high_order_cases
            @test IntegrateUnitary.check_library(expr, dO(d)) !== nothing
        end
        @test IntegrateUnitary.check_library(abs(Sp[1, 1])^2 * abs(Sp[1, 2])^4 * abs(Sp[1, 3])^6, dSp(d)) !==
              nothing

        @test IntegrateUnitary.check_library(O[1, 1]^6, dO(d)) === nothing
        @test IntegrateUnitary.check_library(abs(S_coe[1, 1])^6, dCOE(d)) === nothing
    end
end
