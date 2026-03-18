using IntU
using Test

@testset verbose=true "Weingarten Calculus" begin
    @testset verbose=true "Weingarten Function Values" begin
        # Wg(1^2, d) = 1/(d^2-1)
        @test IntU.weingarten([1, 1], 3) == 1//8

        # Wg(2, d) = -1/(d(d^2-1))
        @test IntU.weingarten([2], 3) == -1//24
    end

    @testset verbose=true "Weingarten Unit Tests" begin
        @testset verbose=true "conjugate_partition" begin
            @test IntU.conjugate_partition([1]) == [1]
            @test IntU.conjugate_partition([2]) == [1, 1]
            @test IntU.conjugate_partition([1, 1]) == [2]
            @test IntU.conjugate_partition([2, 1]) == [2, 1]
            @test IntU.conjugate_partition([3, 1]) == [2, 1, 1]
            @test IntU.conjugate_partition([4]) == [1, 1, 1, 1]
            @test IntU.conjugate_partition(Int[]) == Int[]
        end

        @testset verbose=true "irrep_dimension (Dim of U(d) irrep)" begin
            # s_{1}(1^d) = d
            @test IntU.irrep_dimension([1], 3) == 3//1
            # s_{2}(1^d) = d(d+1)/2 symmetric tensor
            @test IntU.irrep_dimension([2], 3) == 3*4//2 # 6
            # s_{1,1}(1^d) = d(d-1)/2 antisymmetric
            @test IntU.irrep_dimension([1, 1], 3) == 3*2//2 # 3
        end

        @testset verbose=true "character_at_id (Dim of S_n irrep)" begin
            # S3
            # [3] -> 1 (trivial)
            @test IntU.character_at_id([3]) == 1
            # [1,1,1] -> 1 (sign)
            @test IntU.character_at_id([1, 1, 1]) == 1
            # [2,1] -> 2 (standard)
            @test IntU.character_at_id([2, 1]) == 2

            # S4
            # [4] -> 1
            @test IntU.character_at_id([4]) == 1
            # [3,1] -> 3
            @test IntU.character_at_id([3, 1]) == 3
            # [2,2] -> 2
            @test IntU.character_at_id([2, 2]) == 2
            # [2,1,1] -> 3
            @test IntU.character_at_id([2, 1, 1]) == 3
        end

        @testset verbose=true "murnaghan_nakayama (Character table values)" begin
            # S3 Character Table
            # Lambda [3] (Trivial): 1, 1, 1 everywhere.
            @test IntU.murnaghan_nakayama([3], [1, 1, 1]) == 1
            @test IntU.murnaghan_nakayama([3], [2, 1]) == 1
            @test IntU.murnaghan_nakayama([3], [3]) == 1

            # Lambda [1,1,1] (Sign): 1, -1, 1
            @test IntU.murnaghan_nakayama([1, 1, 1], [1, 1, 1]) == 1
            @test IntU.murnaghan_nakayama([1, 1, 1], [2, 1]) == -1
            @test IntU.murnaghan_nakayama([1, 1, 1], [3]) == 1

            # Lambda [2,1] (Standard): 2, 0, -1
            @test IntU.murnaghan_nakayama([2, 1], [1, 1, 1]) == 2
            @test IntU.murnaghan_nakayama([2, 1], [2, 1]) == 0
            @test IntU.murnaghan_nakayama([2, 1], [3]) == -1
        end

        @testset verbose=true "Weingarten Function consistency" begin
            # Wg([1,1], d) should be 1/(d^2-1)
            d = 3
            @test IntU.weingarten([1, 1], d) == 1//(d^2-1)
        end
    end
end
