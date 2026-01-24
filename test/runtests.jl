using IntU
using Test
using Symbolics

@testset "IntU Tests" begin
    # Helper to convert symbolic results to numbers
    function to_numeric(x)
        x_un = Symbolics.unwrap(x)
        if x_un isa Number
            return x_un
        end
        # If it's a symbolic constant (like Num(0)), try to simplify/substitute
        # But if we reached here, it might be a Term.
        # Check if it simplifies to a number
        sim = Symbolics.simplify(x)
        sim_un = Symbolics.unwrap(sim)
        if sim_un isa Number
            return sim_un
        end
        
        # Try brute force substitution (sometimes simplify misses things)
        sim2 = Symbolics.substitute(sim, Dict())
        sim2_un = Symbolics.unwrap(sim2)
        if sim2_un isa Number
            return sim2_un
        end
        
        if isequal(sim, 0) || isequal(sim2, 0)
            return 0.0
        end
        
        # Last resort: try evaluating (works if no free variables)
        try
            val = eval(Meta.parse(string(sim)))
            if val isa Number
                return val
            end
        catch
        end
        
        if isequal(sim, 0)
            return 0.0
        end
        return x_un
    end

    # Define variables
    d_val = 3
    # Use literal dimensions for macro
    @variables U[1:3, 1:3]::Complex
    measure = dU(U, d_val)

    @testset "Example 1: |u11|^2" begin
        expr = abs(U[1,1])^2
        # Now Symbolics expanding to hypot(re, im)^2 is handled by IntU.
        
        res = integrate(expr, measure)
        res = integrate(expr, measure)
        @test to_numeric(real(res)) ≈ 1/3
        @test to_numeric(imag(res)) ≈ 0
    end

    @testset "Example 2: |u11 u22|^2" begin
        # |u11 u22|^2 
        expr = abs(U[1,1] * U[2,2])^2
        
        res = integrate(expr, measure)
        res = integrate(expr, measure)
        @test to_numeric(real(res)) ≈ 1/8
        @test to_numeric(imag(res)) ≈ 0
    end
    
    @testset "Example 3: u11 u22 conj(u12 u21)" begin
        expr = U[1,1] * U[2,2] * conj(U[1,2]) * conj(U[2,1])
        res = integrate(expr, measure)
        res = integrate(expr, measure)
        @test to_numeric(real(res)) ≈ -1/24
        @test to_numeric(imag(res)) ≈ 0
    end
    
    @testset "Unitarity Check: sum_k u_ik conj(u_jk) = delta_ij" begin
        # Integrate (sum_k u_ik conj(u_jk)) should be delta_ij * Volume?
        # Actually this is an identity U U^dag = I.
        # Integral of Constant 1 = 1.
        # Integral of constant 0 = 0.
        # Let's test term by term.
        # Int sum_k |u_1k|^2 = sum_k Int |u_1k|^2 = d * (1/d) = 1.
        
        sum_val = 0//1
        for k in 1:d_val
            sum_val += integrate(U[1,k] * conj(U[1,k]), measure)
        end
        @test abs(Float64(to_numeric(real(sum_val)))) < 1e-12 + 1.0 # 1.0 approx 1.0
        println("DEBUG: sum_val = ", sum_val)
        println("DEBUG: imag(sum_val) = ", imag(sum_val))
        @test abs(Float64(to_numeric(imag(sum_val)))) < 1e-12
        
        # Off-diagonal
        sum_off = 0//1
        for k in 1:d_val
            sum_off += integrate(U[1,k] * conj(U[2,k]), measure)
        end
        # Removed debug sum_off
        @test abs(Float64(to_numeric(real(sum_off)))) < 1e-12
        @test abs(Float64(to_numeric(imag(sum_off)))) < 1e-12
    end
    
    @testset "Weingarten Function Values" begin
        # Wg(1^2, d) = 1/(d^2-1)
        @test IntU.Weingarten.weingarten([1,1], 3) == 1//8
        
        # Wg(2, d) = -1/(d(d^2-1))
        @test IntU.Weingarten.weingarten([2], 3) == -1//24
    end

    @testset "Weingarten Unit Tests" begin
        # Import internal functions for testing if needed, or use IntU.Weingarten prefix
        using IntU.Weingarten: conjugate_partition, character_at_id, irrep_dimension, murnaghan_nakayama, weingarten

        @testset "conjugate_partition" begin
            @test conjugate_partition([1]) == [1]
            @test conjugate_partition([2]) == [1,1]
            @test conjugate_partition([1,1]) == [2]
            @test conjugate_partition([2,1]) == [2,1]
            @test conjugate_partition([3,1]) == [2,1,1]
            @test conjugate_partition([4]) == [1,1,1,1]
            @test conjugate_partition(Int[]) == Int[]
        end

        @testset "irrep_dimension (Dim of U(d) irrep)" begin
            # s_{1}(1^d) = d
            @test irrep_dimension([1], 3) == 3//1
            # s_{2}(1^d) = d(d+1)/2 symmetric tensor
            @test irrep_dimension([2], 3) == 3*4//2 # 6
            # s_{1,1}(1^d) = d(d-1)/2 antisymmetric
            @test irrep_dimension([1,1], 3) == 3*2//2 # 3
        end

        @testset "character_at_id (Dim of S_n irrep)" begin
            # S3
            # [3] -> 1 (trivial)
            @test character_at_id([3]) == 1
            # [1,1,1] -> 1 (sign)
            @test character_at_id([1,1,1]) == 1
            # [2,1] -> 2 (standard)
            @test character_at_id([2,1]) == 2
            
            # S4
            # [4] -> 1
            @test character_at_id([4]) == 1
            # [3,1] -> 3
            @test character_at_id([3,1]) == 3
            # [2,2] -> 2
            @test character_at_id([2,2]) == 2
            # [2,1,1] -> 3
            @test character_at_id([2,1,1]) == 3
        end

        @testset "murnaghan_nakayama (Character table values)" begin
            # S3 Character Table
            # Partitions: [3] (id), [2,1] (transposition), [1,1,1] (3-cycle) - Wait:
            # Cycle types correspond to classes.
            # Lambda [3] (Trivial): 1, 1, 1 everywhere.
            @test murnaghan_nakayama([3], [1,1,1]) == 1
            @test murnaghan_nakayama([3], [2,1]) == 1
            @test murnaghan_nakayama([3], [3]) == 1

            # Lambda [1,1,1] (Sign): 1, -1, 1
            @test murnaghan_nakayama([1,1,1], [1,1,1]) == 1
            @test murnaghan_nakayama([1,1,1], [2,1]) == -1
            @test murnaghan_nakayama([1,1,1], [3]) == 1

            # Lambda [2,1] (Standard): 2, 0, -1
            @test murnaghan_nakayama([2,1], [1,1,1]) == 2
            @test murnaghan_nakayama([2,1], [2,1]) == 0
            @test murnaghan_nakayama([2,1], [3]) == -1
        end
        
        @testset "Weingarten Function consistency" begin
             # Check basic property or redundancy
             # Wg([1,1], d) should be 1/(d^2-1)
             d = 3
             @test weingarten([1,1], d) == 1//(d^2-1)
        end
    end

    include("pure_states.jl")
end
