using IntU
using Symbolics

println("--- GUE Integration Tests ---")

# We use explicit matrix expansion for testing basic moments.
# This validates the combinatorial logic even if we don't do full symbolic-d trace reduction here.

N = 3
println("Using Matrix Size N = $N")

# Create explicit symbolic matrix
# Note: Symbolics.variable(:H, i, j) creates independent variables.
# We interpret them as entries of H.
H = [Symbolics.variable(:H, i, j) for i in 1:N, j in 1:N]

# Measure: GUE with dimension N
meas = dGUE(H, N)


# Helper for safe comparison of symbolic results
function safe_check(val, expected)
    # 1. Direct isequal
    if isequal(val, expected); return true; end
    # 2. Unwrap
    u = Symbolics.unwrap(val)
    if isequal(u, expected); return true; end
    
    # 3. Check if unwrap is number
    if u isa Number && u == expected; return true; end
    
    # 4. Try conversion to Number
    # Symbolics usually provides a way to extract value if it is a constant term
    try
        n = Symbolics.substitute(val, Dict())
        if n isa Number && n == expected; return true; end
        if Symbolics.unwrap(n) isa Number && Symbolics.unwrap(n) == expected; return true; end
    catch; end

    # 5. String based comparison as last resort (risky but effective for verification)
    if string(Symbolics.unwrap(val)) == string(expected) || string(val) == string(expected)
         return true
    end
    
    # 6. Rational handling
    # If val is 9//1 and expected is 9
    str_val = string(val)
    if endswith(str_val, "//1")
        try
            num_part = parse(Int, split(str_val, "//")[1])
            if num_part == expected; return true; end
        catch; end
    end
    # Handle "Term(9//1)"
    if occursin("//1", str_val)
         # Extract number
         # very hacky but works for the "(9//1)" case printed in log
         # Log output was: Result: (9//1)
         cleaned = replace(str_val, "(" => "", ")" => "", " " => "")
         if cleaned == string(expected) * "//1"
             return true
         end
    end
    
    println("DEBUG: comparison failed. val='$val' (type $(typeof(val))), expected='$expected'")
    return false
end

println("\n1. < Tr(H) >")
expr1 = tr(H)
res1 = integrate(expr1, meas)
println("Result: ", res1)
if !safe_check(res1, 0)
    error("Test 1 Failed: Expected 0, got $res1")
end

println("\n2. < Tr(H^2) >")
expr2 = tr(H^2)
res2 = integrate(expr2, meas)
println("Result: ", res2, " Type: ", typeof(res2))
expected2 = N^2
if !safe_check(res2, expected2)
    error("Test 2 Failed: Expected $expected2, got $res2")
end

println("\n3. < Tr(H^4) >")
# Expected: 2N^3 + N (Planar terms d^3, crossing term d)
# For N=3: 2(27) + 3 = 57
expr3 = tr(H^4)
# Note: matrix multiplication of symbolic variables can be slow for large N, but for N=3 it's fine.
res3 = integrate(expr3, meas)
println("Result: ", res3)
expected3 = 2*N^3 + N
if !safe_check(res3, expected3)
    error("Test 3 Failed: Expected $expected3, got $res3")
end

# 4. Mixed component test
println("\n4. Component-wise checks")
# < H_11^2 > = 1
e4a = H[1,1]^2
r4a = integrate(e4a, meas)
println("< H_{11}^2 > = ", r4a)
if !safe_check(r4a, 1)
    error("Test 4a Failed: Expected 1")
end

# < H_12 H_21 > = 1
e4b = H[1,2] * H[2,1]
r4b = integrate(e4b, meas)
println("< H_{12} H_{21} > = ", r4b)
if !safe_check(r4b, 1)
    error("Test 4b Failed: Expected 1")
end

# < H_12^2 > = 0
e4c = H[1,2]^2
r4c = integrate(e4c, meas)
println("< H_{12}^2 > = ", r4c)
if !safe_check(r4c, 0)
    error("Test 4c Failed: Expected 0")
end

println("\nAll GUE tests passed!")
