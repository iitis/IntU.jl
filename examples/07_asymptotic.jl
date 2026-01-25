# examples/07_asymptotic.jl
using IntU
using Symbolics

# 1. Define symbolic dimension and variables
@variables d
@variables U[1:2, 1:2]::Complex
measure = dU(U, d)

println("=== Asymptotic Expansions in 1/d ===")

# 2. Simple expansion: |U_{11}|^2
# Exact value: 1/d
# Expansion: 1/d
expr1 = abs(U[1,1])^2
println("\n1. Expanding |U_{11}|^2")
res1 = asymptotic(expr1, measure, 2)
println("Expansion (order 2): ", res1)

# 3. Higher power: |U_{11}|^4
# Exact value: 2 / (d^2 - 1)
# Expansion should be: 2/d^2 + 2/d^4 + ...
expr2 = abs(U[1,1])^4
println("\n2. Expanding |U_{11}|^4")
res2 = asymptotic(expr2, measure, 4)
println("Expansion (order 4): ", res2)

# 4. Pure state expansion
println("\n3. Pure State Expansion: |psi_1|^2")
@variables psi[1:2]::Complex
measure_psi = dPsi(psi, d)
expr3 = abs(psi[1])^2
res3 = asymptotic(expr3, measure_psi, 2)
println("Expansion (order 2): ", res3)

# 5. Numeric dimension handling
println("\n4. Numeric Dimension Handling (d=3)")
# Even with numeric d, asymptotic returns a symbolic expansion in a dummy variable
measure_numeric = dU(U, 3)
res4 = asymptotic(expr2, measure_numeric, 4)
println("Expansion (order 4): ", res4)
