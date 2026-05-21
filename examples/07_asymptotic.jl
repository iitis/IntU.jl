using IntU
using Symbolics

# 1. Define symbolic dimension and Matrix
@variables d
U = SymbolicMatrix(:U, :U)

println("=== Asymptotic Expansions in 1/d ===")

# 2. Simple expansion: |U_{11}|^2
# Exact value: 1/d
# Expansion: 1/d
println("\n1. Expanding |U_{11}|^2 over U(d)")
res1 = asymptotic(abs(U[1, 1])^2, dU(d), 2)
println("Expansion (order 2): ", res1)

# 3. Higher power: |U_{11}|^4
# Exact value: 2 / (d^2 + d)
# Expansion should be: 2/d^2 - 2/d^3 + ...
println("\n2. Expanding |U_{11}|^4 over U(d)")
res2 = asymptotic(abs(U[1, 1])^4, dU(d), 4)
println("Expansion (order 4): ", res2)

# 4. Pure state expansion
println("\n3. Pure State Expansion: |psi_1|^2")
# dPsi(d) uses the first column of a Haar unitary.
psi = symbolic_pure_state(:psi, d)
res3 = asymptotic(abs(psi[1, 1])^2, dPsi(d), 2)
println("Expansion (order 2): ", res3)

# 5. Numeric dimension handling
println("\n4. Numeric Dimension Handling (d=3)")
res4 = asymptotic(abs(U[1, 1])^4, dU(3), 4)
println("Expansion (order 4): ", res4)
