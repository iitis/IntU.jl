using IntegrateUnitary
using Symbolics
using LinearAlgebra

println("Testing symbolic dimension d with unified clean interface...")

# 1. Integrate |U_11|^2
println("\n1. Integrating |U_11|^2 ...")
res1 = @integrate abs(U[1, 1])^2 dU(d)
println("Result: $res1 (Expected: 1/d)")

# 2. Integrate |U_11|^4
println("\n2. Integrating |U_11|^4 ...")
res2 = @integrate abs(U[1, 1])^4 dU(d)
# Expected: 2 / (d*(d+1))
println("Result: $res2")

# 3. Integrate a 2x2 minor
println("\n3. Integrating |U_11*U_22 - U_12*U_21|^2 ...")
res3 = @integrate abs(U[1, 1]*U[2, 2] - U[1, 2]*U[2, 1])^2 dU(d)
# Theoretical expectation: 2 / (d*(d-1))
println("Result: $res3")

# 4. High-degree moment (Example A: different rows/cols)
println("\n4. Integrating |U_11|^2 * |U_22|^4 * |U_13|^6 (Example A) ...")
res4 = @integrate abs(U[1, 1])^2 * abs(U[2, 2])^4 * abs(U[1, 3])^6 dU(d)
println("Result: $res4")

# 5. High-degree moment (Example B: same row, Dirichlet behavior)
println("\n5. Integrating |U_11|^2 * |U_12|^4 * |U_13|^6 (Example B) ...")
res5 = @integrate abs(U[1, 1])^2 * abs(U[1, 2])^4 * abs(U[1, 3])^6 dU(d)
println("Result: $res5")

# 6. "2-design style" correlator (Example C)
println("\n6. Integrating U_11 * conj(U_12) * U_22 * conj(U_21) (Example C) ...")
res6 = @integrate U[1, 1] * conj(U[1, 2]) * U[2, 2] * conj(U[2, 1]) dU(d)
println("Result: $res6")

# 7. Orthogonal group (Example O1: high powers, single row)
println("\n7. Integrating O_11^2 * O_12^4 * O_13^6 (Example O1) ...")
res7 = @integrate O[1, 1]^2 * O[1, 2]^4 * O[1, 3]^6 dO(d)
println("Result: $res7")

# 8. Orthogonal group (Example O2: mixed rows/cols)
println("\n8. Integrating O_11^2 * O_22^4 * O_13^6 (Example O2) ...")
res8 = @integrate O[1, 1]^2 * O[2, 2]^4 * O[1, 3]^6 dO(d)
println("Result: $res8")

# 9. Symplectic group (Example S1: Symplectic mixed moments)
println("\n9. Integrating |Sp_11|^2 * |Sp_12|^2 (Example S1) ...")
res9 = @integrate abs(Sp[1, 1])^2 * abs(Sp[1, 2])^2 dSp(d)
println("Result: $res9")

# 10. Symplectic group (Example Sp1: high powers, single row)
println("\n10. Integrating |Sp_11|^2 * |Sp_12|^4 * |Sp_13|^6 (Example Sp1) ...")
res10 = @integrate abs(Sp[1, 1])^2 * abs(Sp[1, 2])^4 * abs(Sp[1, 3])^6 dSp(d)
println("Result: $res10")

# 11. GUE Trace Moments
println("\n11. Integrating tr(H^4) and tr(H^6) (GUE) ...")
res_gue4 = @integrate tr(H^4) dGUE(d)
res_gue6 = @integrate tr(H^6) dGUE(d)
println("Tr(H^4): $res_gue4")
println("Tr(H^6): $res_gue6")

# 12. GinUE Identiy
println("\n12. Integrating tr(G * G' * G * G') (GinUE) ...")
res_ginue_sq = @integrate tr(G * G' * G * G') dGinUE(d)
println("Tr(G G' G G'): $res_ginue_sq")

# 13. Wishart-style GinUE Moments
println("\n13. Integrating Wishart-style GinUE moments ...")
x2 = @integrate tr(G * G')^2 dGinUE(d)
y2 = @integrate tr((G * G')^2) dGinUE(d)
println("Tr(G G')^2: $x2")
println("Tr((G G')^2): $y2")

# 14. Circular ensembles (COE)
println("\n14. Integrating Circular Orthogonal Ensemble (COE) moments ...")
diag2 = @integrate abs(S[1, 1])^2 dCOE(d)
off2 = @integrate abs(S[1, 2])^2 dCOE(d)
diag4 = @integrate abs(S[1, 1])^4 dCOE(d)
off4 = @integrate abs(S[1, 2])^4 dCOE(d)
println("COE Diag 2nd: $diag2")
println("COE Off-diag 2nd: $off2")
println("COE Diag 4th: $diag4")
println("COE Off-diag 4th: $off4")

# 15. COE Correlation moment
println("\n15. Integrating COE correlation moment ...")
mix = @integrate abs(S[1, 1])^2 * abs(S[1, 2])^2 dCOE(d)
println("COE Correlation: $mix")

# 16. Circular ensembles (CSE)
println("\n16. Integrating Circular Symplectic Ensemble (CSE) moments ...")
cse_diag2 = @integrate abs(S[1, 1])^2 dCSE(d)
cse_diag4 = @integrate abs(S[1, 1])^4 dCSE(d)
println("CSE Diag 2nd: $cse_diag2")
println("CSE Diag 4th: $cse_diag4")

println("\nSimplified results:")
println("1. $(Symbolics.simplify(res1))")
println("2. $(Symbolics.simplify(res2))")
println("3. $(Symbolics.simplify(res3))")

println("\nDefining another symbolic matrix V using factory function:")
V = symbolic_unitary(:V, d)
res_v = integrate(abs(V[1, 1])^2, dU(d))
println("Result for V: $res_v")
