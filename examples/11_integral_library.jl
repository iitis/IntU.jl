using IntU
using Symbolics

println("=== Pre-computed Integral Library Examples ===\n")

@variables d
# 1. Haar Unitary Trace (tr(U A U' B))
# This is a very common integral in QI and RMT.
println("1. Haar Unitary Trace")
U = SymbolicMatrix(:U, :U)
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)

println("Integrating tr(U * A * U' * B) over dU(d)...")
# This returns instantly because it's in the library
res = @integrate tr(U * A * U' * B) dU(d)
println("Result: ", res)
println(" (Expected: (tr(A) * tr(B)) / d)\n")


# 2. Gaussian Unitary Ensemble (GUE) Moments
println("2. GUE Moments")
H = SymbolicMatrix(:H, :GUE)

m2 = @integrate tr(H^2) dGUE(H, d)
m4 = @integrate tr(H^4) dGUE(H, d)
m6 = @integrate tr(H^6) dGUE(H, d)

println("<Tr(H^2)>_GUE = ", simplify(m2))
println("<Tr(H^4)>_GUE = ", simplify(m4))
println("<Tr(H^6)>_GUE = ", simplify(m6), "\n")


# 3. Gaussian Orthogonal (GOE) and Symplectic (GSE) Moments
println("3. GOE and GSE Moments")

goe2 = @integrate tr(H^2) dGOE(H, d)
gse2 = @integrate tr(H^2) dGSE(H, d)

println("<Tr(H^2)>_GOE = ", simplify(goe2))
println("<Tr(H^2)>_GSE = ", simplify(gse2))

println("\nDone.")
