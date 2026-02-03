# Example 13: ITensors.jl Integration Showcase
# ============================================
# This example demonstrates all major features of the IntU.jl + ITensors.jl bridge.
# It covers Haar, Orthogonal, Symplectic, t-designs, and symbolic dimensions.

using IntU
using ITensors
using Symbolics

# Helper to print ITensors nicely
function print_itensor(label, T)
    println("\n--- $label ---")
    println(T)
    if order(T) == 0
        println("Value: ", scalar(T))
    end
end

println("Starting IntU + ITensors Showcase...\n")

# ==========================================================
# 1. Haar Unitary Integration: E[Tr(U A U' B)]
# ==========================================================
println("### 1. Haar Unitary Integration ###")
i = Index(2, "Out,Index1")
j = Index(2, "In,Index1")
i2 = Index(2, "Out,Index2")
j2 = Index(2, "In,Index2")

U_it = randomITensor(i, j)
U_dag_it = randomITensor(j2, i2)

# Constant tensors A and B to form a trace: Tr(U A U' B)
A = randomITensor(j, j2)
B = randomITensor(i2, i)

# Wrap unitaries
U = ITensorUnitary(U_it; out_indices=[i], in_indices=[j])
U_dag = ITensorUnitary(U_dag_it; out_indices=[j2], in_indices=[i2], is_adj=true)

res1 = integrate([U, A, U_dag, B], dU(2))
print_itensor("Haar Result (scalar)", res1)


# ==========================================================
# 2. Orthogonal Group Integration: E[O_11 O_22]
# ==========================================================
println("\n### 2. Orthogonal Group Integration ###")
o1 = Index(3, "Out1")
i1 = Index(3, "In1")
o2 = Index(3, "Out2")
i2_in = Index(3, "In2")

O1_it = randomITensor(o1, i1)
O2_it = randomITensor(o2, i2_in)

O1 = ITensorUnitary(O1_it; out_indices=[o1], in_indices=[i1])
O2 = ITensorUnitary(O2_it; out_indices=[o2], in_indices=[i2_in])

# Integrate over O(3)
res2 = integrate([O1, O2], dO(3))
print_itensor("Orthogonal Result", res2)


# ==========================================================
# 3. Symplectic Group Integration
# ==========================================================
println("\n### 3. Symplectic Group Integration ###")
# Sp(d) requires d to be even
d = 2
s1 = Index(d, "S1_Out")
e1 = Index(d, "S1_In")
s2 = Index(d, "S2_Out")
e2 = Index(d, "S2_In")

S1_it = randomITensor(s1, e1)
S2_it = randomITensor(s2, e2)

S1 = ITensorUnitary(S1_it; out_indices=[s1], in_indices=[e1])
S2 = ITensorUnitary(S2_it; out_indices=[s2], in_indices=[e2])

res3 = integrate([S1, S2], dSp(d))
print_itensor("Symplectic Result", res3)


# ==========================================================
# 4. Unitary t-designs: E[|U_11|^4] with Design Case
# ==========================================================
println("\n### 4. Unitary t-designs ###")
u_out = Index(2, "Out")
u_in = Index(2, "In")
U_it = randomITensor(u_out, u_in)

U_wrap = ITensorUnitary(U_it; out_indices=[u_out], in_indices=[u_in])
U_dag_wrap = ITensorUnitary(U_it; out_indices=[u_in], in_indices=[u_out], is_adj=true) # Simplified dag for example

# 2-design should match Haar for degree 2 polynomials
measure_2 = dDesign(nothing, 2, 2)
res_haar = integrate([U_wrap, U_dag_wrap], dU(2))
res_design = integrate([U_wrap, U_dag_wrap], measure_2)

println("Haar Integral E[|U_11|^2]: ", scalar(res_haar))
println("2-Design Integral: ", scalar(res_design))


# ==========================================================
# 5. Symbolic Dimensions
# ==========================================================
println("\n### 5. Symbolic Dimensions ###")
@variables d_sym
i_sym = Index(2, "Out") # ITensors needs a numeric size for the object
j_sym = Index(2, "In")
U_it_sym = randomITensor(i_sym, j_sym)

U_sym_wrap = ITensorUnitary(U_it_sym; out_indices=[i_sym], in_indices=[j_sym])
U_dag_sym_wrap = ITensorUnitary(U_it_sym; out_indices=[j_sym], in_indices=[i_sym], is_adj=true)

# Integrate over U(d_sym) - just 1 moment for simplicity
res_sym = integrate([U_sym_wrap, U_dag_sym_wrap], dDesign(nothing, d_sym, 1))
print_itensor("Symbolic Dimension Result", res_sym)


# ==========================================================
# 6. Nested Integration: E_V [ E_U [ U U' V V' ] ]
# ==========================================================
println("\n### 6. Nested Integration ###")
# Setup indices for two independent unitaries U and V
idx_u_out = Index(2, "U,Out")
idx_u_in = Index(2, "U,In")
idx_v_out = Index(2, "V,Out")
idx_v_in = Index(2, "V,In")

U_it = randomITensor(idx_u_out, idx_u_in)
V_it = randomITensor(idx_v_out, idx_v_in)

U_w = ITensorUnitary(U_it; out_indices=[idx_u_out], in_indices=[idx_u_in])
U_dag_w = ITensorUnitary(U_it; out_indices=[idx_u_in], in_indices=[idx_u_out], is_adj=true)

V_w = ITensorUnitary(V_it; out_indices=[idx_v_out], in_indices=[idx_v_in])
V_dag_w = ITensorUnitary(V_it; out_indices=[idx_v_in], in_indices=[idx_v_out], is_adj=true)

println("Integrating E_V [ E_U [ (U ⊗ U') ⊗ (V ⊗ V') ] ]")

# Step 1: Integrate over U first. V parts are treated as constants.
res_inner = integrate([U_w, U_dag_w, V_it, V_dag_w.tensor], dU(2))
# res_inner is an ITensor since it contains the unintegrated V parts.

# Step 2: Integrate the result over V
# Since res_inner grew into a larger tensor after Step 1, we pull out the V unitaries.
res_final = integrate([V_w, V_dag_w, res_inner], dU(2))
print_itensor("Nested Integration Result", res_final)

println("\nShowcase completed.")
