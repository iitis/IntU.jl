# Example 13: ITensors.jl Integration Showcase
# ============================================
# This example demonstrates the IntegrateUnitary.jl + ITensors.jl bridge.

using IntegrateUnitary
using ITensors
using Symbolics

# Helper to print ITensors nicely
function print_itensor(label, T)
    println("\n--- $label ---")
    println(T)
    if T isa ITensor
        if order(T) == 0
            println("Value: ", scalar(T))
        end
    elseif T isa Number
        println("Value: ", T)
    end
end

println("Starting IntegrateUnitary + ITensors Showcase...\n")

# ==========================================================
# 1. Haar Unitary Integration: E[Tr(U A U' B)]
# ==========================================================
println("### 1. Haar Unitary Integration ###")
i = Index(2, "Out,Index1")
j = Index(2, "In,Index1")
i2 = Index(2, "Out,Index2")
j2 = Index(2, "In,Index2")

# Constant tensors A and B to form a trace: Tr(U A U' B)
A = randomITensor(j, j2)
B = randomITensor(i2, i)

# Wrap unitaries
U = ITensorUnitary(out_indices = [i], in_indices = [j])
U_dag = ITensorUnitary(out_indices = [j2], in_indices = [i2], is_adj = true)

# The bridge supports integrate([ITensorUnitary...], measure)
# Note: we use dU(2) for a fixed dimension of 2
res1 = integrate([U, A, U_dag, B], dU(2))
print_itensor("Haar Result (scalar)", res1)


# ==========================================================
# 2. Orthogonal Group Integration
# ==========================================================
println("\n### 2. Orthogonal Group Integration ###")
o1 = Index(3, "Out1")
i1 = Index(3, "In1")
O1 = ITensorUnitary(out_indices = [o1], in_indices = [i1])

# Integrate single entry over O(3)
res2 = integrate([O1], dO(3))
print_itensor("Orthogonal Result", res2)


# ==========================================================
# 3. Symbolic Dimensions in ITensors
# ==========================================================
println("\n### 3. Symbolic Dimensions ###")
@variables d_sym
idx_out = Index(2, "Out")
idx_in = Index(2, "In")

U_sym_wrap = ITensorUnitary(out_indices = [idx_out], in_indices = [idx_in])
U_dag_sym_wrap =
    ITensorUnitary(out_indices = [idx_out], in_indices = [idx_in], is_adj = true)

# Integrate over U(d_sym)
res_sym = integrate([U_sym_wrap, U_dag_sym_wrap], dU(d_sym))
print_itensor("Symbolic Dimension Result", res_sym)

println("\nShowcase completed.")
