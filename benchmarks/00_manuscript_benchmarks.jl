using IntU
using Symbolics
using LinearAlgebra
using BenchmarkTools
using Memoization
using Printf
import ITensors

# Helper to run benchmark and return (median time in ms, memory in MiB)
function measure_median_func(f)
    # Benchmark the function call. We interpolate f to avoid overhead of finding it.
    b =
        @benchmark $f() evals=1 samples=30 seconds=120 setup=(Memoization.empty_all_caches!())
    m = median(b)
    return m.time / 1e6, m.memory / (1024 * 1024)
end

println("========================================================================")
println("Reproducing Manuscript Benchmarks (IntU.jl)")
println("========================================================================")
println("")

# --- TABLE 1: Symbolic vs Numeric ---

println("Table 1: Benchmark results for symbolic vs. numeric integration")
println("------------------------------------------------------------------------")
@printf("%-18s %-25s %-15s %-10s\n", "Group", "Integrand", "Dimension", "Time (ms)")
println("------------------------------------------------------------------------")

# Unitary
@variables d
U = SymbolicMatrix(:U, :U)
measure_sym = dU(d)

# |U_11|^6, Symbolic
t, _ = measure_median_func(() -> integrate(abs(U[1, 1])^6, measure_sym))
@printf("%-18s %-25s %-15s %10.2f\n", "Unitary", "|U_11|^6", "Symbolic", t)

# |U_11|^8, Symbolic
t, _ = measure_median_func(() -> integrate(abs(U[1, 1])^8, measure_sym))
@printf("%-18s %-25s %-15s %10.2f\n", "Unitary", "|U_11|^8", "Symbolic", t)

# |U_11|^10, Symbolic
t, _ = measure_median_func(() -> integrate(abs(U[1, 1])^10, measure_sym))
@printf("%-18s %-25s %-15s %10.2f\n", "Unitary", "|U_11|^10", "Symbolic", t)

# |U_11|^10, d=10
U10 = SymbolicMatrix(:U, :U, 10)
m10 = dU(10)
t, _ = measure_median_func(() -> integrate(abs(U10[1, 1])^10, m10))
@printf("%-18s %-25s %-15s %10.2f\n", "Unitary", "|U_11|^10", "d=10", t)

# |U_11|^10, d=50
U50 = SymbolicMatrix(:U, :U, 50)
m50 = dU(50)
t, _ = measure_median_func(() -> integrate(abs(U50[1, 1])^10, m50))
@printf("%-18s %-25s %-15s %10.2f\n", "Unitary", "|U_11|^10", "d=50", t)

# |tr(U)|^4, Symbolic
t, _ = measure_median_func(() -> integrate(abs(tr(U))^4, measure_sym))
@printf("%-18s %-25s %-15s %10.2f\n", "Unitary", "|tr(U)|^4", "Symbolic", t)

# |tr(U)|^10, d=10
t, _ = measure_median_func(() -> integrate(abs(tr(U10))^10, m10))
@printf("%-18s %-25s %-15s %10.2f\n", "Unitary", "|tr(U)|^10", "d=10", t)

# |tr(U)|^14, d=10
t, _ = measure_median_func(() -> integrate(abs(tr(U10))^14, m10))
@printf("%-18s %-25s %-15s %10.2f\n", "Unitary", "|tr(U)|^14", "d=10", t)

println("------------------------------------------------------------------------")

# Orthogonal
O = SymbolicMatrix(:O, :O)
mO_sym = dO(d)

# O_11^2, Symbolic
t, _ = measure_median_func(() -> integrate(O[1, 1]^2, mO_sym))
@printf("%-18s %-25s %-15s %10.2f\n", "Orthogonal", "O_11^2", "Symbolic", t)

# O_11^4, Symbolic
t, _ = measure_median_func(() -> integrate(O[1, 1]^4, mO_sym))
@printf("%-18s %-25s %-15s %10.2f\n", "Orthogonal", "O_11^4", "Symbolic", t)


# O_11^6, d=10
O10 = SymbolicMatrix(:O, :O, BigInt(10))
mO10 = dO(BigInt(10))
t, _ = measure_median_func(() -> integrate(O10[1, 1]^6, mO10))
@printf("%-18s %-25s %-15s %10.2f\n", "Orthogonal", "O_11^6", "d=10", t)

# O_11^8, d=20
O20 = SymbolicMatrix(:O, :O, BigInt(20))
mO20 = dO(BigInt(20))
t, _ = measure_median_func(() -> integrate(O20[1, 1]^8, mO20))
@printf("%-18s %-25s %-15s %10.2f\n", "Orthogonal", "O_11^8", "d=20", t)

# O_11^10, d=20
t, _ = measure_median_func(() -> integrate(O20[1, 1]^10, mO20))
@printf("%-18s %-25s %-15s %10.2f\n", "Orthogonal", "O_11^10", "d=20", t)


# O_11^10, d=50
O50 = SymbolicMatrix(:O, :O, BigInt(50))
mO50 = dO(BigInt(50))
t, _ = measure_median_func(() -> integrate(O50[1, 1]^10, mO50))
@printf("%-18s %-25s %-15s %10.2f\n", "Orthogonal", "O_11^10", "d=50", t)

println("------------------------------------------------------------------------")

# Symplectic
mSp10 = dSp(BigInt(10))
S10 = SymbolicMatrix(:S, :Sp, BigInt(10)) # Even dimension

# |S_11|^8, d=10. Integration of unitary symplectic matrix elements.
t, _ = measure_median_func(() -> integrate(abs(S10[1, 1])^8, mSp10))
@printf("%-18s %-25s %-15s %10.2f\n", "Symplectic", "|S_11|^8", "d=10", t)

# |S_11|^10, d=10
t, _ = measure_median_func(() -> integrate(abs(S10[1, 1])^10, mSp10))
@printf("%-18s %-25s %-15s %10.2f\n", "Symplectic", "|S_11|^10", "d=10", t)

# |S_11|^10, d=20
S20 = SymbolicMatrix(:S, :Sp, BigInt(20))
mSp20 = dSp(BigInt(20))
t, _ = measure_median_func(() -> integrate(abs(S20[1, 1])^10, mSp20))
@printf("%-18s %-25s %-15s %10.2f\n", "Symplectic", "|S_11|^10", "d=20", t)

println("------------------------------------------------------------------------")

# GinUE
G = SymbolicMatrix(:G, :GinUE)
mG_sym = dGinUE(d)
mG4 = dGinUE(4)

# tr(G G') Symbolic
# tr_lazy?
t, _ = measure_median_func(() -> integrate(tr_lazy(G * G'), mG_sym))
@printf("%-18s %-25s %-15s %10.2f\n", "GinUE", "tr(GG')", "Symbolic", t)

# tr(G G') d=4
t, _ = measure_median_func(() -> integrate(tr_lazy(G * G'), mG4))
@printf("%-18s %-25s %-15s %10.2f\n", "GinUE", "tr(GG')", "d=4", t)

println("------------------------------------------------------------------------")

# Circular Ensembles
# Symmetric S (COE)
S_coe = SymbolicMatrix(:S, :COE)
mCOE = dCOE(d)

t, _ = measure_median_func(() -> integrate(abs(S_coe[1, 1])^2, mCOE))
@printf("%-18s %-25s %-15s %10.2f\n", "Circ. Orthogonal", "|S_11|^2", "Symbolic", t)
t, _ = measure_median_func(() -> integrate(abs(S_coe[1, 1])^4, mCOE))
@printf("%-18s %-25s %-15s %10.2f\n", "Circ. Orthogonal", "|S_11|^4", "Symbolic", t)
# |S_11|^6, Symbolic
t, _ = measure_median_func(() -> integrate(abs(S_coe[1, 1])^6, mCOE))
@printf("%-18s %-25s %-15s %10.2f\n", "Circ. Orthogonal", "|S_11|^6", "Symbolic", t)
# @printf("%-18s %-25s %-15s %10s\n", "Circ. Orthogonal", "|S_11|^6", "Symbolic", "Skipped")


# Self-dual S (CSE)
S_cse = SymbolicMatrix(:S, :CSE)
mCSE = dCSE(d)
t, _ = measure_median_func(() -> integrate(abs(S_cse[1, 1])^2, mCSE))
@printf("%-18s %-25s %-15s %10.2f\n", "Circ. Symplectic", "|S_11|^2", "Symbolic", t)
t, _ = measure_median_func(() -> integrate(abs(S_cse[1, 1])^4, mCSE))
@printf("%-18s %-25s %-15s %10.2f\n", "Circ. Symplectic", "|S_11|^4", "Symbolic", t)
# |S_11|^6, Symbolic
t, _ = measure_median_func(() -> integrate(abs(S_cse[1, 1])^6, mCSE))
@printf("%-18s %-25s %-15s %10.2f\n", "Circ. Symplectic", "|S_11|^6", "Symbolic", t)
# @printf("%-18s %-25s %-15s %10s\n", "Circ. Symplectic", "|S_11|^6", "Symbolic", "Skipped")


println("------------------------------------------------------------------------")

# Permutation
mP100 = dPerm(100)
P100 = SymbolicMatrix(:P, :Perm, 100) # Permutation matrices are real
t, _ = measure_median_func(() -> integrate(P100[1, 1]^10, mP100))
@printf("%-18s %-25s %-15s %10.2f\n", "Permutation", "P_11^10", "d=100", t)

# tr(PA)^2, d=4
mP4 = dPerm(4)
P4 = SymbolicMatrix(:P, :Perm, 4)
A = SymbolicMatrix(:A, :Constant, 4)
t_trpa, _ = measure_median_func(() -> integrate(tr_lazy(P4*A)^2, mP4))
@printf("%-18s %-25s %-15s %10.2f\n", "Permutation", "tr(PA)^2", "d=4", t_trpa)

# CP (Centered)
Y10 = SymbolicMatrix(:Y, :CPerm, 10)
mCP10 = dCPerm(10)
t, _ = measure_median_func(() -> integrate(Y10[1, 1]^4, mCP10))
@printf("%-18s %-25s %-15s %10.2f\n", "Centered Perm.", "Y_11^4", "d=10", t)

println("------------------------------------------------------------------------")

function benchmark_bipartite()
    d_total = 6
    dims = (3, 2)
    U = SymbolicMatrix(:U, :U, d_total)
    rho_fixed = zeros(Num, d_total, d_total)
    rho_fixed[1, 1] = 1
    rho_random = U * rho_fixed * U'
    measure = dU(d_total)

    t, _ = measure_median_func(() -> begin
        rho_A = partial_trace(rho_random, dims, 2)
        integrate(tr(rho_A * rho_A), measure)
    end)
    return t
end
t_pure = benchmark_bipartite()
@printf("%-18s %-25s %-15s %10.2f\n", "Application", "Bipartite", "d=6", t_pure)


# --- TABLE 2: ITensor Scaling ---
println("\n\nTable 2: ITensor network integration (Haar measure)")
println("------------------------------------------------------------------------")
@printf("%-20s %-10s %-10s %-10s\n", "Scaling Type", "Degree k", "Dim d", "Time (ms)")
println("------------------------------------------------------------------------")

function run_itensor_bench(k, d_val)
    # Placeholder for consistency, actual calls moved to loop
end

# Copy create_trace_network from 10_itensor_integration.jl
function create_trace_network(dim, k, measure_type = :U)
    out_indices = [ITensors.Index(dim, "Out,$i") for i = 1:k]
    in_indices = [ITensors.Index(dim, "In,$i") for i = 1:k]

    tensors = Any[]

    if measure_type == :U
        for i = 1:k
            U = IntU.ITensorUnitary(
                out_indices = [out_indices[i]],
                in_indices = [in_indices[i]],
                is_adj = false,
            )

            U_dag = IntU.ITensorUnitary(
                out_indices = [in_indices[i]],
                in_indices = [out_indices[i]],
                is_adj = true,
            )

            push!(tensors, U)
            push!(tensors, U_dag)
        end
        # Create a trace by connecting Out[i] to In[i+1] and In[i] to Out[i+1] etc.
        for i = 1:k
            next_i = (i % k) + 1
            A = ITensors.randomITensor(in_indices[i], out_indices[next_i])
            B = ITensors.randomITensor(in_indices[next_i], out_indices[i])
            push!(tensors, A)
            push!(tensors, B)
        end
        return tensors, dU(dim)
    elseif measure_type == :O
        for i = 1:k
            O = IntU.ITensorUnitary(
                out_indices = [out_indices[i]],
                in_indices = [in_indices[i]],
                is_adj = false,
            )
            push!(tensors, O)
        end
        # Create a trace by connecting Out[i] to In[i+1]
        for i = 1:k
            next_i = (i % k) + 1
            A = ITensors.randomITensor(in_indices[i], out_indices[next_i])
            push!(tensors, A)
        end
        return tensors, dO(dim)
    elseif measure_type == :Sp
        for i = 1:k
            S = IntU.ITensorUnitary(
                out_indices = [out_indices[i]],
                in_indices = [in_indices[i]],
                is_adj = false,
            )
            push!(tensors, S)
        end
        # Create a trace for Symplectic
        for i = 1:k
            next_i = (i % k) + 1
            A = ITensors.randomITensor(in_indices[i], out_indices[next_i])
            push!(tensors, A)
        end
        return tensors, dSp(dim)
    end
end

# Restore ITensor scaling benchmarks
for k_val in [1, 2, 3, 4]
    local t_it, _ = measure_median_func(() -> begin
        tensors, measure = create_trace_network(2, k_val, :U)
        integrate(tensors, measure)
    end)
    @printf("%-20s %-10d %-10d %10.2f\n", "Degree k (U)", k_val, 2, t_it)
end

println("------------------------------------------------------------------------")
for k_val in [2, 4, 6]
    local t_it_o, _ = measure_median_func(() -> begin
        tensors, measure = create_trace_network(3, k_val, :O)
        integrate(tensors, measure)
    end)
    @printf("%-20s %-10d %-10d %10.2f\n", "Degree k (O)", k_val, 3, t_it_o)
end

println("------------------------------------------------------------------------")
for d_val in [2, 10, 50, 100]
    local t_it, _ = measure_median_func(() -> begin
        tensors, measure = create_trace_network(d_val, 2, :U)
        integrate(tensors, measure)
    end)
    @printf("%-20s %-10d %-10d %10.2f\n", "Dimension d (k=2)", 2, d_val, t_it)
end

println("------------------------------------------------------------------------")
for d_val in [2, 10, 20, 30]
    local t_it, _ = measure_median_func(() -> begin
        tensors, measure = create_trace_network(d_val, 3, :U)
        integrate(tensors, measure)
    end)
    @printf("%-20s %-10d %-10d %10.2f\n", "Dimension d (k=3)", 3, d_val, t_it)
end

println("------------------------------------------------------------------------")
# Orthogonal k=6, d=3
t_it_o6, _ = measure_median_func(() -> begin
    tensors, measure = create_trace_network(3, 6, :O)
    integrate(tensors, measure)
end)
@printf("%-20s %-10d %-10d %10.2f\n", "Orthogonal", 6, 3, t_it_o6)


# --- TABLE 3: Matrix Integration ---
println("\n\nTable 3: Matrix integration E[U U'] over U(d)")
println("------------------------------------------------------------------------")
@printf("%-15s %-15s %-20s\n", "Matrix Size", "Time (ms)", "Allocations (MiB)")
println("------------------------------------------------------------------------")

function bench_matrix(N)
    d_val = N
    U_sym = SymbolicMatrix(:U, :U, d_val)
    m = dU(d_val)
    expr = U_sym * U_sym' # This creates an N x N matrix of expressions if d_val is concrete

    t, mem = measure_median_func(() -> integrate(expr, m))
    @printf("%-15s %10.2f %20.2f\n", "$N x $N", t, mem)
end

bench_matrix(2)
bench_matrix(3)
bench_matrix(4)

println("------------------------------------------------------------------------")
