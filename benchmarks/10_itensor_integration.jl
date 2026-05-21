using IntU
using ITensors
using BenchmarkTools
using LinearAlgebra

function create_trace_network(dim, k, measure_type = :U)
    out_indices = [Index(dim, "Out,$i") for i = 1:k]
    in_indices = [Index(dim, "In,$i") for i = 1:k]

    tensors = Any[]

    if measure_type == :U
        for i = 1:k
            U = ITensorUnitary(
                out_indices = [out_indices[i]],
                in_indices = [in_indices[i]],
                is_adj = false,
            )

            U_dag = ITensorUnitary(
                out_indices = [in_indices[i]],
                in_indices = [out_indices[i]],
                is_adj = true,
            )

            push!(tensors, U)
            push!(tensors, U_dag)
        end
        for i = 1:k
            next_i = (i % k) + 1
            A = randomITensor(in_indices[i], out_indices[next_i])
            B = randomITensor(in_indices[next_i], out_indices[i])
            push!(tensors, A)
            push!(tensors, B)
        end
        return tensors, dU(dim)
    elseif measure_type == :O
        for i = 1:k
            O = ITensorUnitary(
                out_indices = [out_indices[i]],
                in_indices = [in_indices[i]],
                is_adj = false,
            )
            push!(tensors, O)
        end
        if isodd(k)
            push!(
                tensors,
                ITensorUnitary(
                    out_indices = [out_indices[1]],
                    in_indices = [in_indices[1]],
                    is_adj = false,
                ),
            )
        end
        return tensors, dO(dim)
    elseif measure_type == :Sp
        for i = 1:k
            S = ITensorUnitary(
                out_indices = [out_indices[i]],
                in_indices = [in_indices[i]],
                is_adj = false,
            )
            push!(tensors, S)
        end
        if isodd(k)
            push!(
                tensors,
                ITensorUnitary(
                    out_indices = [out_indices[1]],
                    in_indices = [in_indices[1]],
                    is_adj = false,
                ),
            ) # Placeholder
        end
        return tensors, dSp(dim)
    end
end

function run_benchmarks()
    println("=== ITensors Integration Scaling Benchmarks ===\n")

    # 1. Scaling with Degree k (Haar Unitary)
    println("--- Haar Unitary Scaling (Degree k) ---")
    d = 2
    for k in [1, 2, 3, 4]
        println("Benchmark: U($d), degree=$k")
        tensors, measure = create_trace_network(d, k, :U)
        integrate(tensors, measure)
        t = @benchmark integrate($tensors, $measure)
        display(t)
        println()
    end

    # 2. Scaling with Dimension d (Haar Unitary)
    println("\n--- Haar Unitary Scaling (Dimension d) ---")
    # k=2 Scaling
    println("--- k = 2 ---")
    for d in [2, 10, 50, 100]
        println("Benchmark: U($d), degree=2")
        tensors, measure = create_trace_network(d, 2, :U)
        integrate(tensors, measure)
        t = @benchmark integrate($tensors, $measure)
        display(t)
        println()
    end

    # k=3 Scaling (Lower max d due to memory)
    println("--- k = 3 ---")
    for d in [2, 10, 20, 30]
        println("Benchmark: U($d), degree=3")
        tensors, measure = create_trace_network(d, 3, :U)
        integrate(tensors, measure)
        t = @benchmark integrate($tensors, $measure)
        display(t)
        println()
    end

    # 3. Orthogonal Scaling
    println("\n--- Orthogonal Scaling (Degree k) ---")
    d = 3
    for k in [2, 4, 6]
        println("Benchmark: O($d), degree=$k")
        tensors, measure = create_trace_network(d, k, :O)
        integrate(tensors, measure)
        t = @benchmark integrate($tensors, $measure)
        display(t)
        println()
    end
end

run_benchmarks()
