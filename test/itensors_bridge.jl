
using Test
using IntU
using Symbolics

struct MockIndex
    id::Int
    dim::Int
end

struct MockTensor
    indices::Vector{MockIndex}
    name::Symbol
end

function IntU._contract_all(cs::Vector)
    if isempty(cs)
        return 1
    end
    names = [c.name for c in cs]
    return Symbol("Contracted(", join(names, ","), ")")
end

function IntU._create_deltas(idxs1::Vector, idxs2::Vector)
    return [(idxs1[i], idxs2[i]) for i = 1:length(idxs1)]
end

const mock_mode = Ref(:Unitary)

function IntU._contract_with_deltas(cs, ds, wg)
    if mock_mode[] == :Unitary
        return string(wg) * " * Tr(A)Tr(B)"
    elseif mock_mode[] == :Orthogonal
        return string(wg) * " * Tr(O1 O2)"
    elseif mock_mode[] == :Symplectic
        return string(wg) * " * Tr(O1 O2)"
    else
        return string(wg) * " * ResultWithDeltas"
    end
end

@testset "ITensors Bridge Mock" begin
    # Test Tr(U A U' B) integration
    # Legs:
    # U: out=1, in=2
    # A: out=2, in=3
    # U_dag: out=3, in=4
    # B: out=4, in=1

    idx1 = MockIndex(1, 2)
    idx2 = MockIndex(2, 2)
    idx3 = MockIndex(3, 2)
    idx4 = MockIndex(4, 2)

    U = GraphicalUnitary([idx1], [idx2], false)
    A = MockTensor([idx2, idx3], :A)
    U_dag = GraphicalUnitary([idx3], [idx4], true)
    B = MockTensor([idx4, idx1], :B)

    constants = [A, B]
    unitaries = [U, U_dag]

    mock_mode[] = :Unitary

    res = integrate_graphical(constants, unitaries, dU(2)) # dim=2

    @test res == "1//2 * Tr(A)Tr(B)"

    # 2. Test Orthogonal
    # Tr(O A O B)
    idx1 = MockIndex(1, 2)
    idx2 = MockIndex(2, 2)
    O1 = GraphicalUnitary([idx1], [idx2], false)
    O2 = GraphicalUnitary([idx2], [idx1], false) # Loop

    mock_mode[] = :Orthogonal

    # Orthogonal integration for n=2
    res_o = integrate_graphical([], [O1, O2], dO(2))
    @test contains(res_o, "Tr(O1 O2)")

    # 3. Test Symplectic
    function IntU._create_deltas_symplectic(idxs1, idxs2, dim)
        return ["J($dim)"]
    end

    mock_mode[] = :Symplectic

    res_sp = integrate_graphical([], [O1, O2], dSp(2))
    @test contains(res_sp, "Tr(O1 O2)")
end
