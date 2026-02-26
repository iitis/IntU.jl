using Test
using Aqua

@testset "Aqua.jl Quality Assurance" begin
    Aqua.test_all(IntU)
end
