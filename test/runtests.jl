using Test
using TestItemRunner
using TestItems

@testset "Package tests: CompositionalNetworks" begin
    include("Aqua.jl")
    include("TestItemRunner.jl")
    include("qubo.jl")
end
