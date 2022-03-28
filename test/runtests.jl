using QUBOConstraints
using Test

@testset "QUBOConstraints.jl" begin
    include("encoding.jl")
    include("learn.jl")
end
