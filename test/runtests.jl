using QUBOConstraints
using Test

@testset "QUBOConstraints.jl" begin
    include("conversion.jl")
    include("learn.jl")
end
