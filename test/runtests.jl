using ConstraintDomains
using DataFrames
using Flux
using PrettyTables
using QUBOConstraints
using Test

import Flux.Optimise: update!
import Flux: params

include("gradient_descent.jl")

@testset "QUBOConstraints.jl" begin
    include("encoding.jl")
    include("learn.jl")
end
