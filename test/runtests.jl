using TestItemRunner
using TestItems

@run_package_tests

@testitem "Learn" tags = [:qubo, :learn, :regression] default_imports=false begin
    using ConstraintDomains
    using DataFrames
    using Flux
    using PrettyTables
    using QUBOConstraints
    using Test

    import Flux.Optimise: update!
    import Flux: params

    include("gradient_descent.jl")
    include("learn.jl")
end
