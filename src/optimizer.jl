abstract type AbstractOptimizer end

struct GradientDescentOptimizer <: AbstractOptimizer
    binarization::Symbol
    η::Float64
    precision::Int
    oversampling::Bool
end

function GradientDescentOptimizer(;
    binarization = :one_hot,
    η = .001,
    precision = 5,
    oversampling = false,
)
    return GradientDescentOptimizer(binarization, η, precision, oversampling)
end
