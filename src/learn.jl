"""
    AbstractOptimizer

An abstract type (interface) used to learn QUBO matrices from constraints. Only a `train` method is required.
"""
abstract type AbstractOptimizer end

"""
    train(args...)

Default `train` method for any AbstractOptimizer.
"""
train(args...) = nothing
