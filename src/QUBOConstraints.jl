module QUBOConstraints

# SECTION - Usings and imports
using ConstraintDomains
using DataFrames
using Flux
using LinearAlgebra
using PrettyTables

import Flux.Optimise: update!
import Flux: params
import DataFrames: describe

# SECTION - exports
export binarize
export debinarize
export is_valid
export train

export GradientDescentOptimizer
export QUBO_linear_sum

# SECTION - includes
include("base.jl")

include("handmade/linear_sum.jl")

include("encoding/domain_wall.jl")
include("encoding/one_hot.jl")
include("encoding/conversion.jl")

include("optimizer.jl")
include("learn.jl")

end
