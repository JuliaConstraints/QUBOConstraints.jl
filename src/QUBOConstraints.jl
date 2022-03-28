module QUBOConstraints

# SECTION - Usings and imports
using ConstraintDomains
using DataFrames
using Flux
using LinearAlgebra
using PrettyTables

import Flux.Optimise: update!

# SECTION - exports
export binarize
export integerize
export train!

# SECTION - includes
include("base.jl")

include("encoding/domain_wall.jl")
include("encoding/one_hot.jl")
include("encoding/conversion.jl")

include("learn.jl")

end
