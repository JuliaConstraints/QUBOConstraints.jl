module QUBOConstraints

# SECTION - Usings and imports
using ConstraintDomains
using LinearAlgebra
using TestItemRunner
using TestItems

# SECTION - exports
export binarize
export debinarize
export is_valid
export train

export QUBO_linear_sum

# SECTION - includes
include("base.jl")

include("handmade/linear_sum.jl")

include("encoding/domain_wall.jl")
include("encoding/one_hot.jl")
include("encoding/conversion.jl")

include("learn.jl")

end
