module ConstraintsExt
using QUBOConstraints
import Constraints

function QUBOConstraints.constraint_oracle(name::Symbol; kwargs...)
    predicate = Constraints.concept(name)
    parameters = deepcopy((; kwargs...))
    return x -> predicate(x; parameters...)
end
end
