function default_sat(concept; encoding = :domain_wall)
    # Invalid binary representation
    check_validity(x, Val(encoding)) && return sqrt(length(X))
    # Return -1 if the configuration is a solution, 1 otherwise
    return concept(integerize(X)) ? -1. : 1.
end

default_opt(x...) = 0.0
