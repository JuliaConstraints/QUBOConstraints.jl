is_valid(x, encoding::Symbol = :none) = is_valid(x, Val(encoding))

is_valid(x, ::Val) = true

function binarize(x::Vector{T}, d = length(x); binarization = :one_hot) where {T <: Number}
    return collect(Iterators.flatten(map(i -> binarize(i, d; binarization), x)))
end

function binarize(x::T, domain_size::Int; binarization = :one_hot) where {T <: Number}
    return binarize(x, domain(0:(domain_size - 1)), Val(binarization))
end

function binarize(x::T, d::D; binarization=:one_hot
) where {T <: Number, D <: DiscreteDomain{T}}
    return binarize(x, d, Val(binarization))
end

function binarize(x::Vector{T}, d::Vector{D}; binarization=:one_hot
) where {T <: Number, D <: DiscreteDomain{T}}
    return collect(Iterators.flatten(map(elt -> binarize(elt[1], elt[2]; binarization), zip(x,d))))
end

function debinarize(x; binarization = :one_hot)
    ds::Int = if binarization == :domain_wall
        (-1 + sqrt(1 + 4 * length(x))) ÷ 2 + 1
    else
        sqrt(length(x))
    end
    return debinarize(x, ds; binarization)
end

function debinarize(x, domain_size::Int; binarization = :one_hot)
    return debinarize(x, domain(0:(domain_size - 1)); binarization)
end

function debinarize(x, domains_sizes::Vector{Int}; binarization = :one_hot)
    domains = map(ds -> domain(0:(ds - 1)), domains_sizes)
    return debinarize(x, domains; binarization)
end

function debinarize(x, d::D; binarization = :one_hot) where {D <: DiscreteDomain}
    k::Int = length(x) / length(d)
    if binarization == :domain_wall
        typeof(d) <: RangeDomain && first(get_domain(d)) == 0 && (k += 1)
    end
    domains = fill(d, k)
    return debinarize(x, domains; binarization)
end

function debinarize(x, domains::Vector{D}; binarization = :one_hot
) where {D <: DiscreteDomain}
    return debinarize(x, domains, Val(binarization))
end
