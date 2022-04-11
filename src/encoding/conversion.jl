function binarize(x::Vector{T}, d = length(x); encoding = :one_hot) where {T <: Number}
    return Iterators.flatten(map(i -> binarize(i, d; encoding), x))
end

function binarize(x::T, domain_size::Int; encoding = :one_hot) where {T <: Number}
    return binarize(x, domain(0:(domain_size - 1)), Val(encoding))
end

function binarize(x::T, d::D; encoding=:one_hot) where {T <: Number, D <: DiscreteDomain{T}}
    return binarize(x, d, Val(encoding))
end


function integerize(x; encoding = :one_hot)
    ds::Int = if encoding == :domain_wall
        (-1 + sqrt(1 + 4 * length(x))) ÷ 2 + 1
    else
        sqrt(length(x))
    end
    return integerize(x, ds; encoding)
end

function integerize(x, domain_size; encoding = :one_hot)
    return integerize(x, domain(0:(domain_size - 1)); encoding)
end

function integerize(x, domains_sizes::Vector{Int}; encoding = :one_hot)
    domains = map(ds -> domain(0:(ds - 1)), domains_sizes)
    return integerize(x, domains; encoding)
end

function integerize(x, d::D; encoding = :one_hot) where {D <: DiscreteDomain}
    k::Int = length(x) / length(d)
    if encoding == :domain_wall
        typeof(d) <: RangeDomain && first(get_domain(d)) == 0 && (k += 1)
    end
    domains = fill(d, k)
    return integerize(x, domains; encoding)
end

function integerize(x, domains::Vector{D}; encoding = :one_hot) where {D <: DiscreteDomain}
    return integerize(x, domains, Val(encoding))
end
