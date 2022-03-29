binarize(n, domain_size; encoding = :one_hot) = binarize(n, domain_size, Val(encoding))

function binarize(x::Vector{Int}, dim = length(x); encoding = :one_hot)
    return binarize(x, dim, Val(encoding))
end

function integerize(x; encoding = :one_hot)
    ds = sqrt(length(x))
    if encoding == :domain_wall
        ds = (-1 + sqrt(1 + 4 * length(x))) ÷ 2 + 1
    end
    return integerize(x, Int(ds); encoding)
end

function integerize(x, dom_size; encoding = :one_hot)
    return integerize(x, dom_size, Val(encoding))
end

function integerize(x, domains_sizes::Vector{Int}; encoding = :one_hot)
    i = 1
    y = Vector{Int}()
    for ds in domains_sizes
        push!(y, first(integerize(@view(x[i:i+ds-1]), ds; encoding)))
        i += ds
    end
    return y
end
