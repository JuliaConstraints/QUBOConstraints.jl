binarize(n, domain_size; encoding = :one_hot) = binarize(n, domain_size, Val(encoding))

binarize(x, dim = length(x); encoding = :one_hot) = binarize(x, dim, Val(encoding))

function integerize(x, dom_size = Int(sqrt(length(x))); encoding = :one_hot)
    return integerize(x, dom_size, Val(encoding))
end
