function binarize(n, domain_size, ::Val{:one_hot})
	y = zeros(Int, domain_size)
    y[n+1] = 1
    return y
end

function binarize(x::Vector{Int}, dim, ::Val{:one_hot})
    return Iterators.flatten(map(i -> binarize(i, dim), x))
end

function integerize(x, dom_size, ::Val{:one_hot})
    function aux(b)
		a = findall(x -> x == 1, b)
		return length(a) != 1 ? dom_size : first(a) - 1
	end
    return map(i -> aux(@view(x[((i-1)*dom_size+1):(i*dom_size)])), 1:dom_size)
end
