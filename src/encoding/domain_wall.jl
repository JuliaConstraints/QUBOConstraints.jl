function binarize(n, domain_size, ::Val{:domain_wall})
	y = zeros(Int, domain_size)
    for i in 1:n
        y[i] = 1
    end
    return y
end

function binarize(x::Vector{Int}, dim, ::Val{:domain_wall})
    return Iterators.flatten(map(i -> binarize(i, dim, Val(:domain_wall)), x))
end

function integerize(x, dom_size, ::Val{:domain_wall})
    function aux(b)
        σ = sum(b)
        if σ == dom_size || σ ≥ findfirst(x -> x == 0, b) - 1
            return σ
        end
        return dom_size + 1
	end
    n = length(x) ÷ dom_size
    return map(i -> aux(@view(x[((i-1)*dom_size+1):(i*dom_size)])), 1:n)
end
