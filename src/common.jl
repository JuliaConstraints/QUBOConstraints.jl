function δ(X, discrete = true)
    mn, mx = extrema(Iterators.flatten(X))
    return mx - mn + discrete
end

function δ(X, Y, discrete = true)
    mnx, mxx = extrema(Iterators.flatten(X))
    mny, mxy = extrema(Iterators.flatten(Y))
    return max(mxx, mxy) - min(mnx, mny) + discrete
end

to_domains(domain_sizes::Vector{Int}) = map(ds -> domain(0:ds), domain_sizes)

function to_domains(X, ds::Int = δ(X))
    d = domain(0:ds-1)
    return fill(d, length(first(X)))
end

function to_domains(X, d::D) where {D <: DiscreteDomain}
    n::Int = length(first(X)) / domain_size(d)
    return fill(d, n)
end
