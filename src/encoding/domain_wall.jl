function binarize(x, d::D, ::Val{:domain_wall}) where {T <: Number, D <: DiscreteDomain{T}}
    b0 = typeof(d) <: RangeDomain && first(get_domain(d)) == 0
    y = zeros(T, length(d) - b0)
    is_in = false
    for (i, v) in enumerate(get_domain(d))
        i > b0 && (y[i - b0] = 1)
        (is_in = x == v) && break
    end
    @assert is_in "The value $x is not in the domain $d"
    return y
end

function integerize(x, d::D, ::Val{:domain_wall}) where {T <: Number, D <: DiscreteDomain{T}}
    b0 = typeof(d) <: RangeDomain && first(get_domain(d)) == 0
    itr = Iterators.flatten(b0 ? (1, x) : (x))
    at_one = true
    val = typemax(T)
    @info "Debug: start" b0 itr at_one val d collect(itr)
    for (b, v) in Iterators.zip(itr, get_domain(d))
        @info "Debug: for" at_one val b v
        if at_one
            b == 1 && (val = v)
            b == 0 && (at_one = false)
            @info "Debug: if" at_one val b v
        else

            @info "Debug: else" at_one val b v
            b == 1 && return typemax(T)
        end
    end
    return val
end

function integerize(x, domains, ::Val{:domain_wall})
    start = 0
    stop = 0
    function aux(d)
        b0 = typeof(d) <: RangeDomain && first(get_domain(d)) == 0
        start = stop + 1
        stop += length(d) - b0
        return integerize(@view(x[start:stop]), d, Val(:domain_wall))
    end
    return map(aux, domains)
end
