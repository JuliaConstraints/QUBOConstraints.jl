function is_valid(x, ::Val{:domain_wall})
    for (i, b) in enumerate(x)
        iszero(b) && (return iszero(x[i+1:end]))
    end
    return true
end

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

function debinarize(x, d::D, ::Val{:domain_wall}) where {T <: Number, D <: DiscreteDomain{T}}
    b0::Bool = typeof(d) <: RangeDomain && first(get_domain(d)) == 0
    at_one = true
    val = typemax(T)
    for (i, v) in zip((b0 ? 0 : 1):length(x), get_domain(d))
        b = i > 0 ? x[i] : 1
        if at_one
            b == 1 ? (val = v) : (at_one = false)
        else
            b == 1 && return typemax(T)
        end
    end
    return val
end

function debinarize(x, domains::Vector{D}, ::Val{:domain_wall}
) where {T <: Number, D <: DiscreteDomain{T}}
    start = 0
    stop = 0
    function aux(d)
        b0 = typeof(d) <: RangeDomain && first(get_domain(d)) == 0
        start = stop + 1
        stop += length(d) - b0

        return debinarize(x[start:stop], d, Val(:domain_wall))
    end
    return map(aux, domains)
end

## SECTION - Test Items
@testitem "Domain Wall" tags = [:encoding, :domain_wall] default_imports=false begin
    using ConstraintDomains
    using QUBOConstraints
    using Test

    x = [0,1,2,3,4,5]

    d = domain([0,1,2,3,4])
    y = [0,1,2,3,4]

    a = [0,0,0,0,0, 1,0,0,0,0, 1,1,0,0,0, 1,1,1,0,0, 1,1,1,1,0, 1,1,1,1,1]

    @test collect(binarize(x; binarization = :domain_wall)) == a
    @test debinarize(a; binarization = :domain_wall) == x

    @test debinarize(collect(binarize(y, d; binarization = :domain_wall)), d; binarization = :domain_wall) == y
end
