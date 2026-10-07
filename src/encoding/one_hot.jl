is_valid(x, ::Val{:one_hot}) = sum(x) == 1

function binarize(x, d::D, ::Val{:one_hot}) where {T <: Number, D <: DiscreteDomain{T}}
    y = zeros(T, length(d))
    is_in = false
    for (i, v) in enumerate(get_domain(d))
        if x == v
            y[i] = 1
            is_in = true
            break
        end
    end
    @assert is_in "The value $x is not in the domain $d"
    return y
end

function debinarize(x, d::D, ::Val{:one_hot}) where {T <: Number, D <: DiscreteDomain{T}}
    valid = false
    val = typemax(T)
    for (b, v) in Iterators.zip(x, get_domain(d))
        if b == 1
            valid && break
            valid = true
            val = v
        end
    end
    return val
end

function debinarize(x, domains, ::Val{:one_hot})
    start = 0
    stop = 0
    function aux(d)
        start = stop + 1
        stop += length(d)
        return debinarize(@view(x[start:stop]), d, Val(:one_hot))
    end
    return map(aux, domains)
end

## SECTION - Test Items
@testitem "One-Hot" tags=[:encoding, :domain_wall] default_imports=false begin
    using ConstraintDomains
    using QUBOConstraints
    using Test

    x = [0, 1, 2, 3, 4]
    b = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]

    @test collect(binarize(x)) == b
    @test debinarize(b) == x

    d = domain(x)
    @test debinarize(collect(binarize(x, d)), d) == x
end
