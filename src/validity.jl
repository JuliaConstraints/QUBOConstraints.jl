"""Unit-gap one-hot validity polynomial; rejects codebooks that are not canonical one-hot."""
function one_hot_validity(books::AbstractVector{<:AbstractCodebook})
    for book in books
        n = length(book.values)
        if book isa StructuredCodebook
            book.encoding===:one_hot || throw(ArgumentError("canonical one-hot codebook required"))
            continue
        end
        length(book.bits)==n && size(book.codes,2)==n && book.value_indices==collect(1:n) ||
            throw(ArgumentError("canonical one-hot codebook required"))
        all(book.codes[i,j] == (i==j) for i in 1:n for j in 1:n) ||
            throw(ArgumentError("canonical one-hot codebook required"))
    end
    bits = reduce(vcat,(b.bits for b in books);init=BitID[])
    positions = Dict(b=>i for (i,b) in enumerate(bits))
    linear = [(i,big(-1)) for i in eachindex(bits)]
    quadratic = [(positions[book.bits[i]],positions[book.bits[j]],big(2))
        for book in books for i in eachindex(book.bits) for j in i+1:length(book.bits)]
    return QUBOComponent{BigInt}(bits;linear,quadratic,offset=big(length(books)),codebooks=books,
        provenance="one-hot validity sum of squares v1")
end

"""Add a sufficient one-hot validity weight, preserving energies on every valid code.

For exact coefficients, B = offset + sum(min(0,coefficient)) bounds every binary
energy below. Weight max(0,gap-B) makes every invalid code's energy at least gap.
This proves validity ONLY: it does not establish the constraint's zero set on valid codes.
"""
function guard_one_hot(q::QUBOComponent{T}; gap::Integer=1) where T
    T in (BigInt,Rational{BigInt}) || throw(ArgumentError("exact coefficients required"))
    gap>0 || throw(ArgumentError("positive integer gap required"))
    covered = BitID[b for book in q.codebooks for b in book.bits]
    Set(covered)==Set(filter(b->b.role===:primary,q.bits)) || throw(ArgumentError("incomplete primary codebooks"))
    validity = one_hot_validity(q.codebooks)
    lower = q.offset + sum(v->min(zero(T),v),nonzeros(q.linear);init=zero(T)) +
        sum(v->min(zero(T),v),nonzeros(q.quadratic);init=zero(T))
    weight = max(zero(T),T(gap)-lower)
    weighted = QUBOComponent{T}(validity.bits;
        linear=[(i,weight*v) for (i,v) in _linear_terms(validity)],
        quadratic=[(i,j,weight*v) for (i,j,v) in _quadratic_terms(validity)],
        offset=weight*validity.offset,codebooks=validity.codebooks,
        provenance="validity bound v1; B=$(lower); gap=$(gap); weight=$(weight)")
    return (;component=compose(q,weighted),lower_bound=lower,weight)
end
