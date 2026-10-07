"""Stable logical bit identity; owner namespaces auxiliaries and names variables."""
struct BitID
    owner::Symbol
    index::Int
    role::Symbol
    function BitID(owner::Symbol, index::Integer; role::Symbol = :primary)
        index > 0 || throw(ArgumentError("bit index must be positive"))
        role in (:primary, :semantic_auxiliary, :quadratization_auxiliary) ||
            throw(ArgumentError("unknown logical bit role: $role"))
        return new(owner, Int(index), role)
    end
end
Base.:(==)(a::BitID, b::BitID) = a.owner == b.owner && a.index == b.index && a.role == b.role
Base.isequal(a::BitID, b::BitID) = a == b
Base.hash(a::BitID, h::UInt) = hash((a.owner, a.index, a.role), h)
Base.isless(a::BitID, b::BitID) =
    isless((a.role, a.owner, a.index), (b.role, b.owner, b.index))

abstract type AbstractCodebook end

# Typed protocol accessors retain field contracts across heterogeneous codebook kinds.
# In particular, oracle position and bit buffers must not become Vector{Any}.
_book_bits(book::AbstractCodebook)::Vector{BitID} = book.bits
_book_values(book::AbstractCodebook)::Vector = book.values

"""Finite codebook. Each column is a valid code mapped to one value index.

Multiple columns may encode the same value; a code may never denote two values.
Inputs are copied. Treat stored arrays as read-only.
"""
struct Codebook{V} <: AbstractCodebook
    variable::Symbol
    bits::Vector{BitID}
    values::Vector{V}
    codes::Matrix{Bool}
    value_indices::Vector{Int}
    function Codebook(variable::Symbol, bits::AbstractVector{BitID}, values::AbstractVector{V},
            codes::AbstractMatrix, value_indices::AbstractVector{<:Integer}) where V
        isempty(values) && throw(ArgumentError("empty semantic domain"))
        allunique(values) || throw(ArgumentError("duplicate semantic values"))
        allunique(bits) || throw(ArgumentError("duplicate bits"))
        all(b -> b.role === :primary, bits) || throw(ArgumentError("codebook bits must be primary"))
        size(codes) == (length(bits), length(value_indices)) ||
            throw(DimensionMismatch("codebook dimensions do not match"))
        all(x -> x == 0 || x == 1, codes) || throw(ArgumentError("codes must be binary"))
        all(i -> 1 <= i <= length(values), value_indices) || throw(ArgumentError("value index out of range"))
        Set(value_indices) == Set(eachindex(values)) || throw(ArgumentError("unrepresented semantic value"))
        allunique(Tuple(col) for col in eachcol(codes)) || throw(ArgumentError("duplicate code"))
        return new{V}(variable, collect(bits), collect(values), Matrix{Bool}(codes), Int.(value_indices))
    end
end

"""Explicit decoding result. `valid` distinguishes `nothing` values from invalid codes."""
struct DecodeResult{V}
    valid::Bool
    value::Union{Nothing,V}
    reason::Symbol
end

function _code_index(book::Codebook, z)
    length(z) == length(book.bits) || return 0
    all(b -> b == 0 || b == 1, z) || return 0
    for j in axes(book.codes, 2)
        all(i -> book.codes[i, j] == z[i], axes(book.codes, 1)) && return j
    end
    return 0
end

function decode_code(book::Codebook{V}, z) where V
    j = _code_index(book, z)
    j == 0 && return DecodeResult{V}(false, nothing, :invalid_code)
    return DecodeResult{V}(true, book.values[book.value_indices[j]], :valid)
end

"""Return the first registered valid code for a value; use `encoded_codes` for all codes."""
function encode_code(book::Codebook, value)
    i = findfirst(v -> isequal(v, value), book.values)
    i === nothing && throw(ArgumentError("value outside codebook domain"))
    return book.codes[:, findfirst(==(i), book.value_indices)]
end

function encoded_codes(book::Codebook, value)
    i = findfirst(v -> isequal(v, value), book.values)
    i === nothing && throw(ArgumentError("value outside codebook domain"))
    return book.codes[:, findall(==(i), book.value_indices)]
end

"""Reference explicit codebook for one-hot or domain-wall over an arbitrary finite domain."""
function codebook(variable::Symbol, values; encoding::Symbol = :one_hot)
    domain = collect(values)
    n = length(domain)
    n > 0 || throw(ArgumentError("empty semantic domain"))
    encoding in (:one_hot, :domain_wall) || throw(ArgumentError("unsupported codebook encoding"))
    width = encoding === :one_hot ? n : n - 1
    codes = falses(width, n)
    for j in 1:n, i in 1:width
        codes[i, j] = encoding === :one_hot ? i == j : i < j
    end
    return Codebook(variable, [BitID(variable, i) for i in 1:width], domain, codes, collect(1:n))
end

function _compatible(a::Codebook, b::Codebook)
    return a.variable == b.variable && a.bits == b.bits && isequal(a.values, b.values) &&
        a.codes == b.codes && a.value_indices == b.value_indices
end

_compatible(a::AbstractCodebook, b::AbstractCodebook) = false
_value_index(book::Codebook,z) = ((j=_code_index(book,z)); j==0 ? 0 : book.value_indices[j])
