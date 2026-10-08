"""Canonical polynomial `offset + Σ linear[i]z[i] + Σ(i<j) quadratic[i,j]z[i]z[j]`.

Indices in constructor terms refer to the supplied bit order. Stored bits are
sorted by identity; use that order for `energy`. Diagonals fold into linear terms.
Default BigInt coefficients avoid overflow. Rational{BigInt} and Float64 are also
supported; floating evaluations receive numerical, not exact, certificates.
Arrays are owned by the component and must be treated as read-only.
"""
struct QUBOComponent{T<:Real}
    bits::Vector{BitID}
    linear::SparseVector{T,Int}
    quadratic::SparseMatrixCSC{T,Int}
    offset::T
    codebooks::Vector{AbstractCodebook}
    auxiliary_meanings::Dict{BitID,String}
    applicability::String
    provenance::String
    function QUBOComponent{T}(bits::AbstractVector{BitID}; linear = [], quadratic = [],
            offset = 0, codebooks = Codebook[],
            auxiliary_meanings = Dict{BitID,String}(), applicability = "",
            provenance = "unspecified") where {T<:Real}
        T in (BigInt, Rational{BigInt}, Float64) || throw(ArgumentError(
            "supported coefficients: BigInt, Rational{BigInt}, Float64"))
        allunique(bits) || throw(ArgumentError("duplicate bit identity"))
        ordered = sort(collect(bits))
        positions = Dict(b => i for (i, b) in enumerate(ordered))
        remap = [positions[b] for b in bits]
        n = length(bits)
        ri, cj, vv = Int[], Int[], T[]
        for (i, v) in linear
            1 <= i <= n || throw(ArgumentError("linear index out of bounds"))
            push!(ri, remap[i]); push!(cj, remap[i]); push!(vv, T(v))
        end
        for (i, j, v) in quadratic
            1 <= i <= n && 1 <= j <= n || throw(ArgumentError("quadratic index out of bounds"))
            a, b = minmax(remap[i], remap[j])
            push!(ri, a); push!(cj, b); push!(vv, T(v))
        end
        all(isfinite, vv) && isfinite(T(offset)) || throw(ArgumentError("nonfinite coefficient"))
        perm = sortperm(eachindex(vv); by = k -> (cj[k], ri[k], vv[k]))
        # Sparse assembly reads the sorted coordinates without retaining them;
        # views avoid three copied permutation buffers while preserving order.
        matrix = dropzeros!(sparse(view(ri, perm), view(cj, perm), view(vv, perm), n, n))
        all(isfinite, nonzeros(matrix)) || throw(ArgumentError("coefficient overflow"))
        li, lv = Int[], T[]
        for j in 1:n, p in nzrange(matrix, j)
            i, v = rowvals(matrix)[p], nonzeros(matrix)[p]
            if i == j
                push!(li, i); push!(lv, v)
            end
        end
        # The canonical matrix already owns sorted, combined sparse storage.
        # Remove its diagonal in place instead of rebuilding off-diagonal terms.
        SparseArrays.fkeep!((i, j, _) -> i != j, matrix)
        books = AbstractCodebook[deepcopy(b) for b in codebooks]
        allunique(b.variable for b in books) || throw(ArgumentError("duplicate codebook variable"))
        covered = BitID[]
        for book in books
            _append_codebook_bits!(covered, book)
        end
        allunique(covered) || throw(ArgumentError("overlapping codebooks"))
        all(b -> haskey(positions, b), covered) || throw(ArgumentError("codebook bit missing from component"))
        meanings = Dict{BitID,String}(b => get(auxiliary_meanings, b, "uninterpreted")
            for b in ordered if b.role !== :primary)
        all(b -> b in keys(meanings), keys(auxiliary_meanings)) ||
            throw(ArgumentError("meaning assigned to an absent or primary bit"))
        return new{T}(ordered, sparsevec(li, lv, n), matrix, T(offset),
            books, meanings, String(applicability), String(provenance))
    end
end

# Specialize iteration once per concrete codebook, instead of carrying an
# abstract inner iterator through the flattened metadata comprehension.
function _append_codebook_bits!(covered::Vector{BitID}, book::AbstractCodebook)
    append!(covered, book.bits)
    return covered
end

function QUBOComponent(bits::AbstractVector{BitID}; coefficient_type::Type{T} = BigInt,
        kwargs...) where {T<:Real}
    return QUBOComponent{T}(bits; kwargs...)
end

function _energy(q::QUBOComponent, z)
    value = q.offset
    for p in eachindex(nonzeros(q.linear))
        z[q.linear.nzind[p]] == 1 && (value += nonzeros(q.linear)[p])
    end
    for j in axes(q.quadratic, 2)
        z[j] == 1 || continue
        for p in nzrange(q.quadratic, j)
            z[rowvals(q.quadratic)[p]] == 1 && (value += nonzeros(q.quadratic)[p])
        end
    end
    return value
end

function energy(q::QUBOComponent, z::AbstractVector{Bool})
    length(z) == length(q.bits) || throw(DimensionMismatch("wrong bit count"))
    return _energy(q, z)
end
function energy(q::QUBOComponent, z::AbstractVector)
    length(z) == length(q.bits) || throw(DimensionMismatch("wrong bit count"))
    all(b -> b == 0 || b == 1, z) || throw(ArgumentError("energy requires binary values"))
    # Validation establishes the same binary predicate used by the sparse
    # traversal; a temporary Bool vector is unnecessary for numeric inputs.
    return _energy(q, z)
end

"""Return `(matrix, offset)` using the upper-triangular `z' * matrix * z` convention."""
function dense_qubo(q::QUBOComponent{T}) where T
    out = fill(zero(T), length(q.bits), length(q.bits))
    for p in eachindex(nonzeros(q.linear))
        i = q.linear.nzind[p]
        out[i, i] = nonzeros(q.linear)[p]
    end
    for j in axes(q.quadratic, 2), p in nzrange(q.quadratic, j)
        out[rowvals(q.quadratic)[p], j] = nonzeros(q.quadratic)[p]
    end
    return (; matrix = out, offset = q.offset)
end

"""Polynomial equality on the same identified bits; does not prove semantic equivalence."""
canonical_equal(a::QUBOComponent, b::QUBOComponent) = a.bits == b.bits &&
    a.offset == b.offset && a.linear == b.linear && a.quadratic == b.quadratic

_linear_terms(q) = [(q.linear.nzind[p], nonzeros(q.linear)[p]) for p in eachindex(nonzeros(q.linear))]
_quadratic_terms(q) = [(rowvals(q.quadratic)[p], j, nonzeros(q.quadratic)[p])
    for j in axes(q.quadratic, 2) for p in nzrange(q.quadratic, j)]

"""Rename every auxiliary into an explicit namespace, preserving primary identities."""
function rename_auxiliaries(q::QUBOComponent{T}, namespace::Symbol) where T
    rename(b) = b.role === :primary ? b :
        BitID(Symbol(namespace, "/", b.owner), b.index; role = b.role)
    bits = rename.(q.bits)
    return QUBOComponent{T}(bits; linear = _linear_terms(q), quadratic = _quadratic_terms(q),
        offset = q.offset, codebooks = q.codebooks,
        auxiliary_meanings = Dict(rename(b) => v for (b, v) in q.auxiliary_meanings),
        applicability = q.applicability, provenance = q.provenance)
end

"""Sum polynomials, sharing compatible primary codebooks and renaming right auxiliaries.

The result carries no inherited exactness certificate. Shared auxiliary witnesses
are deliberately unsupported in Q1; certify the resulting component explicitly.
"""
function compose(a::QUBOComponent{T}, original::QUBOComponent{S};
        namespace::Symbol = :right) where {T,S}
    b = rename_auxiliaries(original, namespace)
    collision = intersect(filter(x -> x.role !== :primary, a.bits), b.bits)
    isempty(collision) || throw(ArgumentError("auxiliary namespace collision"))
    books = copy(a.codebooks)
    for book in b.codebooks
        old = findfirst(x -> x.variable == book.variable, books)
        if old === nothing
            push!(books, book)
        else
            _compatible(books[old], book) || throw(ArgumentError("incompatible shared codebooks"))
        end
    end
    bits = sort!(unique(vcat(a.bits, b.bits)))
    pos = Dict(bit => i for (i, bit) in enumerate(bits))
    linear = [(pos[q.bits[i]], v) for q in (a, b) for (i, v) in _linear_terms(q)]
    quadratic = [(pos[q.bits[i]], pos[q.bits[j]], v)
        for q in (a, b) for (i, j, v) in _quadratic_terms(q)]
    R = T === Float64 || S === Float64 ? Float64 : promote_type(T, S)
    return QUBOComponent{R}(bits; linear, quadratic, offset = a.offset + b.offset,
        codebooks = books,
        auxiliary_meanings = merge(a.auxiliary_meanings, b.auxiliary_meanings),
        applicability = "(" * a.applicability * ") AND (" * b.applicability * ")",
        provenance = "add(" * a.provenance * ", " * b.provenance * ")")
end
