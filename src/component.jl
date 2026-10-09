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
        ordered = sort!(collect(bits))
        # Equal immutable identities are adjacent in the canonical order.
        # Check the owned sorted buffer instead of allocating a hash set.
        all(i -> ordered[i-1] != ordered[i], 2:length(ordered)) ||
            throw(ArgumentError("duplicate bit identity"))
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
        # Canonical upper-triangular columns place a stored diagonal last.
        # Its sorted unique indices can use compact owned buffers directly.
        diagonal_count = 0
        for j in 1:n
            p = matrix.colptr[j+1] - 1
            diagonal_count += p >= matrix.colptr[j] && rowvals(matrix)[p] == j
        end
        li = Vector{Int}(undef, diagonal_count)
        lv = Vector{T}(undef, diagonal_count)
        diagonal_index = 0
        for j in 1:n
            p = matrix.colptr[j+1] - 1
            if p >= matrix.colptr[j] && rowvals(matrix)[p] == j
                diagonal_index += 1
                li[diagonal_index] = j
                lv[diagonal_index] = nonzeros(matrix)[p]
            end
        end
        # The canonical matrix already owns sorted, combined sparse storage.
        # Remove its diagonal in place instead of rebuilding off-diagonal terms.
        SparseArrays.fkeep!((i, j, _) -> i != j, matrix)
        books = AbstractCodebook[deepcopy(b) for b in codebooks]
        # Iterate larger owned book lists directly, retaining the original
        # key equality and support for heterogeneous variable key types.
if length(books) < 64
            allunique(b.variable for b in books) ||
                throw(ArgumentError("duplicate codebook variable"))
        else
            seen_variables = Set{Any}()
            sizehint!(seen_variables, length(books))
            for book in books
                in!(book.variable, seen_variables) &&
                    throw(ArgumentError("duplicate codebook variable"))
            end
        end
        covered = BitID[]
        for book in books
            _append_codebook_bits!(covered, book)
        end
        # Canonical positions bijectively identify the covered bits. Larger
        # valid lists can retain integer keys; small or missing-bit lists keep
        # the original overlap-before-missing validation.
        if length(covered) >= 128 && all(b -> haskey(positions, b), covered)
            seen = Set{Int}()
            sizehint!(seen, min(length(covered), length(positions)))
            for b in covered
                position = positions[b]
                position in seen && throw(ArgumentError("overlapping codebooks"))
                push!(seen, position)
            end
        else
            allunique(covered) || throw(ArgumentError("overlapping codebooks"))
            all(b -> haskey(positions, b), covered) ||
                throw(ArgumentError("codebook bit missing from component"))
        end
        meanings = Dict{BitID,String}(b => get(auxiliary_meanings, b, "uninterpreted")
            for b in ordered if b.role !== :primary)
        all(b -> b in keys(meanings), keys(auxiliary_meanings)) ||
            throw(ArgumentError("meaning assigned to an absent or primary bit"))
        return new{T}(ordered, SparseVector(n, li, lv), matrix, T(offset),
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
    R = T === Float64 || S === Float64 ? Float64 : promote_type(T, S)
    linear = Tuple{Int,R}[]
    quadratic = Tuple{Int,Int,R}[]
    sizehint!(linear, nnz(a.linear) + nnz(b.linear))
    sizehint!(quadratic, nnz(a.quadratic) + nnz(b.quadratic))
    _append_composed_terms!(linear, quadratic, pos, a)
    _append_composed_terms!(linear, quadratic, pos, b)
    return QUBOComponent{R}(bits; linear, quadratic, offset = a.offset + b.offset,
        codebooks = books,
        auxiliary_meanings = merge(a.auxiliary_meanings, b.auxiliary_meanings),
        applicability = "(" * a.applicability * ") AND (" * b.applicability * ")",
        provenance = "add(" * a.provenance * ", " * b.provenance * ")")
end

# Keep sparse traversal specialized for each input coefficient type while
# converting into the final promoted term buffers.
function _append_composed_terms!(linear::Vector{Tuple{Int,R}},
        quadratic::Vector{Tuple{Int,Int,R}}, positions, q::QUBOComponent) where R
    for p in eachindex(nonzeros(q.linear))
        i = q.linear.nzind[p]
        push!(linear, (positions[q.bits[i]], nonzeros(q.linear)[p]))
    end
    for j in axes(q.quadratic, 2), p in nzrange(q.quadratic, j)
        push!(quadratic, (positions[q.bits[rowvals(q.quadratic)[p]]],
            positions[q.bits[j]], nonzeros(q.quadratic)[p]))
    end
    return nothing
end
