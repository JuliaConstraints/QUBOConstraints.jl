const _ENCODINGS = (:native,:one_hot,:zero_one_hot,:domain_wall,:unary,
    :binary,:bounded_binary,:gray,:arithmetic,:bounded_coefficient,:fibonacci)

"""Compact finite-domain encoding. Stores O(domain + bits), never all redundant codes.

Use `structured_codebook` to construct. Indices 0:length(values)-1 are encoded, then
mapped to semantic values. Stored arrays are owned and must be treated as read-only.
"""
struct StructuredCodebook{V} <: AbstractCodebook
    variable::Symbol
    bits::Vector{BitID}
    values::Vector{V}
    encoding::Symbol
    weights::Vector{BigInt}
    coefficient_bound::BigInt
    function StructuredCodebook(variable::Symbol,domain::Vector{V},encoding::Symbol,mu::Integer) where V
        isempty(domain) && throw(ArgumentError("empty semantic domain"))
        allunique(domain) || throw(ArgumentError("duplicate semantic values"))
        encoding in _ENCODINGS || throw(ArgumentError("unsupported structured encoding"))
        mu>0 || throw(ArgumentError("coefficient bound must be positive"))
        K = length(domain)-1
        weights = BigInt[]
        if encoding===:native
            K<=1 || throw(ArgumentError("native encoding needs at most two values"))
            push!(weights,1)
        elseif encoding===:one_hot
            weights = BigInt.(0:K)
        elseif encoding===:zero_one_hot
            weights = BigInt.(1:K)
        elseif encoding in (:domain_wall,:unary)
            weights = ones(BigInt,K)
        elseif encoding in (:binary,:gray)
            weights = BigInt[big(1)<<i for i in 0:ndigits(K;base=2)-1 if K>0]
        elseif encoding===:fibonacci
            a,b = big(1),big(2)
            while a<=K
                push!(weights,a)
                a,b = b,a+b
            end
        else
            total = big(0)
            while total<K
                next = encoding===:arithmetic ? big(length(weights)+1) : total+1
                encoding===:bounded_coefficient && (next=min(next,big(mu)))
                push!(weights,min(next,K-total))
                total += last(weights)
            end
            # Sorting the final residual preserves subset sums and greedy completeness.
            sort!(weights)
        end
        bits = [BitID(variable,i) for i in eachindex(weights)]
        return new{V}(variable,bits,copy(domain),encoding,weights,big(mu))
    end
end

"""E00-E10 registry v1. E02 uses the FIRST domain value as the all-zero code.

`coefficient_bound` is used by E09 only; it bounds decoding weights, not squared-QUBO
coefficients. Signed offsets/arbitrary dictionaries are represented by `values`.
"""
structured_codebook(variable::Symbol,values;encoding::Symbol=:one_hot,coefficient_bound::Integer=2) =
    StructuredCodebook(variable,collect(values),encoding,coefficient_bound)

function _value_index(book::StructuredCodebook,z)
    length(z)==length(book.bits) && all(b->b==0 || b==1,z) || return 0
    mode = book.encoding
    if mode===:one_hot
        count(!iszero,z)==1 || return 0
    elseif mode===:zero_one_hot
        count(!iszero,z)<=1 || return 0
    elseif mode===:domain_wall
        any(!iszero(z[i+1]) && iszero(z[i]) for i in 1:length(z)-1) && return 0
    elseif mode===:fibonacci
        any(!iszero(z[i]) && !iszero(z[i+1]) for i in 1:length(z)-1) && return 0
    end
    rank = big(0)
    if mode===:gray
        binary = false
        for i in reverse(eachindex(z))
            binary = xor(binary,!iszero(z[i]))
            binary && (rank+=book.weights[i])
        end
    else
        for i in eachindex(z)
            !iszero(z[i]) && (rank+=book.weights[i])
        end
    end
    return rank<length(book.values) ? Int(rank)+1 : 0
end

function decode_code(book::StructuredCodebook{V},z) where V
    i = _value_index(book,z)
    return i==0 ? DecodeResult{V}(false,nothing,:invalid_code) :
        DecodeResult{V}(true,book.values[i],:valid)
end

function encode_code(book::StructuredCodebook,value)
    i = findfirst(v->isequal(v,value),book.values)
    i===nothing && throw(ArgumentError("value outside codebook domain"))
    rank = big(i-1)
    z = falses(length(book.bits))
    mode = book.encoding
    if mode===:one_hot
        z[i] = true
    elseif mode===:zero_one_hot
        i>1 && (z[i-1]=true)
    elseif mode===:domain_wall
        z[1:i-1] .= true
    elseif mode===:gray
        rank = xor(rank,rank>>1)
        for j in eachindex(z)
            z[j] = !iszero(rank & book.weights[j])
        end
    else
        for j in reverse(eachindex(z))
            if book.weights[j]<=rank
                z[j] = true
                rank -= book.weights[j]
            end
        end
        iszero(rank) || error("incomplete encoding weights")
    end
    return z
end

"""Materialize all representations of ONE value under an explicit full-cube budget.

Large redundant books remain usable for encode/decode/validity without this operation.
"""
function encoded_codes(book::StructuredCodebook,value;max_states::Integer=1_048_576)
    i = findfirst(v->isequal(v,value),book.values)
    i===nothing && throw(ArgumentError("value outside codebook domain"))
    width = length(book.bits)
    width<8sizeof(Int)-1 && big(2)^width<=max_states || throw(ArgumentError("representation enumeration exceeds max_states"))
    found = Vector{Bool}[]
    for mask in 0:(1<<width)-1
        z = [!iszero(mask & (1<<(j-1))) for j in 1:width]
        _value_index(book,z)==i && push!(found,z)
    end
    return isempty(found) ? falses(width,0) : hcat(found...)
end

_compatible(a::StructuredCodebook,b::StructuredCodebook) = a.variable==b.variable &&
    a.bits==b.bits && isequal(a.values,b.values) && a.encoding==b.encoding &&
    a.weights==b.weights && a.coefficient_bound==b.coefficient_bound

"""Versioned inventory. Planned entries remain explicit and are never silently aliased."""
function encoding_registry()
    entries = [(id="E"*lpad(i-1,2,'0'),encoding=e,status=:implemented) for (i,e) in enumerate(_ENCODINGS)]
    append!(entries,[(id="E"*string(i),encoding=e,status=:planned) for (i,e) in
        ((11,:mixed_radix),(12,:digit_binary),(13,:digit_one_hot_domain_wall),(14,:semi_domain),
         (16,:sign_magnitude_twos_complement),(17,:signed_digit))])
    push!(entries,(id="E15",encoding=:signed_offset,status=:implemented))
    return (;schema_version="qubo-encodings/1",entries=sort!(entries;by=x->x.id))
end
