"""Integer condition with a scalar or finite-set/interval operand (no implicit variable lookup)."""
struct IntegerCondition{T}
    operator::Symbol
    operand::T
    function IntegerCondition(operator::Symbol,operand)
        operator in (:eq,:ne,:lt,:le,:ge,:gt,:in,:notin) || throw(ArgumentError("unknown integer condition"))
        if operator in (:in,:notin)
            operand isa Union{AbstractVector,AbstractSet,AbstractRange} && eltype(operand)<:Integer ||
                throw(ArgumentError("integer set or interval operand required"))
            operand = operand isa AbstractRange ? operand : sort!(unique(BigInt.(collect(operand))))
        else
            operand isa Integer || throw(ArgumentError("integer scalar operand required"))
            operand = big(operand)
        end
        return new{typeof(operand)}(operator,operand)
    end
end
function satisfies(c::IntegerCondition,x::Integer)
    op,y = c.operator,c.operand
    op===:eq && return x==y
    op===:ne && return x!=y
    op===:lt && return x<y
    op===:le && return x<=y
    op===:ge && return x>=y
    op===:gt && return x>y
    return (x in y)==(op===:in)
end

mutable struct _ArithmeticBuilder
    bits::Vector{BitID}
    linear::Vector{Tuple{Int,BigInt}}
    quadratic::Vector{Tuple{Int,Int,BigInt}}
    offset::BigInt
    meanings::Dict{BitID,String}
    max_terms::Int
end
_ArithmeticBuilder(bits,max_terms) = _ArithmeticBuilder(copy(bits),Tuple{Int,BigInt}[],Tuple{Int,Int,BigInt}[],big(0),Dict{BitID,String}(),max_terms)
function _arithmetic_budget(b,extra=0)
    requested = big(length(b.linear))+length(b.quadratic)+extra
    requested<=b.max_terms || throw(CompilationLimit(:arithmetic_terms,requested,big(b.max_terms)))
end
function _arithmetic_add!(b,q;semantic_owner=nothing)
    identities = [semantic_owner===nothing || bit.role!==:primary ? bit :
        BitID(semantic_owner,bit.index;role=:semantic_auxiliary) for bit in q.bits]
    _arithmetic_budget(b,nnz(q.linear)+nnz(q.quadratic))
    pos = Dict(bit=>i for (i,bit) in enumerate(b.bits))
    for (old,bit) in zip(q.bits,identities)
        if !haskey(pos,bit)
            push!(b.bits,bit); pos[bit]=length(b.bits)
        end
        if bit.role!==:primary
            b.meanings[bit] = old.role===:primary ? "encoded semantic witness for $(semantic_owner), bit $(bit.index)" : q.auxiliary_meanings[old]
        end
    end
    ids = [pos[bit] for bit in identities]
    append!(b.linear,((ids[i],v) for (i,v) in _linear_terms(q)))
    append!(b.quadratic,((ids[i],ids[j],v) for (i,j,v) in _quadratic_terms(q)))
    b.offset += q.offset
    return Dict(old=>ids[i] for (i,old) in enumerate(q.bits))
end
function _arithmetic_square!(b,offset,terms)
    combined = Dict{Int,BigInt}()
    for (i,w) in terms
        combined[i] = get(combined,i,big(0))+w
    end
    filter!(p->!iszero(last(p)),combined)
    rows = sort!(collect(combined);by=first)
    _arithmetic_budget(b,big(length(rows))*(length(rows)+1)÷2)
    b.offset += offset^2
    append!(b.linear,((i,w*w+2offset*w) for (i,w) in rows))
    append!(b.quadratic,((rows[i][1],rows[j][1],2rows[i][2]*rows[j][2]) for i in eachindex(rows) for j in i+1:length(rows)))
end

# Affine semantic value, not merely the rank. Compact encodings require an arithmetic dictionary.
function _arithmetic_affine(book::StructuredCodebook,values,ids)
    all(v->v isa Integer,values) || throw(ArgumentError("integer decoded values required"))
    length(values)==length(book.values) || throw(DimensionMismatch("value map size"))
    v = BigInt.(values)
    mode = book.encoding
    mode===:one_hot && return (big(0),[(ids[book.bits[k]],v[k]) for k in eachindex(v)])
    if mode===:domain_wall
        return (first(v),[(ids[book.bits[k]],v[k+1]-v[k]) for k in eachindex(book.bits)])
    elseif mode===:zero_one_hot
        return (first(v),[(ids[book.bits[k]],v[k+1]-first(v)) for k in eachindex(book.bits)])
    end
    step = length(v)==1 ? big(0) : v[2]-v[1]
    all(v[i]==v[1]+(i-1)*step for i in eachindex(v)) || throw(ArgumentError("non-affine compact value map; use one_hot/domain_wall/zero_one_hot or truth_component"))
    mode===:gray && length(book.bits)>1 && !iszero(step) && throw(ArgumentError("Gray decoding is not affine; use truth_component"))
    return (first(v),[(ids[bit],step*book.weights[k]) for (k,bit) in enumerate(book.bits)])
end

"""Compile an integer separable sum subject to a condition, with exact nonnegative witnesses.

`value_maps[i][k]` is the contribution of semantic value k of codebook i. Defaults to
the semantic integer itself. Repeated source occurrences can be combined in coefficients.
No minimization over primary representations. Slack encoding is independently selected.
Coefficients are constants; variable coefficients require the explicit finite truth fallback.
The result preserves zeros with unit gap, not necessarily a 0/1 violation indicator.
"""
function sum_component(books,condition::IntegerCondition;coefficients=ones(Int,length(books)),
        value_maps=nothing,slack_encoding::Symbol=:bounded_binary,coefficient_bound::Integer=2,
        max_witness_values::Integer=100000,max_terms::Integer=1000000,
        namespace::Symbol=:sum_witness)
    books = collect(books)
    all(b->b isa StructuredCodebook,books) || throw(ArgumentError("structured codebooks required"))
    allunique(b.variable for b in books) || throw(ArgumentError("duplicate variable identity"))
    max_witness_values>0 && 0<max_terms<=typemax(Int) || throw(ArgumentError("positive finite budgets required"))
    coefficients = collect(coefficients)
    length(coefficients)==length(books) && all(c->c isa Integer,coefficients) || throw(ArgumentError("one integer coefficient per codebook required"))
    maps = value_maps===nothing ? [book.values for book in books] : collect(value_maps)
    length(maps)==length(books) || throw(DimensionMismatch("one value map per codebook"))
    bits = reduce(vcat,(book.bits for book in books);init=BitID[])
    builder = _ArithmeticBuilder(bits,Int(max_terms))
    ids = Dict(bit=>i for (i,bit) in enumerate(bits))
    offset,lower,upper = big(0),big(0),big(0)
    terms = Tuple{Int,BigInt}[]
    for (book,map,c) in zip(books,maps,coefficients)
        length(map)==length(book.values) && all(v->v isa Integer,map) || throw(ArgumentError("integer map with one contribution per value required"))
        contributions = [big(c)*big(v) for v in map]
        lower += minimum(contributions); upper += maximum(contributions)
        constant,rows = _arithmetic_affine(book,contributions,ids)
        offset += constant; append!(terms,rows)
        _arithmetic_add!(builder,encoding_validity(book;max_terms))
    end
    # Bounds over independent domains are conservative even when a higher-level scope repeats.
    op,rhs = condition.operator,condition.operand
    lo,hi = lower,upper
    op===:eq && (lo=max(lo,rhs); hi=min(hi,rhs))
    op===:le && (hi=min(hi,rhs))
    op===:lt && (hi=min(hi,rhs-1))
    op===:ge && (lo=max(lo,rhs))
    op===:gt && (lo=max(lo,rhs+1))
    if op===:in && rhs isa AbstractUnitRange
        isempty(rhs) ? (hi=lo-1) : (lo=max(lo,first(rhs)); hi=min(hi,last(rhs)))
    end
    width = max(big(0),hi-lo+1)
    width<=max_witness_values || throw(CompilationLimit(:witness_values,width,big(max_witness_values)))
    allowed = [v for v in lo:hi if satisfies(condition,v)]
    if isempty(allowed)
        builder.offset += 1
    elseif length(allowed)==1
        _arithmetic_square!(builder,offset-only(allowed),terms)
    else
        any(book->book.variable===namespace,books) && throw(ArgumentError("witness namespace collides with a primary variable"))
        witness = structured_codebook(namespace,allowed;encoding=slack_encoding,coefficient_bound)
        witness_ids = _arithmetic_add!(builder,encoding_validity(witness;max_terms);semantic_owner=namespace)
        constant,rows = _arithmetic_affine(witness,allowed,witness_ids)
        _arithmetic_square!(builder,offset-constant,[terms;[(i,-w) for (i,w) in rows]])
    end
    return QUBOComponent(builder.bits;linear=builder.linear,quadratic=builder.quadratic,offset=builder.offset,
        codebooks=books,auxiliary_meanings=builder.meanings,
        applicability="integer separable sum; exact zero set; all valid primary representations; private result witnesses",
        provenance="Q5 sum v1; condition=$(condition); coefficients=$(coefficients); maps=$(maps); witness_encoding=$(slack_encoding); bounds=[$(lower),$(upper)]")
end

"""Count membership in a constant value set, then apply an integer condition."""
function count_component(books,values,condition::IntegerCondition;kwargs...)
    maps = [[Int(any(v->isequal(x,v),values)) for x in book.values] for book in books]
    return sum_component(books,condition;value_maps=maps,kwargs...)
end
