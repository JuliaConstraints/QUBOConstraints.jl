"""Allowed unary rows in the order of the referenced codebook's semantic values."""
struct UnaryRelation
    variable::Symbol
    allowed::BitVector
    UnaryRelation(variable::Symbol,allowed::AbstractVector{Bool}) = new(variable,BitVector(allowed))
end

"""Allowed binary table: rows reference left values, columns reference right values.

Repeated scope variables are allowed; only diagonal assignments are then reachable.
Stored tables are copied and must be treated as read-only.
"""
struct BinaryRelation
    left::Symbol
    right::Symbol
    allowed::BitMatrix
    BinaryRelation(left::Symbol,right::Symbol,allowed::AbstractMatrix{Bool}) = new(left,right,BitMatrix(allowed))
end

"""Explicit short-table wildcard. Ordinary strings such as \"*\" remain ordinary values."""
struct AnyValue end
const ANY_VALUE = AnyValue()

# Affine value indicators, valid on the specified codebook. Signed outside it.
function _value_indicator(book::StructuredCodebook,k,ids)
    n = length(book.values)
    if book.encoding===:one_hot
        return (big(0),[(ids[k],big(1))])
    elseif book.encoding===:zero_one_hot
        return k==1 ? (big(1),[(i,big(-1)) for i in ids]) : (big(0),[(ids[k-1],big(1))])
    elseif book.encoding===:domain_wall
        n==1 && return (big(1),Tuple{Int,BigInt}[])
        k==1 && return (big(1),[(first(ids),big(-1))])
        k==n && return (big(0),[(last(ids),big(1))])
        return (big(0),[(ids[k-1],big(1)),(ids[k],big(-1))])
    elseif book.encoding===:native
        return k==1 ? (big(1),[(first(ids),big(-1))]) : (big(0),[(first(ids),big(1))])
    end
    throw(ArgumentError("affine value indicators unsupported for $(book.encoding)"))
end

"""Validity for canonical explicit Q1 one-hot/domain-wall books; arbitrary tables are rejected."""
function encoding_validity(book::Codebook;max_terms::Integer=1_000_000)
    for mode in (:one_hot,:domain_wall)
        reference = codebook(book.variable,book.values;encoding=mode)
        if _compatible(book,reference)
            raw = encoding_validity(structured_codebook(book.variable,book.values;encoding=mode);max_terms)
            return QUBOComponent(raw.bits;linear=_linear_terms(raw),quadratic=_quadratic_terms(raw),
                offset=raw.offset,codebooks=[book],provenance="canonical explicit validity v1; $(mode)")
        end
    end
    throw(ArgumentError("validity construction requires a supported canonical codebook"))
end

"""Protect all invalid primary codes using independent validity witnesses and an exact bound.

For B=offset+sum(min(0,c)), weight=max(0,gap-B) suffices. Validity witnesses are private:
every valid code retains its original minimum energy over existing auxiliaries. This
does not certify the original constraint semantics or a minimal validity weight.
"""
function guard_encodings(q::QUBOComponent{T};gap::Integer=1,max_terms::Integer=1_000_000,
        namespace::Symbol=:encoding_validity) where T
    T in (BigInt,Rational{BigInt}) || throw(ArgumentError("exact coefficients required"))
    max_terms>0 || throw(ArgumentError("positive term budget required"))
    gap>0 || throw(ArgumentError("positive integer gap required"))
    covered = BitID[b for book in q.codebooks for b in _book_bits(book)]
    Set(covered)==Set(filter(b->b.role===:primary,q.bits)) || throw(ArgumentError("incomplete primary codebooks"))
    lower = q.offset + sum(v->min(zero(T),v),nonzeros(q.linear);init=zero(T)) +
        sum(v->min(zero(T),v),nonzeros(q.quadratic);init=zero(T))
    weight = max(zero(T),T(gap)-lower)
    (iszero(weight) || isempty(q.codebooks)) && return (;component=q,lower_bound=lower,weight)
    parts = [encoding_validity(book;max_terms) for book in q.codebooks]
    all(p->iszero(p.offset) && nnz(p.linear)==0 && nnz(p.quadratic)==0,parts) &&
        return (;component=q,lower_bound=lower,weight)
    sum(nnz(p.linear)+nnz(p.quadratic) for p in parts)<=max_terms || throw(ArgumentError("combined validity exceeds max_terms"))
    bits = reduce(vcat,(p.bits for p in parts);init=BitID[])
    positions = Dict(b=>i for (i,b) in enumerate(bits))
    weighted = QUBOComponent{T}(bits;
        linear=[(positions[p.bits[j]],weight*v) for p in parts for (j,v) in _linear_terms(p)],
        quadratic=[(positions[p.bits[j]],positions[p.bits[k]],weight*v) for p in parts for (j,k,v) in _quadratic_terms(p)],
        offset=weight*sum(p.offset for p in parts),codebooks=q.codebooks,
        auxiliary_meanings=Dict(b=>meaning for p in parts for (b,meaning) in p.auxiliary_meanings),
        provenance="encoding guard v1; B=$(lower); weight=$(weight); gap=$(gap)")
    return (;component=compose(q,weighted;namespace),lower_bound=lower,weight)
end

"""Compile sums of unary/binary violations using affine value indicators and global validity.

Supported structured books: one_hot, domain_wall, zero_one_hot, native. No semantic or
quadratization auxiliaries are required. Signed indicator products are protected by an
exact coefficient bound. Valid-code energy counts violated relations, not their indicator.
"""
function local_component(books;unary=UnaryRelation[],binary=BinaryRelation[],max_terms::Integer=1_000_000,
        provenance::AbstractString="local relation tables v1")
    books = collect(books)
    max_terms>0 || throw(ArgumentError("positive term budget required"))
    all(b->b isa StructuredCodebook && b.encoding in (:one_hot,:domain_wall,:zero_one_hot,:native),books) ||
        throw(ArgumentError("local affine compiler requires supported structured books"))
    allunique(b.variable for b in books) || throw(ArgumentError("duplicate codebook variable"))
    validity_size = sum((b.encoding===:one_hot ? big(length(b.bits))*(length(b.bits)+1)÷2 :
        b.encoding===:zero_one_hot ? big(length(b.bits))*(length(b.bits)-1)÷2 :
        b.encoding===:domain_wall ? 2big(max(0,length(b.bits)-1)) : big(length(b.values)==1)) for b in books;init=big(0))
    validity_size<=max_terms || throw(CompilationLimit(:affine_validity_terms,validity_size,big(max_terms)))
    byname = Dict(b.variable=>b for b in books)
    bits = reduce(vcat,(b.bits for b in books);init=BitID[])
    positions = Dict(b=>i for (i,b) in enumerate(bits))
    indicators = Dict{Symbol,Vector{Tuple{BigInt,Vector{Tuple{Int,BigInt}}}}}()
    for book in books
        ids = [positions[b] for b in book.bits]
        indicators[book.variable] = [_value_indicator(book,k,ids) for k in eachindex(book.values)]
    end
    linear = Tuple{Int,BigInt}[]
    quadratic = Tuple{Int,Int,BigInt}[]
    offset = big(0)
    function add_affine!(indicator,weight=big(1))
        offset += weight*indicator[1]
        append!(linear,((i,weight*c) for (i,c) in indicator[2]))
    end
    function budget!()
        length(linear)+length(quadratic)<=max_terms || throw(ArgumentError("local relation generation exceeds max_terms"))
    end
    for relation in unary
        haskey(byname,relation.variable) || throw(ArgumentError("unknown unary variable"))
        length(relation.allowed)==length(byname[relation.variable].values) || throw(DimensionMismatch("unary table size"))
        for a in eachindex(relation.allowed)
            relation.allowed[a] || add_affine!(indicators[relation.variable][a])
            budget!()
        end
    end
    for relation in binary
        haskey(byname,relation.left) && haskey(byname,relation.right) || throw(ArgumentError("unknown binary variable"))
        size(relation.allowed)==(length(byname[relation.left].values),length(byname[relation.right].values)) ||
            throw(DimensionMismatch("binary table size"))
        for b in axes(relation.allowed,2), a in axes(relation.allowed,1)
            relation.allowed[a,b] && continue
            left,right = indicators[relation.left][a],indicators[relation.right][b]
            add_affine!(left,right[1])
            append!(linear,((i,left[1]*c) for (i,c) in right[2]))
            append!(quadratic,((i,j,c*d) for (i,c) in left[2] for (j,d) in right[2]))
            budget!()
        end
    end
    raw = QUBOComponent(bits;linear,quadratic,offset,codebooks=books,provenance)
    return guard_encodings(raw;max_terms).component
end

function _local_books(domains,variables,encoding)
    length(variables)==length(domains) && allunique(variables) || throw(ArgumentError("invalid variable identities"))
    modes = encoding isa Symbol ? fill(encoding,length(domains)) : collect(encoding)
    length(modes)==length(domains) || throw(DimensionMismatch("one encoding per variable required"))
    return [structured_codebook(variables[i],domains[i];encoding=modes[i]) for i in eachindex(domains)]
end

function _local_table_budget(domains,family,limit)
    limit>0 || throw(ArgumentError("positive table budget required"))
    n = length(domains)
    sizes = [big(length(d)) for d in domains]
    entries = sum(sizes;init=big(0))
    if family===:table && n==2
        entries += sizes[1]*sizes[2]
    elseif family!==:instantiation && family!==:table
        for i in 1:n, j in i+1:n
            family in (:all_equal,:ordered) && j!=i+1 && continue
            entries += sizes[i]*sizes[j]
        end
    end
    entries<=limit || throw(CompilationLimit(:semantic_table_entries,entries,big(limit)))
end

"""Exact scalar local families with explicit finite domains and XCSP-style parameters.

Families: instantiation, all_equal, all_different (except), ordered (lt/le/ge/gt and fixed
distances), channel (single list, start_index=0), no_overlap (fixed lengths, zero_ignored=true).
This is a compiler API, not an XML parser or coverage claim for other variants.
"""
function local_constraint(family::Symbol,domains;variables=[Symbol(:x,i) for i in eachindex(domains)],
        encoding=:one_hot,values=nothing,except=[],operator::Symbol=:le,distances=nothing,
        lengths=nothing,zero_ignored::Bool=true,start_index::Integer=0,max_terms::Integer=1_000_000,
        max_table_entries::Integer=1_000_000)
    family in (:instantiation,:all_equal,:all_different,:ordered,:channel,:no_overlap) || throw(ArgumentError("unsupported local family"))
    values===nothing || family===:instantiation || throw(ArgumentError("values only applies to instantiation"))
    isempty(except) || family===:all_different || throw(ArgumentError("except only applies to all_different"))
    (operator===:le && distances===nothing) || family===:ordered || throw(ArgumentError("order parameters only apply to ordered"))
    (lengths===nothing && zero_ignored) || family===:no_overlap || throw(ArgumentError("duration parameters only apply to no_overlap"))
    iszero(start_index) || family===:channel || throw(ArgumentError("start_index only applies to channel"))
    _local_table_budget(domains,family,max_table_entries)
    books = _local_books(domains,variables,encoding)
    n = length(books)
    unary,binary = UnaryRelation[],BinaryRelation[]
    if family===:instantiation
        values!==nothing && length(values)==n || throw(ArgumentError("one target value per variable required"))
        append!(unary,(UnaryRelation(variables[i],[isequal(v,values[i]) for v in books[i].values]) for i in 1:n))
    else
        if family in (:ordered,:channel,:no_overlap)
            all(b->all(v->v isa Integer,b.values),books) || throw(ArgumentError("integer domains required"))
        end
        if family===:ordered
            operator in (:lt,:le,:ge,:gt) || throw(ArgumentError("unsupported order operator"))
            distances = distances===nothing ? zeros(Int,max(0,n-1)) : collect(distances)
            length(distances)==max(0,n-1) && all(d->d isa Integer,distances) || throw(ArgumentError("fixed integer distances required"))
        elseif family===:no_overlap
            lengths!==nothing && length(lengths)==n && all(l->l isa Integer && l>=0,lengths) ||
                throw(ArgumentError("one fixed nonnegative length per task required"))
        elseif family===:channel
            append!(unary,(UnaryRelation(variables[i],[start_index<=v<big(start_index)+n for v in books[i].values]) for i in 1:n))
        end
        for i in 1:n, j in i+1:n
            family in (:ordered,:all_equal) && j!=i+1 && continue
            allowed = [if family===:all_equal
                isequal(x,y)
            elseif family===:all_different
                !isequal(x,y) || any(v->isequal(x,v),except)
            elseif family===:ordered
                a,b = big(x)+big(distances[i]),big(y)
                operator===:lt ? a<b : operator===:le ? a<=b : operator===:ge ? a>=b : a>b
            elseif family===:channel
                (x==big(start_index)+j-1)==(y==big(start_index)+i-1)
            else
                (zero_ignored && (iszero(lengths[i]) || iszero(lengths[j]))) ||
                    big(x)+big(lengths[i])<=big(y) || big(y)+big(lengths[j])<=big(x)
            end for x in books[i].values, y in books[j].values]
            push!(binary,BinaryRelation(variables[i],variables[j],allowed))
        end
    end
    return local_component(books;unary,binary,max_terms,
        provenance="local family v1; $(family); values=$(values); except=$(except); operator=$(operator); distances=$(distances); lengths=$(lengths); zero_ignored=$(zero_ignored); start_index=$(start_index)")
end

"""Unary/binary supports or conflicts with explicit ANY_VALUE short-tuple wildcard."""
function table_constraint(domains,tuples;supports::Bool=true,
        variables=[Symbol(:x,i) for i in eachindex(domains)],encoding=:one_hot,max_terms::Integer=1_000_000,
        max_table_entries::Integer=1_000_000)
    n = length(domains)
    1<=n<=2 || throw(ArgumentError("local table compiler supports arity one or two"))
    _local_table_budget(domains,:table,max_table_entries)
    rows = collect(tuples)
    all(t->length(t)==n,rows) || throw(DimensionMismatch("table tuple arity"))
    books = _local_books(domains,variables,encoding)
    all(b->all(v->!(v isa AnyValue),b.values),books) || throw(ArgumentError("wildcard cannot be a semantic domain value"))
    matches(v) = any(t->all(t[i] isa AnyValue || isequal(t[i],v[i]) for i in 1:n),rows)
    allowed(v) = matches(v)==supports
    if n==1
        relation = UnaryRelation(only(variables),[allowed((v,)) for v in only(books).values])
        return local_component(books;unary=[relation],max_terms,provenance="unary table v1; supports=$(supports); rows=$(rows)")
    end
    relation = BinaryRelation(variables[1],variables[2],[allowed((a,b)) for a in books[1].values, b in books[2].values])
    return local_component(books;binary=[relation],max_terms,provenance="binary table v1; supports=$(supports); rows=$(rows)")
end
