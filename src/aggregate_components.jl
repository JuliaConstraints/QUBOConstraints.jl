# Q5 factors are nonnegative on every bit string before existential projection.
# Shared semantic witnesses are primary during factor assembly, then hidden ONCE.
function _hide_witnesses(q::QUBOComponent{T},names) where T
    names = Set(names)
    rename(bit) = bit.role===:primary && bit.owner in names ?
        BitID(bit.owner,bit.index;role=:semantic_auxiliary) : bit
    bits = rename.(q.bits)
    allunique(bits) || throw(ArgumentError("hidden witness identity collision"))
    meanings = Dict(rename(b)=>v for (b,v) in q.auxiliary_meanings)
    for b in q.bits
        rename(b)!=b && (meanings[rename(b)]="shared semantic witness $(b.owner), bit $(b.index)")
    end
    return QUBOComponent{T}(bits;linear=_linear_terms(q),quadratic=_quadratic_terms(q),offset=q.offset,
        codebooks=[b for b in q.codebooks if !(b.variable in names)],auxiliary_meanings=meanings,
        applicability=q.applicability*"; existential semantic witnesses",provenance=q.provenance)
end

function _aggregate_books(books,namespace)
    books = collect(books)
    all(b->b isa StructuredCodebook,books) || throw(ArgumentError("structured codebooks required"))
    allunique(b.variable for b in books) || throw(ArgumentError("duplicate variable identity"))
    prefix = string(namespace)*"/"
    any(b->b.variable===namespace || startswith(string(b.variable),prefix),books) &&
        throw(ArgumentError("reserved witness namespace overlaps a primary variable"))
    return books
end

function _aggregate_validity(books;max_terms)
    # A zero contribution introduces no result witness and keeps every source codebook.
    return sum_component(books,IntegerCondition(:eq,0);
        value_maps=[zeros(Int,length(b.values)) for b in books],max_terms)
end

function _aggregate_finish(q;max_terms,max_auxiliaries)
    count = big(nnz(q.linear))+nnz(q.quadratic)
    count<=max_terms || throw(CompilationLimit(:aggregate_terms,count,big(max_terms)))
    count = big(Base.count(b->b.role!==:primary,q.bits))
    count<=max_auxiliaries || throw(CompilationLimit(:aggregate_auxiliaries,count,big(max_auxiliaries)))
    return q
end

function _arithmetic_product!(b,left,right,weight)
    a,rows = left; c,cols = right
    _arithmetic_budget(b,big(length(rows))*length(cols)+length(rows)+length(cols))
    b.offset += weight*a*c
    append!(b.linear,((i,weight*c*w) for (i,w) in rows))
    append!(b.linear,((j,weight*a*v) for (j,v) in cols))
    append!(b.quadratic,((i,j,weight*w*v) for (i,w) in rows for (j,v) in cols))
end
_bit_affine(i) = (big(0),[(i,big(1))])
_constant_affine(v) = (big(v),Tuple{Int,BigInt}[])

function _boolean_gate!(b,u,v,y,op)
    one = _constant_affine(1)
    if op===:or
        for a in (u,v,y); _arithmetic_product!(b,a,one,1); end
        _arithmetic_product!(b,u,v,1)
        _arithmetic_product!(b,u,y,-2); _arithmetic_product!(b,v,y,-2)
    elseif op===:and
        _arithmetic_product!(b,u,v,1)
        _arithmetic_product!(b,u,y,-2); _arithmetic_product!(b,v,y,-2)
        _arithmetic_product!(b,y,one,3)
    else
        throw(ArgumentError("unsupported Boolean gate"))
    end
end

# Return an exact relation between source indicators and native Boolean outputs.
# Signed affine indicators outside their codebooks require a global validity guard.
function _aggregate_signals(books,outputs,maps,op;namespace,max_terms,max_auxiliaries)
    max_auxiliaries>=0 || throw(ArgumentError("nonnegative auxiliary budget required"))
    needed = big(length(outputs))*max(0,length(books)-2)+sum(length(b.bits) for b in outputs;init=0)
    needed<=max_auxiliaries || throw(CompilationLimit(:aggregate_witnesses,needed,big(max_auxiliaries)))
    allbooks = [books;outputs]
    bits = reduce(vcat,(b.bits for b in allbooks);init=BitID[])
    b = _ArithmeticBuilder(bits,Int(max_terms))
    ids = Dict(bit=>i for (i,bit) in enumerate(bits))
    gate = 0
    for (output,group) in zip(outputs,maps)
        inputs = [_arithmetic_affine(book,map,ids) for (book,map) in zip(books,group)]
        target = _bit_affine(ids[only(output.bits)])
        if length(inputs)<=1
            source = isempty(inputs) ? _constant_affine(op===:and) : only(inputs)
            _arithmetic_square!(b,source[1]-target[1],[source[2];[(i,-w) for (i,w) in target[2]]])
        else
            accumulator = first(inputs)
            for k in 2:length(inputs)
                next = if k==length(inputs)
                    target
                else
                    gate+=1
                    bit = BitID(Symbol(namespace,"/gate"),gate;role=:quadratization_auxiliary)
                    push!(b.bits,bit); b.meanings[bit]="$(op) prefix gate $(gate)"
                    _bit_affine(length(b.bits))
                end
                _boolean_gate!(b,accumulator,inputs[k],next,op)
                accumulator = next
            end
        end
    end
    q = QUBOComponent(b.bits;linear=b.linear,quadratic=b.quadratic,offset=b.offset,codebooks=allbooks,
        auxiliary_meanings=b.meanings,provenance="Q5 Boolean aggregate $(op); guarded affine indicators")
    return guard_encodings(q;max_terms,namespace=Symbol(namespace,"/validity")).component
end

"""Number of distinct semantic values, excluding a constant `except` set.

Uses OR chains for value occurrence and a separable count condition. General value
indicators require one_hot/domain_wall/zero_one_hot (or affine two-value books).
The witness bound is constructive, not minimal; all factors have a unit zero-set gap.
"""
function nvalues_component(books,condition::IntegerCondition;except=Int[],namespace::Symbol=:nvalues,
        max_terms::Integer=1000000,max_auxiliaries::Integer=10000,kwargs...)
    books = _aggregate_books(books,namespace)
    universe = unique([v for b in books for v in b.values if !(v in except)])
    outputs = [structured_codebook(Symbol(namespace,"/occurs/",i),[0,1];encoding=:native) for i in eachindex(universe)]
    maps = [[[Int(isequal(v,target)) for v in b.values] for b in books] for target in universe]
    relation = _aggregate_signals(books,outputs,maps,:or;namespace,max_terms,max_auxiliaries)
    condition_part = sum_component(outputs,condition;max_terms,namespace=Symbol(namespace,"/result"),kwargs...)
    q = _hide_witnesses(compose(relation,condition_part;namespace=Symbol(namespace,"/condition")),[b.variable for b in outputs])
    return _aggregate_finish(q;max_terms,max_auxiliaries)
end

"""Cardinality with constant, distinct counted values and one integer condition per value.

`closed=false` leaves unlisted values unrestricted. Variable values/occurrences use the
explicit finite compiler, not this constant-parameter interface.
"""
function cardinality_component(books,values,conditions;closed::Bool=false,
        namespace::Symbol=:cardinality,max_terms::Integer=1000000,max_auxiliaries::Integer=10000,kwargs...)
    books = _aggregate_books(books,namespace)
    values,conditions = collect(values),collect(conditions)
    allunique(values) && length(values)==length(conditions) && all(c->c isa IntegerCondition,conditions) ||
        throw(ArgumentError("distinct values and one IntegerCondition per value required"))
    q = _aggregate_validity(books;max_terms)
    for (i,(v,c)) in enumerate(zip(values,conditions))
        part = count_component(books,[v],c;max_terms,namespace=Symbol(namespace,"/count/",i),kwargs...)
        q = compose(q,part;namespace=Symbol(namespace,"/factor/",i))
    end
    if closed
        part = count_component(books,values,IntegerCondition(:eq,length(books));max_terms)
        q = compose(q,part;namespace=Symbol(namespace,"/closed"))
    end
    return _aggregate_finish(q;max_terms,max_auxiliaries)
end

"""Minimum/maximum through threshold AND/OR chains and a domain-wall result witness.

Nonempty integer scope; no arg-index/rank semantics are implied by this interface.
"""
function extremum_component(books,kind::Symbol,condition::IntegerCondition;
        namespace::Symbol=:extremum,max_terms::Integer=1000000,max_auxiliaries::Integer=10000,kwargs...)
    books = _aggregate_books(books,namespace)
    kind in (:minimum,:maximum) || throw(ArgumentError("minimum or maximum required"))
    !isempty(books) && all(b->all(v->v isa Integer,b.values),books) || throw(ArgumentError("nonempty integer scope required"))
    universe = sort!(unique(BigInt[v for b in books for v in b.values]))
    outputs = [structured_codebook(Symbol(namespace,"/threshold/",i),[0,1];encoding=:native) for i in 2:length(universe)]
    maps = [[[Int(v>=t) for v in b.values] for b in books] for t in universe[2:end]]
    relation = _aggregate_signals(books,outputs,maps,kind===:minimum ? :and : :or;namespace,max_terms,max_auxiliaries)
    # min/max = first(universe) + sum successive differences * threshold indicators.
    shifted = condition.operator in (:in,:notin) ?
        IntegerCondition(condition.operator,condition.operand .- first(universe)) :
        IntegerCondition(condition.operator,condition.operand-first(universe))
    part = sum_component(outputs,shifted;coefficients=diff(universe),max_terms,
        namespace=Symbol(namespace,"/result"),kwargs...)
    q = _hide_witnesses(compose(relation,part;namespace=Symbol(namespace,"/condition")),[b.variable for b in outputs])
    return _aggregate_finish(q;max_terms,max_auxiliaries)
end

"""Constant-weight/profit knapsack: conjunction of the two specified integer conditions."""
function knapsack_component(books,weights,profits,weight_condition::IntegerCondition,profit_condition::IntegerCondition;
        namespace::Symbol=:knapsack,max_terms::Integer=1000000,max_auxiliaries::Integer=10000,kwargs...)
    books = _aggregate_books(books,namespace)
    a = sum_component(books,weight_condition;coefficients=weights,namespace=Symbol(namespace,"/weight"),max_terms,kwargs...)
    b = sum_component(books,profit_condition;coefficients=profits,namespace=Symbol(namespace,"/profit"),max_terms,kwargs...)
    return _aggregate_finish(compose(a,b;namespace=Symbol(namespace,"/profit_factor"));max_terms,max_auxiliaries)
end

"""Bin packing with nonnegative constant sizes and an explicit empty-bin policy.

Single `IntegerCondition`: only occupied bins (XCSP3 common-condition form).
Vector of conditions plus `bins`: every declared bin, including empty, by default
(`profile=:validator`). `profile=:pdf_occupied` retains the alternative literal
occupied-only reading for declared conditions. Out-of-declared-bin domains are rejected.
Zero-size objects still occupy bins. Load witnesses use domain-wall encoding.
"""
function binpacking_component(books,sizes,conditions;bins=nothing,profile::Symbol=:validator,
        namespace::Symbol=:binpacking,max_terms::Integer=1000000,max_auxiliaries::Integer=10000,
        max_witness_values::Integer=100000,kwargs...)
    books = _aggregate_books(books,namespace)
    sizes = collect(sizes)
    length(sizes)==length(books) && all(s->s isa Integer && s>=0,sizes) || throw(ArgumentError("one nonnegative integer size per object required"))
    all(b->all(v->v isa Integer,b.values),books) || throw(ArgumentError("integer bin domains required"))
    profile in (:validator,:pdf_occupied) || throw(ArgumentError("unknown bin-packing profile"))
    common = conditions isa IntegerCondition
    common && bins!==nothing && throw(ArgumentError("single condition derives occupied bins from domains"))
    !common && bins===nothing && throw(ArgumentError("declared conditions require explicit bins"))
    labels = common ? sort!(unique(BigInt[v for b in books for v in b.values])) : collect(bins)
    all(v->v isa Integer,labels) && allunique(labels) || throw(ArgumentError("distinct integer bins required"))
    cs = common ? fill(conditions,length(labels)) : collect(conditions)
    length(cs)==length(labels) && all(c->c isa IntegerCondition,cs) || throw(ArgumentError("one condition per bin required"))
    all(b->all(in(labels),b.values),books) || throw(ArgumentError("source domain outside declared bins"))
    q = _aggregate_validity(books;max_terms)
    for (i,(label,c)) in enumerate(zip(labels,cs))
        maps = [[Int(v==label) for v in b.values] for b in books]
        if !common && profile===:validator
            part = sum_component(books,c;coefficients=sizes,value_maps=maps,max_terms,max_witness_values,
                namespace=Symbol(namespace,"/load/",i),kwargs...)
        else
            limit = sum((big(s) for (b,s) in zip(books,sizes) if label in b.values);init=big(0))
            limit+1<=max_witness_values || throw(CompilationLimit(:bin_load_values,limit+1,big(max_witness_values)))
            limit+1+max(0,length(books)-2)<=max_auxiliaries ||
                throw(CompilationLimit(:bin_witnesses,limit+1+max(0,length(books)-2),big(max_auxiliaries)))
            load = structured_codebook(Symbol(namespace,"/load/",i),collect(big(0):limit);encoding=:domain_wall)
            used = structured_codebook(Symbol(namespace,"/used/",i),[0,1];encoding=:native)
            link = sum_component([books;[load]],IntegerCondition(:eq,0);
                coefficients=[sizes;[-1]],value_maps=[maps;[load.values]],max_terms,max_witness_values)
            signal = _aggregate_signals(books,[used],[maps],:or;
                namespace=Symbol(namespace,"/occupancy/",i),max_terms,max_auxiliaries)
            allowed = [u==0 || satisfies(c,v) for u in used.values,v in load.values]
            rule = local_component([used,load];binary=[BinaryRelation(used.variable,load.variable,allowed)],max_terms)
            part = compose(compose(link,signal;namespace=Symbol(namespace,"/signal/",i)),rule;
                namespace=Symbol(namespace,"/rule/",i))
            part = _hide_witnesses(part,[load.variable,used.variable])
        end
        q = compose(q,part;namespace=Symbol(namespace,"/bin/",i))
        _aggregate_finish(q;max_terms,max_auxiliaries)
    end
    return q
end

"""Explicit bounded fallback for arithmetic variants with variable parameters.

The pure callback must encode the complete declared semantics (including repeated scopes,
variable coefficients/conditions, and index ranks). This is a finite truth construction,
NOT a scalable XCSP3 adapter or a claim that unspecified variants were implemented.
"""
function finite_arithmetic_component(books,predicate;semantic_id::AbstractString,kwargs...)
    return truth_component(books,predicate;oracle_id=semantic_id,kwargs...)
end
