# Q7 uses the explicitly selected whole-integer-time operational profile.
function _task_terms(books, terms; nonnegative=false)
    positions=Dict(b.variable=>i for (i,b) in enumerate(books))
    terms=collect(terms)
    for t in terms
        t isa Integer || t isa Symbol && haskey(positions,t) || throw(ArgumentError("task term must be an integer or known variable"))
        values=t isa Symbol ? _book_values(books[positions[t]]) : (t,)
        nonnegative && any(<(0),values) && throw(ArgumentError("negative duration or height domain"))
    end
    return terms,positions
end
_task_value(t,v,positions)=t isa Symbol ? big(v[positions[t]]) : big(t)

"""Cumulative over ALL integer times, including the exterior zero load.

Fixed nonnegative lengths/heights use event-indexed separable sums, without scanning
the numeric horizon. Variable lengths/heights or explicit ends use a bounded finite
truth construction. Terms are explicit integers or variable Symbols, not XML tokens.
"""
function cumulative_component(books,origins,lengths,heights,condition::IntegerCondition;
        ends=nothing,route::Symbol=:auto,time_domain::Symbol=:integers,
        namespace::Symbol=:cumulative_witness,max_events::Integer=10000,
        max_terms::Integer=100000,max_auxiliaries::Integer=1024,kwargs...)
    time_domain===:integers || throw(ArgumentError("Q7 operational profile requires time_domain=:integers"))
    route in (:auto,:events,:finite) || throw(ArgumentError("unknown cumulative route"))
    max_events>=0 && max_terms>0 && max_auxiliaries>=0 || throw(ArgumentError("invalid cumulative budget"))
    books,_,_=_finite_scope(books,Symbol[])
    origins,pos=_task_terms(books,origins)
    lengths,_=_task_terms(books,lengths;nonnegative=true)
    heights,_=_task_terms(books,heights;nonnegative=true)
    length(origins)==length(lengths)==length(heights) || throw(DimensionMismatch("task vectors"))
    if ends!==nothing
        ends,_=_task_terms(books,ends)
        length(ends)==length(origins) || throw(DimensionMismatch("task ends"))
    end
    fixed=all(t->t isa Integer,[lengths;heights]) && ends===nothing
    route===:events && !fixed && throw(ArgumentError("event route requires fixed lengths/heights and no explicit ends"))
    if route===:finite || !fixed
        function predicate(v)
            o=[_task_value(t,v,pos) for t in origins]
            l=[_task_value(t,v,pos) for t in lengths]
            h=[_task_value(t,v,pos) for t in heights]
            ends===nothing || all(i->_task_value(ends[i],v,pos)==o[i]+l[i],eachindex(o)) || return false
            satisfies(condition,0) || return false
            events=unique([o;o.+l])
            return all(t->satisfies(condition,sum((h[i] for i in eachindex(o) if o[i]<=t<o[i]+l[i]);init=big(0))),events)
        end
        return truth_component(books,predicate;oracle_id="Q7/cumulative/Z/finite-v1/$(origins)/$(lengths)/$(heights)/$(ends)/$(condition)",
            namespace,max_terms,max_auxiliaries,kwargs...)
    end
    isempty(kwargs) || throw(ArgumentError("finite-only options supplied to event route"))
    books,_,_=_structural_books(books,Symbol[],namespace)
    events=Set{BigInt}()
    for i in eachindex(origins)
        values=origins[i] isa Symbol ? books[pos[origins[i]]].values : (origins[i],)
        for o in values, t in (big(o),big(o)+lengths[i])
            push!(events,t)
            length(events)<=max_events || throw(CompilationLimit(:cumulative_events,big(length(events)),big(max_events)))
        end
    end
    q=_structural_constant(books,Int(!satisfies(condition,0));max_terms)
    for (j,t) in enumerate(sort!(collect(events)))
        maps=[zeros(BigInt,length(b.values)) for b in books]
        constant=big(0)
        for i in eachindex(origins)
            if origins[i] isa Symbol
                k=pos[origins[i]]
                for (a,o) in enumerate(books[k].values)
                    o<=t<big(o)+lengths[i] && (maps[k][a]+=heights[i])
                end
            elseif origins[i]<=t<big(origins[i])+lengths[i]
                constant+=heights[i]
            end
        end
        # Shift the condition operand, preserving finite sets and intervals.
        rhs=condition.operand
        shifted=condition.operator in (:in,:notin) ?
            (rhs isa AbstractRange ? ((first(rhs)-constant):step(rhs):(last(rhs)-constant)) : rhs .- constant) : rhs-constant
        part=sum_component(books,IntegerCondition(condition.operator,shifted);value_maps=maps,
            namespace=Symbol(namespace,"/event/",j),max_terms)
        q=compose(q,part;namespace=Symbol(namespace,"/private/",j))
        _aggregate_finish(q;max_terms,max_auxiliaries)
    end
    return _aggregate_finish(q;max_terms,max_auxiliaries)
end

"""Pairwise finite noOverlap factors for rectangular n-by-k matrices of task terms.

Vector inputs denote 1D. With zero_ignored, an object with ANY zero extent is ignored.
Each pair is quadratized separately, with shared primaries and private auxiliaries.
This is exponential in the bits of a pair, not a scalable orientation formulation.
"""
function nooverlap_component(books,origins,lengths;zero_ignored::Bool=true,
        namespace::Symbol=:nooverlap_witness,max_pairs::Integer=10000,
        max_terms::Integer=100000,max_auxiliaries::Integer=1024,kwargs...)
    max_pairs>=0 && max_terms>0 && max_auxiliaries>=0 || throw(ArgumentError("invalid noOverlap budgets"))
    books,_,_=_structural_books(books,Symbol[],namespace)
    _finite_scope(books,Symbol[])
    o=origins isa AbstractMatrix ? Matrix(origins) : reshape(collect(origins),:,1)
    l=lengths isa AbstractMatrix ? Matrix(lengths) : reshape(collect(lengths),:,1)
    size(o)==size(l) && size(o,2)>0 || throw(DimensionMismatch("matching nonempty dimensions"))
    _task_terms(books,vec(o)); _task_terms(books,vec(l);nonnegative=true)
    n,k=size(o); pairs=big(n)*(n-1)÷2
    pairs<=max_pairs || throw(CompilationLimit(:nooverlap_pairs,pairs,big(max_pairs)))
    q=_structural_constant(books,0;max_terms)
    for i in 1:n, j in i+1:n
        names=Set(t for t in [vec(o[[i,j],:]);vec(l[[i,j],:])] if t isa Symbol)
        subset=[b for b in books if b.variable in names]
        pos=Dict(b.variable=>a for (a,b) in enumerate(subset))
        predicate=v->begin
            a=[_task_value(t,v,pos) for t in o[i,:]]; b=[_task_value(t,v,pos) for t in o[j,:]]
            da=[_task_value(t,v,pos) for t in l[i,:]]; db=[_task_value(t,v,pos) for t in l[j,:]]
            (zero_ignored && (any(iszero,da) || any(iszero,db))) ||
                any(d->a[d]+da[d]<=b[d] || b[d]+db[d]<=a[d],1:k)
        end
        part=truth_component(subset,predicate;oracle_id="Q7/noOverlap/pair-v1/$(o[[i,j],:])/$(l[[i,j],:])/$(zero_ignored)",
            namespace=Symbol(namespace,"/pair/",i,"/",j),max_terms,max_auxiliaries,kwargs...)
        q=compose(q,part;namespace=Symbol(namespace,"/private/",i,"/",j))
        _aggregate_finish(q;max_terms,max_auxiliaries)
    end
    return q
end

"""Bounded finite single-subcircuit predicate, NOT just allDifferent successors.

There is exactly one nontrivial cycle; other vertices are fixed points. Optional
size is an integer or a variable identity. start_index is explicit (core default 0).
"""
function circuit_component(books,scope;size=nothing,start_index::Integer=0,kwargs...)
    books,indices,pos=_finite_scope(books,scope)
    n=length(indices)
    size===nothing || _task_terms(books,[size])
    function predicate(v)
        successors=[big(v[i])-start_index for i in indices]
        all(x->0<=x<n,successors) && allunique(successors) || return false
        active=findall(i->successors[i]!=i-1,1:n)
        length(active)>1 || return false
        size===nothing || _task_value(size,v,pos)==length(active) || return false
        seen=Set{Int}(); node=first(active)
        while !(node in seen)
            push!(seen,node); node=Int(successors[node])+1
        end
        return node==first(active) && length(seen)==length(active)
    end
    return truth_component(books,predicate;oracle_id="Q7/circuit/subcycle-finite-v1/$(indices)/$(size)/$(start_index)",kwargs...)
end
