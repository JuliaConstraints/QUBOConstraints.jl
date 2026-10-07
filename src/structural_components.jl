# Q6 structural constructions. All factors are nonnegative before projection.
function _structural_books(books, scope, namespace)
    books = _aggregate_books(books, namespace)
    all(b->b.encoding in (:one_hot,:domain_wall,:zero_one_hot,:native),books) ||
        throw(ArgumentError("structural affine construction requires affine indicator encodings"))
    scope = Symbol.(collect(scope))
    byname = Dict(b.variable=>b for b in books)
    all(s->haskey(byname,s),scope) || throw(ArgumentError("unknown scope variable"))
    return books,scope,byname
end

function _structural_constant(books, value; max_terms)
    base = local_component(books;max_terms)
    return QUBOComponent(base.bits;linear=_linear_terms(base),quadratic=_quadratic_terms(base),
        offset=base.offset+value,codebooks=base.codebooks,
        auxiliary_meanings=base.auxiliary_meanings,provenance="Q6 constant $(value) with primary validity")
end

"""N-ary supports via a private tuple selector, or conflicts via separable counts.

Repeated scope identities and explicit ANY_VALUE wildcards are supported. Empty supports
are false and empty conflicts true. This preserves the exact zero set with unit gap,
not the indicator value, and introduces no implicit XML semantics. Primary encodings
must have affine indicators. Budgets describe constructive limits, never impossibility.
"""
function table_component(books, scope, tuples;supports::Bool=true,
        selector_encoding::Symbol=:one_hot,namespace::Symbol=:table_witness,
        max_rows::Integer=10000,max_terms::Integer=100000,max_auxiliaries::Integer=1024)
    max_rows>=0 && max_terms>0 && max_auxiliaries>=0 || throw(ArgumentError("invalid structural budgets"))
    selector_encoding in (:one_hot,:domain_wall,:zero_one_hot) || throw(ArgumentError("unsupported selector encoding"))
    books,scope,byname = _structural_books(books,scope,namespace)
    any(b->any(v->v isa AnyValue,b.values),books) && throw(ArgumentError("wildcard cannot be a domain value"))
    rows = Vector{Any}[]
    visited = big(0)
    for tuple in tuples
        visited+=1
        visited<=max_rows || throw(CompilationLimit(:table_rows,visited,big(max_rows)))
        length(tuple)==length(scope) || throw(DimensionMismatch("table arity"))
        row = Any[tuple...]
        # Delete unreachable rows, including incompatible repeated occurrences.
        bindings = Dict{Symbol,Any}()
        possible = true
        for (name,value) in zip(scope,row)
            value isa AnyValue && continue
            if !any(v->isequal(value,v),byname[name].values) ||
                    (haskey(bindings,name) && !isequal(bindings[name],value))
                possible=false; break
            end
            bindings[name]=value
        end
        possible && !any(r->isequal(r,row),rows) && push!(rows,row)
    end
    if isempty(rows) || isempty(scope)
        return _structural_constant(books,Int(isempty(rows)==supports);max_terms)
    end
    if supports
        selector = structured_codebook(Symbol(namespace,"/tuple"),1:length(rows);encoding=selector_encoding)
        length(selector.bits)<=max_auxiliaries ||
            throw(CompilationLimit(:table_selectors,big(length(selector.bits)),big(max_auxiliaries)))
        entries = sum((big(length(rows))*length(byname[name].values) for name in scope);init=big(0))
        entries<=max_terms || throw(CompilationLimit(:table_relation_entries,entries,big(max_terms)))
        relations = [BinaryRelation(selector.variable,name,
            [rows[r][j] isa AnyValue || isequal(rows[r][j],v) for r in eachindex(rows),v in byname[name].values])
            for (j,name) in enumerate(scope)]
        q = local_component([books;[selector]];binary=relations,max_terms,
            provenance="Q6 supports selector v1; scope=$(scope); rows=$(rows); selector=$(selector_encoding)")
        return _aggregate_finish(_hide_witnesses(q,[selector.variable]);max_terms,max_auxiliaries)
    end
    q = _structural_constant(books,0;max_terms)
    for (r,row) in enumerate(rows)
        bindings = Dict(name=>value for (name,value) in zip(scope,row) if !(value isa AnyValue))
        maps = [[Int(haskey(bindings,b.variable) && isequal(v,bindings[b.variable])) for v in b.values] for b in books]
        part = sum_component(books,IntegerCondition(:le,length(bindings)-1);value_maps=maps,
            namespace=Symbol(namespace,"/conflict/",r),max_terms)
        q = compose(q,part;namespace=Symbol(namespace,"/private/",r))
        _aggregate_finish(q;max_terms,max_auxiliaries)
    end
    return q
end

"""Exact path acceptance in a layered labelled graph (MDD), using edge selectors.

Each layer contains (source,label,target) triples. States may be arbitrary immutable
identities, labels are exact semantic values, no epsilon transitions or implicit '*'.
Repeated input variables remain shared; witnesses for distinct layers remain distinct.
One layer per scope position, one start identity and an explicit set of final identities.
"""
function mdd_component(books,scope,layers;start,finals,edge_encoding::Symbol=:one_hot,
        namespace::Symbol=:path_witness,max_edges::Integer=10000,
        max_terms::Integer=100000,max_auxiliaries::Integer=1024)
    max_edges>=0 && max_terms>0 && max_auxiliaries>=0 || throw(ArgumentError("invalid path budgets"))
    edge_encoding in (:one_hot,:domain_wall,:zero_one_hot) || throw(ArgumentError("unsupported edge encoding"))
    books,scope,byname = _structural_books(books,scope,namespace)
    layers = collect(layers); finals=collect(finals)
    length(layers)==length(scope) || throw(DimensionMismatch("one graph layer per scope position"))
    isempty(scope) && return _structural_constant(books,Int(!(start in finals));max_terms)
    edges = Vector{Tuple{Any,Any,Any}}[]
    total=big(0)
    for (i,layer) in enumerate(layers)
        current=Tuple{Any,Any,Any}[]
        for e in layer
            total+=1
            total<=max_edges || throw(CompilationLimit(:path_edges,total,big(max_edges)))
            length(e)==3 || throw(DimensionMismatch("edge requires source, label, target"))
            e[2] isa AnyValue && throw(ArgumentError("path labels must be explicit values"))
            e[2] in byname[scope[i]].values && !(Tuple(e) in current) && push!(current,Tuple(e))
        end
        push!(edges,current)
    end
    # Forward and backward pruning preserves every full accepting path.
    reachable=Set([start])
    for layer in edges
        filter!(e->e[1] in reachable,layer)
        reachable=Set(e[3] for e in layer)
    end
    reachable=Set(finals)
    for layer in Iterators.reverse(edges)
        filter!(e->e[3] in reachable,layer)
        reachable=Set(e[1] for e in layer)
    end
    any(isempty,edges) && return _structural_constant(books,1;max_terms)
    witnesses=[structured_codebook(Symbol(namespace,"/edge/",i),1:length(e);encoding=edge_encoding) for (i,e) in enumerate(edges)]
    count=sum((big(length(b.bits)) for b in witnesses);init=big(0))
    count<=max_auxiliaries || throw(CompilationLimit(:path_witnesses,count,big(max_auxiliaries)))
    entries=sum((big(length(edges[i]))*length(byname[s].values) for (i,s) in enumerate(scope));init=big(0))+
        sum((big(length(edges[i]))*length(edges[i+1]) for i in 1:length(scope)-1);init=big(0))
    entries<=max_terms || throw(CompilationLimit(:path_relation_entries,entries,big(max_terms)))
    relations=BinaryRelation[]
    for (i,name) in enumerate(scope)
        push!(relations,BinaryRelation(witnesses[i].variable,name,[isequal(e[2],v) for e in edges[i],v in byname[name].values]))
        if i<length(scope)
            push!(relations,BinaryRelation(witnesses[i].variable,witnesses[i+1].variable,
                [isequal(e[3],f[1]) for e in edges[i],f in edges[i+1]]))
        end
    end
    q=local_component([books;witnesses];binary=relations,max_terms,
        provenance="Q6 layered path v1; scope=$(scope); edges=$(edges); start=$(start); finals=$(finals); encoding=$(edge_encoding)")
    return _aggregate_finish(_hide_witnesses(q,[b.variable for b in witnesses]);max_terms,max_auxiliaries)
end

"""Unroll a finite automaton into the layered exact path construction. Nondeterminism is allowed."""
function regular_component(books,scope,transitions;start,finals,max_edges::Integer=10000,kwargs...)
    rows=collect(transitions)
    count=big(length(rows))*length(scope)
    count<=max_edges || throw(CompilationLimit(:automaton_unrolled_edges,count,big(max_edges)))
    return mdd_component(books,scope,fill(rows,length(scope));start,finals,max_edges,kwargs...)
end

"""Window scopes with zero-based starts but ordinary Julia vector access.

Circular starts are k*offset<n, not an orbit until modular repetition. All positions
wrap. An empty sequence generates no windows. Positive width/offset are required.
"""
function slide_scopes(scope,width::Integer;offset::Integer=1,circular::Bool=false,max_windows::Integer=10000,
        max_scope_entries::Integer=100000)
    width>0 && offset>0 && max_windows>=0 && max_scope_entries>=0 || throw(ArgumentError("invalid slide dimensions or budget"))
    scope=collect(scope); n=length(scope)
    last_start=circular ? big(n)-1 : big(n)-width
    count=last_start<0 ? big(0) : fld(last_start,offset)+1
    count<=max_windows || throw(CompilationLimit(:slide_windows,count,big(max_windows)))
    count*width<=max_scope_entries || throw(CompilationLimit(:slide_scope_entries,count*width,big(max_scope_entries)))
    width<=typemax(Int) || throw(CompilationLimit(:slide_width,big(width),big(typemax(Int))))
    return [[scope[Int(mod(i+j,n))+1] for j in 0:Int(width)-1] for i in big(0):big(offset):last_start]
end

"""Compose trusted nonnegative, exact window factors with private auxiliary identities.

The builder receives (all_codebooks, window_scope, private_namespace). Its contract
must be established independently: an arbitrary callback is not a certificate. Each
factor must preserve precisely the supplied primary books and use exact coefficients.
"""
function slide_component(books,scope,width::Integer,builder;offset::Integer=1,circular::Bool=false,
        namespace::Symbol=:slide_witness,max_windows::Integer=10000,
        max_scope_entries::Integer=100000,max_terms::Integer=100000,max_auxiliaries::Integer=1024)
    books,scope,_=_structural_books(books,scope,namespace)
    windows=slide_scopes(scope,width;offset,circular,max_windows,max_scope_entries)
    q=_structural_constant(books,0;max_terms)
    for (i,window) in enumerate(windows)
        part=builder(books,window,Symbol(namespace,"/window/",i))
        part isa QUBOComponent{BigInt} || throw(ArgumentError("exact BigInt window component required"))
        length(part.codebooks)==length(books) && all(b->any(c->_compatible(b,c),part.codebooks),books) ||
            throw(ArgumentError("window builder changed primary codebooks"))
        q=compose(q,part;namespace=Symbol(namespace,"/private/",i))
        _aggregate_finish(q;max_terms,max_auxiliaries)
    end
    return q
end
