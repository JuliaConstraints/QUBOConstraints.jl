# Structural recipes. Lists contain integer constants, variable Symbols, or expression nodes.
# They expand at construction time to the closed atomic alphabet; no predicate callbacks.
_ax(op,args...)=IntensionNode(op,args...)
_afold(op,x,identity)=isempty(x) ? identity : length(x)==1 ? only(x) : _ax(op,x...)
_aall(x)=_afold(:and,collect(x),1)
_aany(x)=_afold(:or,collect(x),0)
_asum(x)=_afold(:add,collect(x),0)
_amember(x,values)=_aany(_ax(:eq,x,v) for v in values)
_acount(xs,values)=_asum(_amember(x,values) for x in xs)

"""A condition whose scalar operand may itself be a variable/expression.

For in/notin, an integer range is an interval; a vector is a finite set, whose
members may be expressions in the explicit normalized API (a superset of core).
"""
struct AtomicCondition
    operator::Symbol
    operand::Any
    function AtomicCondition(op::Symbol,operand)
        op in (:eq,:ne,:lt,:le,:gt,:ge,:in,:notin) || throw(ArgumentError("invalid atomic condition"))
        if op in (:in,:notin)
            operand isa Union{AbstractVector,AbstractRange,AbstractSet} || throw(ArgumentError("set/interval condition needs a collection"))
        else
            operand isa Union{Integer,Symbol,IntensionNode} || throw(ArgumentError("scalar condition needs an integer expression"))
        end
        new(op,deepcopy(operand))
    end
end
AtomicCondition(c::IntegerCondition)=AtomicCondition(c.operator,c.operand)
function _acond(x,c)
    c isa IntegerCondition && (c=AtomicCondition(c))
    c isa AtomicCondition || throw(ArgumentError("explicit AtomicCondition or IntegerCondition required"))
    op=c.operator; rhs=c.operand
    if op in (:in,:notin)
        yes=rhs isa AbstractRange && step(rhs)==1 ?
            isempty(rhs) ? 0 : _aall([_ax(:ge,x,first(rhs)),_ax(:le,x,last(rhs))]) : _amember(x,rhs)
        return op===:in ? yes : _ax(:not,yes)
    end
    _ax(op,x,rhs)
end

function _apairs(xs,op,except)
    _aall(_aany([_ax(op,xs[i],xs[j]),_amember(xs[i],except),_amember(xs[j],except)])
        for i in eachindex(xs) for j in i+1:length(xs))
end
function _alex(x,y,op)
    length(x)==length(y) || throw(DimensionMismatch("lex lists must have equal lengths"))
    op in (:lt,:le,:gt,:ge) || throw(ArgumentError("invalid lex order"))
    prefix=1; strict=0; cmp=op in (:lt,:le) ? :lt : :gt
    for (a,b) in zip(x,y)
        strict=_aany([strict,_aall([prefix,_ax(cmp,a,b)])])
        prefix=_aall([prefix,_ax(:eq,a,b)])
    end
    op in (:le,:ge) ? _aany([strict,prefix]) : strict
end
function _arank(matches,rank)
    rank in (:any,:first,:last) || throw(ArgumentError("rank must be any/first/last"))
    rank===:any && return matches
    [_aall([matches[i],_ax(:not,_aany(matches[j] for j in eachindex(matches) if rank===:first ? j<i : j>i))]) for i in eachindex(matches)]
end
function _arect(matrix)
    rows=matrix isa AbstractMatrix ? [collect(matrix[i,:]) for i in axes(matrix,1)] : [collect(row) for row in matrix]
    isempty(rows) && throw(ArgumentError("empty matrix"))
    all(row->length(row)==length(first(rows)),rows) && !isempty(first(rows)) || throw(DimensionMismatch("ragged/empty rows"))
    rows
end
_acolumns(rows)=[[row[j] for row in rows] for j in eachindex(first(rows))]

# Recover the layered MDD interface from the core transitions-only representation.
function _amdd_layers(transitions,n)
    edges=collect(transitions)
    !isempty(edges) && all(e->length(e)==3 && e[2] isa Integer,edges) || throw(ArgumentError("MDD needs labelled transitions"))
    sources=Set(e[1] for e in edges); targets=Set(e[3] for e in edges)
    roots=setdiff(sources,targets); sinks=setdiff(targets,sources)
    length(roots)==length(sinks)==1 || throw(ArgumentError("MDD requires a unique root and terminal"))
    root=only(roots); terminal=only(sinks); levels=Dict{Any,Int}(root=>0)
    layers=[Tuple[] for _ in 1:n]
    for depth in 0:n-1
        for (a,label,b) in edges
            get(levels,a,-1)==depth || continue
            haskey(levels,b) && levels[b]!=depth+1 && throw(ArgumentError("MDD is cyclic or not layered"))
            levels[b]=depth+1
            push!(layers[depth+1],(a,label,b))
        end
    end
    length(levels)==length(union(sources,targets)) && get(levels,terminal,-1)==n &&
        sum(length,layers)==length(edges) || throw(ArgumentError("MDD depth/disconnected transitions mismatch"))
    layers,root,[terminal]
end

"""Normalized XCSP3-core recipes; DOES NOT parse XML or silently infer variants.

Core3.2-v4 semantics and explicitly named operational profiles are recorded in
the Q10 coverage ledger. Unknown parameters are rejected. The returned expression
contains only scalar atoms, not an opaque whole-constraint evaluator.
"""
function core_atomic_expression(family::Symbol, xs; kwargs...)
    p=Dict{Symbol,Any}(kwargs); xs=collect(xs)
    take(key,default=nothing)=pop!(p,key,default)
    required(key)=haskey(p,key) ? pop!(p,key) : throw(ArgumentError("missing $(family) parameter $(key)"))
    result=if family===:intension
        required(:expression)
    elseif family===:extension
        tuples=required(:tuples); supports=take(:supports,true)
        all(t->length(t)==length(xs),tuples) || throw(DimensionMismatch("extension tuple width"))
        match=_aany(_aall(_ax(:eq,x,v) for (x,v) in zip(xs,t) if !(v isa AnyValue)) for t in tuples)
        supports ? match : _ax(:not,match)
    elseif family in (:allDifferent,:allEqual)
        except=take(:except,[]); form=take(:form,:list)
        op=family===:allDifferent ? :ne : :eq
        if form===:list
            _apairs(xs,op,except)
        elseif form===:matrix
            family===:allDifferent || throw(ArgumentError("allEqual matrix not a retained core form"))
            rows=_arect(required(:matrix))
            _aall(_apairs(row,op,except) for row in [rows;_acolumns(rows)])
        elseif form===:lists
            family===:allDifferent || throw(ArgumentError("allEqual lifted lists not retained"))
            rows=_arect(required(:lists))
            all(t->length(t)==length(first(rows)),except) || throw(DimensionMismatch("except tuple width"))
            excluded(row)=_aany(_aall(_ax(:eq,a,b) for (a,b) in zip(row,t)) for t in except)
            _aall(_aany([_aany(_ax(:ne,a,b) for (a,b) in zip(rows[i],rows[j])),excluded(rows[i]),excluded(rows[j])])
                for i in eachindex(rows) for j in i+1:length(rows))
        else
            throw(ArgumentError("unknown list lifting"))
        end
    elseif family===:ordered
        form=take(:form,:list)
        if form in (:lists,:matrix)
            result=core_atomic_expression(:lex,xs;form,p...)
            empty!(p)
            return result
        end
        form===:list || throw(ArgumentError("invalid ordered form"))
        op=take(:operator,:le); op in (:lt,:le,:ge,:gt) || throw(ArgumentError("invalid order"))
        lengths=take(:lengths,zeros(Int,max(0,length(xs)-1)))
        length(lengths)==max(0,length(xs)-1) || throw(DimensionMismatch("ordered lengths"))
        _aall(_ax(op,_ax(:add,xs[i],lengths[i]),xs[i+1]) for i in 1:length(xs)-1)
    elseif family===:lex
        op=take(:operator,:le); form=take(:form,:lists)
        rows=_arect(form===:matrix ? required(:matrix) : required(:lists))
        form in (:lists,:matrix) || throw(ArgumentError("invalid lex form"))
        groups=form===:matrix ? [rows,_acolumns(rows)] : [rows]
        _aall(_alex(rows[i],rows[i+1],op) for rows in groups for i in 1:length(rows)-1)
    elseif family===:precedence
        vals=required(:values); allunique(vals) || throw(ArgumentError("precedence values must be distinct"))
        covered=take(:covered,false)
        # An occurrence of the next value requires an earlier occurrence of its predecessor.
        rules=[_ax(:imp,_ax(:eq,xs[i],vals[k+1]),_amember(vals[k],xs[1:i-1])) for k in 1:length(vals)-1 for i in eachindex(xs)]
        covered && append!(rules,[_amember(v,xs) for v in vals])
        _aall(rules)
    elseif family===:sum
        coeff=take(:coefficients,ones(Int,length(xs)))
        length(coeff)==length(xs) || throw(DimensionMismatch("sum coefficients"))
        _acond(_asum(_ax(:mul,c,x) for (c,x) in zip(coeff,xs)),required(:condition))
    elseif family===:count
        _acond(_acount(xs,required(:values)),required(:condition))
    elseif family===:nValues
        except=take(:except,[])
        # Count first occurrences, not all values of an enumerated Cartesian product.
        distinct_count=_asum(_aall([_ax(:not,_amember(xs[i],except)),_ax(:not,_amember(xs[i],xs[1:i-1]))]) for i in eachindex(xs))
        _acond(distinct_count,required(:condition))
    elseif family===:cardinality
        vals=required(:values); occurs=required(:occurs); closed=take(:closed,false)
        distinct=take(:distinct_values,all(x->x isa Symbol,[xs;collect(vals);collect(occurs)]))
        distinct isa Bool || throw(ArgumentError("distinct_values must be Boolean"))
        length(vals)==length(occurs) || throw(DimensionMismatch("cardinality rows"))
        rules=Any[_acond(_acount(xs,[v]),o isa AbstractRange ? AtomicCondition(:in,o) : AtomicCondition(:eq,o)) for (v,o) in zip(vals,occurs)]
        closed && append!(rules,[_amember(x,vals) for x in xs])
        distinct && push!(rules,_apairs(vals,:ne,[]))
        _aall(rules)
    elseif family in (:minimum,:maximum,:minimumArg,:maximumArg)
        isempty(xs) && throw(ArgumentError("empty extrema scope"))
        op=family in (:minimum,:minimumArg) ? :min : :max
        extremum=_afold(op,xs,0); condition=required(:condition)
        if family in (:minimum,:maximum)
            _acond(extremum,condition)
        else
            # Java-operational Arg profile: condition on the selected INDEX, not value.
            rank=take(:rank,:any); base=take(:start_index,0)
            matches=_arank([_ax(:eq,x,extremum) for x in xs],rank)
            _aany(_aall([matches[i],_acond(big(base)+i-1,condition)]) for i in eachindex(xs))
        end
    elseif family===:element
        value=take(:value)
        condition=value===nothing ? required(:condition) : AtomicCondition(:eq,value)
        form=take(:form,:list)
        if form===:matrix
            rows=_arect(required(:matrix)); ri=required(:row_index); ci=required(:col_index)
            rb=take(:start_row_index,0); cb=take(:start_col_index,0)
            _aany(_aall([_ax(:eq,ri,big(rb)+i-1),_ax(:eq,ci,big(cb)+j-1),_acond(rows[i][j],condition)])
                for i in eachindex(rows) for j in eachindex(rows[i]))
        elseif form===:list
            index=take(:index); base=take(:start_index,0); rank=take(:rank,:any)
            matches=_arank([_acond(x,condition) for x in xs],rank)
            index===nothing ? _aany(matches) : _aany(_aall([matches[i],_ax(:eq,index,big(base)+i-1)]) for i in eachindex(xs))
        else
            throw(ArgumentError("invalid element form"))
        end
    elseif family===:channel
        form=take(:form,:single); base=take(:start_index,0)
        if form===:value
            value=required(:value)
            _aall([_ax(:eq,_asum(xs),1);[_aall([_amember(x,[0,1]),_ax(:iff,_ax(:eq,x,1),_ax(:eq,value,big(base)+i-1))]) for (i,x) in enumerate(xs)]])
        else
            ys=form===:single ? xs : form===:inverse ? collect(required(:inverse)) : throw(ArgumentError("invalid channel form"))
            otherbase=form===:single ? base : take(:inverse_start_index,0)
            length(xs)<=length(ys) || throw(ArgumentError("inverse channel requires first list no longer than second"))
            # If unequal, only the first list has a range restriction (partial inverse).
            _aall(_aany(_aall([_ax(:eq,x,big(otherbase)+j-1),_ax(:eq,y,big(base)+i-1)]) for (j,y) in enumerate(ys)) for (i,x) in enumerate(xs))
        end
    elseif family===:noOverlap
        origins=required(:origins); lengths=required(:lengths); ignored=take(:zero_ignored,true)
        o=origins isa AbstractMatrix || (!isempty(origins) && first(origins) isa AbstractVector) ? _arect(origins) : [[x] for x in origins]
        l=lengths isa AbstractMatrix || (!isempty(lengths) && first(lengths) isa AbstractVector) ? _arect(lengths) : [[x] for x in lengths]
        length(o)==length(l) && all(length(a)==length(b) for (a,b) in zip(o,l)) || throw(DimensionMismatch("noOverlap shapes"))
        all(row->length(row)==length(first(o)),o) || throw(DimensionMismatch("noOverlap dimensions"))
        signs=take(:parameter_signs,:nonnegative)
        signs in (:nonnegative,:algebraic) || throw(ArgumentError("invalid parameter_signs"))
        rules=signs===:nonnegative ? Any[_ax(:ge,x,0) for row in l for x in row] : Any[]
        for i in eachindex(o), j in i+1:length(o)
            separate=_aany(_aany([_ax(:le,_ax(:add,o[i][d],l[i][d]),o[j][d]),_ax(:le,_ax(:add,o[j][d],l[j][d]),o[i][d])]) for d in eachindex(o[i]))
            ignored && (separate=_aany([separate,_amember(0,l[i]),_amember(0,l[j])]))
            push!(rules,separate)
        end
        _aall(rules)
    elseif family===:cumulative
        origins=collect(required(:origins)); lengths=collect(required(:lengths)); heights=collect(required(:heights)); ends=take(:ends)
        condition=required(:condition); profile=take(:time_profile,:integer)
        profile in (:integer,:natural) || throw(ArgumentError("time_profile must be integer or natural"))
        signs=take(:parameter_signs,:nonnegative)
        signs in (:nonnegative,:algebraic) || throw(ArgumentError("invalid parameter_signs"))
        n=length(origins)
        length(lengths)==n && length(heights)==n && (ends===nothing || length(ends)==n) || throw(DimensionMismatch("cumulative fields"))
        finish=[_ax(:add,o,l) for (o,l) in zip(origins,lengths)]
        rules=Any[_acond(0,condition)]
        signs===:nonnegative && append!(rules,[_ax(:ge,x,0) for x in [lengths;heights]])
        ends!==nothing && append!(rules,[_ax(:eq,e,f) for (e,f) in zip(ends,finish)])
        # The actual load is constant between actual endpoints. No horizon enumeration.
        for t in (profile===:natural ? [Any[0];origins;finish] : [origins;finish])
            load=_asum(_ax(Symbol("if"),_aall([_ax(:le,origins[i],t),_ax(:lt,t,finish[i])]),heights[i],0) for i in 1:n)
            rule=_acond(load,condition)
            push!(rules,profile===:natural ? _ax(:imp,_ax(:ge,t,0),rule) : rule)
        end
        _aall(rules)
    elseif family===:binPacking
        sizes=required(:sizes); length(sizes)==length(xs) || throw(DimensionMismatch("binPacking sizes"))
        bins=required(:bins); allunique(bins) && all(x->x isa Integer,bins) || throw(ArgumentError("explicit distinct integer bins required"))
        condition=take(:condition); conditions=take(:conditions); loads=take(:loads); capacities=take(:capacities)
        limits=take(:limits)
        limits!==nothing && capacities!==nothing && throw(ArgumentError("limits and capacities are aliases"))
        limits!==nothing && (capacities=limits)
        count(!isnothing,(condition,conditions,loads,capacities))==1 || throw(ArgumentError("exactly one binPacking condition form required"))
        profile=take(:profile,:declared)
        profile in (:declared,:pdf_occupied) || throw(ArgumentError("invalid binPacking profile"))
        if condition===nothing
            per=conditions!==nothing ? collect(conditions) : [AtomicCondition(loads!==nothing ? :eq : :le,v) for v in (loads!==nothing ? loads : capacities)]
            length(per)==length(bins) || throw(DimensionMismatch("bin conditions"))
        else
            per=fill(condition,length(bins))
        end
        signs=take(:parameter_signs,:nonnegative)
        signs in (:nonnegative,:algebraic) || throw(ArgumentError("invalid parameter_signs"))
        rules=Any[_amember(x,bins) for x in xs]
        signs===:nonnegative && append!(rules,[_ax(:ge,s,0) for s in sizes])
        for (bin,c) in zip(bins,per)
            hit=[_ax(:eq,x,bin) for x in xs]
            load=_asum(_ax(:mul,s,h) for (s,h) in zip(sizes,hit))
            rule=_acond(load,c)
            (condition!==nothing || profile===:pdf_occupied) && (rule=_ax(:imp,_aany(hit),rule))
            push!(rules,rule)
        end
        _aall(rules)
    elseif family===:knapsack
        weights=required(:weights); profits=required(:profits)
        length(xs)==length(weights)==length(profits) || throw(DimensionMismatch("knapsack fields"))
        _aall([_acond(_asum(_ax(:mul,x,w) for (x,w) in zip(xs,weights)),required(:weight_condition)),
            _acond(_asum(_ax(:mul,x,p) for (x,p) in zip(xs,profits)),required(:profit_condition))])
    elseif family===:circuit
        base=take(:start_index,0); size=take(:size); n=length(xs)
        n>=2 || throw(ArgumentError("circuit requires at least two successors"))
        indices=[big(base)+i-1 for i in 1:n]
        active=[_ax(:ne,x,i) for (x,i) in zip(xs,indices)]
        # Select the first active vertex, then unroll n-1 successor steps by MUX.
        root=indices[1]
        for i in n:-1:1; root=_ax(Symbol("if"),active[i],indices[i],root); end
        orbit=Any[root]
        for step in 1:n-1
            target=0
            for (i,x) in zip(indices,xs)
                target=_ax(Symbol("if"),_ax(:eq,last(orbit),i),x,target)
            end
            push!(orbit,target)
        end
        rules=Any[[_apairs(xs,:ne,[]),_ax(:ge,_asum(active),2)];[_amember(x,indices) for x in xs]]
        append!(rules,[_ax(:imp,a,_amember(i,orbit)) for (a,i) in zip(active,indices)])
        size!==nothing && push!(rules,_ax(:eq,_asum(active),size))
        _aall(rules)
    elseif family===:instantiation
        values=required(:values); length(values)==length(xs) || throw(DimensionMismatch("instantiation values"))
        _aall(_ax(:eq,x,v) for (x,v) in zip(xs,values))
    elseif family in (:regular,:mdd)
        if family===:mdd && haskey(p,:transitions)
            layers,start,finals=_amdd_layers(required(:transitions),length(xs))
        else
            start=required(:start); finals=required(:finals)
            layers=family===:regular ? fill(collect(required(:transitions)),length(xs)) : collect(required(:layers))
        end
        length(layers)==length(xs) || throw(DimensionMismatch("one transition layer per symbol"))
        reach=Dict{Any,Any}(start=>1)
        for (x,edges) in zip(xs,layers)
            next=Dict{Any,Vector{Any}}()
            for edge in edges
                length(edge)==3 || throw(DimensionMismatch("transition must be source,label,target"))
                a,label,b=edge
                push!(get!(next,b,Any[]),_aall([get(reach,a,0),_ax(:eq,x,label)]))
            end
            reach=Dict{Any,Any}(b=>_aany(paths) for (b,paths) in next)
        end
        _aany(get(reach,f,0) for f in finals)
    elseif family===:slide
        width=required(:width); offset=take(:offset,1); circular=take(:circular,false)
        template=required(:template)
        # Explicit placeholder nodes avoid arbitrary callback code in the coverage proof.
        replace(e,window)=e isa IntensionNode ? (e.operator===:parameter ? begin
            length(e.arguments)==1 && only(e.arguments) isa Integer || throw(ArgumentError("template parameter syntax"))
            k=only(e.arguments); 0<=k<length(window) || throw(ArgumentError("template parameter outside width")); window[k+1]
        end : _ax(e.operator,(replace(a,window) for a in e.arguments)...)) : e
        scopes=slide_scopes(collect(1:length(xs)),width;offset,circular)
        _aall(replace(template,xs[s]) for s in scopes)
    else
        throw(ArgumentError("no retained atomic recipe for $(family)"))
    end
    isempty(p) || throw(ArgumentError("unused $(family) parameters: $(sort!(collect(keys(p))))"))
    result
end

"""Build a traceable atomic plan for a normalized retained core variant."""
function core_atomic_plan(books,family::Symbol;list=[b.variable for b in books],
        semantic_id::AbstractString="XCSP3-core-3.2-v4/$(family)/normalized",plan_options=NamedTuple(),kwargs...)
    # Precedence's omitted values mean the ordered union of possible source values.
    if family===:precedence && !haskey(kwargs,:values)
        byname=Dict(b.variable=>b for b in books)
        all(x->x isa Symbol && haskey(byname,x),list) || throw(ArgumentError("implicit precedence values require variable list"))
        values=sort!(unique(BigInt[v for x in list for v in byname[x].values]))
        expression=core_atomic_expression(family,list;values,kwargs...)
    elseif family===:binPacking && !haskey(kwargs,:bins)
        if haskey(kwargs,:condition)
            byname=Dict(b.variable=>b for b in books)
            all(x->x isa Symbol && haskey(byname,x),list) || throw(ArgumentError("implicit bins require variable list"))
            bins=sort!(unique(BigInt[v for x in list for v in byname[x].values]))
        else
            field=findfirst(k->haskey(kwargs,k),(:limits,:capacities,:loads,:conditions))
            field===nothing && throw(ArgumentError("binPacking needs a condition or declared bin list"))
            bins=collect(0:length(kwargs[(:limits,:capacities,:loads,:conditions)[field]])-1)
        end
        expression=core_atomic_expression(family,list;bins,kwargs...)
    else
        expression=core_atomic_expression(family,list;kwargs...)
    end
    atomic_plan(books,expression;semantic_id,merge((partial=:defined,),plan_options)...)
end
