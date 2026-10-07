"""A finite, topologically ordered operation. Tables belong to ONE atom, never a constraint.

`rows` contains input-domain indices followed by the output-domain index. Fields
are owned by the plan and must be treated as read-only, like QUBOComponent arrays.
"""
struct AtomicNode
    operator::Symbol
    inputs::Vector{Int}
    data::Any
    domain::Vector{BigInt}
    rows::Vector{Vector{Int}}
end

struct AtomicPlan
    codebooks::Vector{AbstractCodebook}
    nodes::Vector{AtomicNode}
    root::Int
    semantic_id::String
end

const _ATOMIC_UNARY = (:neg,:abs,:sqr,:not)
const _ATOMIC_BINARY = (:add,:sub,:mul,:div,:mod,:pow,:dist,:min,:max,
    :eq,:ne,:lt,:le,:ge,:gt,:and,:or,:xor,:iff,:imp,
    :div_total,:mod_total,:pow_total,:pow_defined)

function _atomic_apply(op, a, max_value_bits)
    if op===:pow_defined
        return BigInt(a[2]>=0 || abs(a[1])==1)
    elseif op in (:div_total,:mod_total)
        return iszero(a[2]) ? big(0) : _atomic_apply(op===:div_total ? :div : :mod,a,max_value_bits)
    elseif op===:pow_total
        if a[2]<0
            return a[1]==1 ? big(1) : a[1]==-1 ? big(isodd(a[2]) ? -1 : 1) : big(0)
        end
        return _atomic_apply(:pow,a,max_value_bits)
    end
    result = if op===:div || op===:mod
        iszero(a[2]) && throw(ArgumentError("undefined $(op): zero denominator; no implicit XCSP3 convention"))
        op===:div ? div(a[1],a[2],RoundToZero) : rem(a[1],a[2])
    elseif op===:pow
        a[2]>=0 || throw(ArgumentError("negative exponent outside the retained integer-power profile"))
        # Bound growth BEFORE exponentiation; special bases do not require Int conversion.
        a[2]==0 ? big(1) : a[1]==0 ? big(0) : a[1]==1 ? big(1) : a[1]==-1 ? big(isodd(a[2]) ? -1 : 1) : begin
            lower=(ndigits(abs(a[1]);base=2)-1)*a[2]+1
            lower<=max_value_bits || throw(CompilationLimit(:atomic_value_bits,lower,big(max_value_bits)))
            a[1]^Int(a[2])
        end
    else
        _intension_apply(op,a,max_value_bits)
    end
    ndigits(abs(result);base=2)<=max_value_bits ||
        throw(CompilationLimit(:atomic_value_bits,big(ndigits(abs(result);base=2)),big(max_value_bits)))
    return BigInt(result)
end

"""Lower an explicit expression into unary/binary atoms and ternary `if`.

Variadic reductions are expanded, not tabulated globally. Only local Cartesian
products are inspected. Partial arithmetic is rejected on the inferred domains;
this may conservatively reject a correlated expression, never certify it falsely.
"""
function atomic_plan(books, expression; semantic_id::AbstractString,
        max_nodes::Integer=100000, max_local_rows::Integer=100000,
        max_total_rows::Integer=1000000, max_value_bits::Integer=256,
        max_depth::Integer=1000, partial::Symbol=:reject)
    partial in (:reject,:defined) || throw(ArgumentError("partial must be reject or defined"))
    if partial===:defined
        expression=_defined_expression(expression;max_depth,max_nodes)
        semantic_id=String(semantic_id)*"/defined-integer-strict-if-selected-v1"
    end
    all(x->x>0,(max_nodes,max_local_rows,max_total_rows,max_value_bits,max_depth)) || throw(ArgumentError("positive atomic budgets required"))
    books=AbstractCodebook[deepcopy(b) for b in books]
    allunique(b.variable for b in books) || throw(ArgumentError("duplicate source variable"))
    all(b->all(v->v isa Integer,b.values),books) || throw(ArgumentError("integer source domains required"))
    byname=Dict(b.variable=>i for (i,b) in enumerate(books))
    nodes=AtomicNode[]; heights=Int[]; intern=Dict{Any,Int}(); seen=IdDict{IntensionNode,Int}(); total=big(0)
    function insert(op, inputs, data)
        key=(op,Tuple(inputs),data)
        haskey(intern,key) && return intern[key]
        length(nodes)<max_nodes || throw(CompilationLimit(:atomic_nodes,big(length(nodes)+1),big(max_nodes)))
        height=isempty(inputs) ? 1 : 1+maximum(heights[i] for i in inputs)
        height<=max_depth || throw(CompilationLimit(:atomic_depth,big(height),big(max_depth)))
        rows=Vector{Int}[]
        domain=if op===:constant
            BigInt[data]
        elseif op===:variable
            BigInt.(books[data].values)
        else
            count=prod((big(length(nodes[i].domain)) for i in inputs);init=big(1))
            count<=max_local_rows || throw(CompilationLimit(:atomic_local_rows,count,big(max_local_rows)))
            total+=count
            total<=max_total_rows || throw(CompilationLimit(:atomic_total_rows,total,big(max_total_rows)))
            values=BigInt[]
            for indices in Iterators.product((eachindex(nodes[i].domain) for i in inputs)...)
                args=BigInt[nodes[i].domain[k] for (i,k) in zip(inputs,indices)]
                value=_atomic_apply(op,args,max_value_bits)
                k=findfirst(==(value),values)
                if k===nothing; push!(values,value); k=length(values); end
                push!(rows,[Int[indices...];k])
            end
            values
        end
        all(v->ndigits(abs(v);base=2)<=max_value_bits,domain) || throw(ArgumentError("source/constant exceeds atomic value budget"))
        push!(nodes,AtomicNode(op,collect(inputs),data,domain,rows))
        push!(heights,height)
        intern[key]=length(nodes)
    end
    function visit(e,depth)
        depth<=max_depth || throw(CompilationLimit(:atomic_depth,big(depth),big(max_depth)))
        e isa Integer && return insert(:constant,Int[],big(e))
        if e isa Symbol
            haskey(byname,e) || throw(ArgumentError("unknown atomic variable $(e)"))
            return insert(:variable,Int[],byname[e])
        end
        e isa IntensionNode || throw(ArgumentError("atomic expressions require IntensionNode, integer or Symbol"))
        haskey(seen,e) && return seen[e]
        op=e.operator; args=e.arguments
        if op in (:in,:notin)
            length(args)==2 || throw(ArgumentError("membership arity"))
            set=args[2]
            set isa IntensionNode && set.operator===:set || throw(ArgumentError("membership requires constant set node"))
            all(x->x isa Integer,set.arguments) || throw(ArgumentError("intension sets must be constant"))
            x=visit(args[1],depth+1)
            ids=[insert(:eq,[x,visit(v,depth+1)],nothing) for v in unique(set.arguments)]
            out=isempty(ids) ? insert(:constant,Int[],big(0)) : foldl((a,b)->insert(:or,[a,b],nothing),ids)
            seen[e]=op===:in ? out : insert(:not,[out],nothing)
            return seen[e]
        end
        n=length(args)
        if op in (:add,:mul,:min,:max,:and,:or,:xor,:eq,:iff)
            n>=2 || throw(ArgumentError("variadic atomic operator needs at least two operands"))
            ids=[visit(a,depth+1) for a in args]
            if op in (:eq,:iff)
                seen[e]=foldl((a,b)->insert(:and,[a,b],nothing),[insert(op,[ids[1],i],nothing) for i in ids[2:end]])
                return seen[e]
            end
            seen[e]=foldl((a,b)->insert(op,[a,b],nothing),ids)
            return seen[e]
        end
        valid=op in _ATOMIC_UNARY ? n==1 : op in _ATOMIC_BINARY ? n==2 : op===Symbol("if") ? n==3 : false
        valid || throw(ArgumentError("unsupported atomic operator/arity $(op)/$(n)"))
        seen[e]=insert(op,[visit(a,depth+1) for a in args],nothing)
        return seen[e]
    end
    # Even unused source variables retain their validity in the final QUBO.
    for b in books; visit(b.variable,1); end
    root=visit(expression,1)
    all(in((0,1)),nodes[root].domain) || throw(ArgumentError("atomic root must be Boolean"))
    return AtomicPlan(books,nodes,root,String(semantic_id))
end

"""Evaluate the explicit atomic graph (a diagnostic, NOT a post-QUBO repair)."""
function atomic_values(plan::AtomicPlan, assignment)
    values=BigInt[]
    for n in plan.nodes
        v=n.operator===:constant ? n.data : n.operator===:variable ? assignment[plan.codebooks[n.data].variable] :
            _atomic_apply(n.operator,BigInt[values[i] for i in n.inputs],typemax(Int)÷4)
        v in n.domain || throw(ArgumentError("assignment outside atomic domains"))
        push!(values,BigInt(v))
    end
    values
end

struct AtomicCompilation
    plan::AtomicPlan
    component::QUBOComponent{BigInt}
    factors::Vector{QUBOComponent{BigInt}}
    labels::Vector{String}
    wires::Vector{Vector{BitID}}
    selectors::Dict{Int,Vector{BitID}}
    input_codes::Dict{Int,Vector{BitVector}}
end

"""Recheck the plan's local truth relations before relying on the compositional proof.

This also protects callers constructing AtomicPlan directly or altering its arrays.
It does not attest that a user-supplied recipe matches an external specification.
"""
function validate_atomic_plan(plan::AtomicPlan;max_rows::Integer=1000000,max_value_bits::Integer=256)
    1<=plan.root<=length(plan.nodes) || throw(ArgumentError("invalid atomic root"))
    allunique(b.variable for b in plan.codebooks) || throw(ArgumentError("duplicate input identity"))
    covered=Int[]; total=big(0)
    for (i,n) in enumerate(plan.nodes)
        !isempty(n.domain) && allunique(n.domain) || throw(ArgumentError("invalid atomic domain"))
        if n.operator in (:variable,:constant)
            isempty(n.inputs) && isempty(n.rows) || throw(ArgumentError("leaf has incoming relation"))
            if n.operator===:variable
                n.data isa Int && 1<=n.data<=length(plan.codebooks) || throw(ArgumentError("invalid source reference"))
                n.domain==BigInt.(plan.codebooks[n.data].values) || throw(ArgumentError("source domain mismatch"))
                push!(covered,n.data)
            else
                n.domain==[n.data] || throw(ArgumentError("constant domain mismatch"))
            end
            continue
        end
        arity=length(n.inputs)
        (n.operator in _ATOMIC_UNARY ? arity==1 : n.operator in _ATOMIC_BINARY ? arity==2 : n.operator===Symbol("if") && arity==3) ||
            throw(ArgumentError("invalid atomic operator/arity"))
        all(j->1<=j<i,n.inputs) || throw(ArgumentError("atomic plan is not topological"))
        expected=prod((big(length(plan.nodes[j].domain)) for j in n.inputs);init=big(1))
        total+=expected
        total<=max_rows || throw(CompilationLimit(:atomic_validation_rows,total,big(max_rows)))
        length(n.rows)==expected || throw(ArgumentError("incomplete atomic truth relation"))
        seen=Set{Tuple}()
        for row in n.rows
            length(row)==arity+1 || throw(ArgumentError("invalid atomic row width"))
            all(1<=k<=length(plan.nodes[j].domain) for (j,k) in zip(n.inputs,row[1:end-1])) && 1<=last(row)<=length(n.domain) ||
                throw(ArgumentError("atomic row index outside domain"))
            key=Tuple(row[1:end-1]); key in seen && throw(ArgumentError("duplicate atomic input row")); push!(seen,key)
            value=_atomic_apply(n.operator,BigInt[plan.nodes[j].domain[k] for (j,k) in zip(n.inputs,key)],max_value_bits)
            n.domain[last(row)]==value || throw(ArgumentError("incorrect atomic truth row"))
        end
    end
    sort!(covered)==collect(eachindex(plan.codebooks)) || throw(ArgumentError("missing/duplicate source nodes"))
    all(in((0,1)),plan.nodes[plan.root].domain) || throw(ArgumentError("non-Boolean atomic root"))
    true
end

"""Compile a graph to an explicit sum of nonnegative square matrices.

Every wire and every local operation tuple has an exactly-one selector. Equality
of marginal selectors links the ports. Source code bridges preserve ALL redundant
codes. No whole-constraint truth table, hidden Boolean check, or optimization solver
is used. Internal wires are one-hot; the chosen PRIMARY encoding is unchanged.
"""
function compile_atomic(plan::AtomicPlan; namespace::Union{Nothing,Symbol}=nothing,
        max_terms::Integer=2000000,max_bits::Integer=100000,max_code_states::Integer=65536,
        max_validation_rows::Integer=1000000,max_value_bits::Integer=256)
    max_terms>0 && max_bits>0 && max_code_states>0 || throw(ArgumentError("positive compilation budgets required"))
    validate_atomic_plan(plan;max_rows=max_validation_rows,max_value_bits)
    if namespace===nothing
        namespace=:atomic; suffix=0
        while any(b->b.variable===namespace || startswith(string(b.variable),string(namespace)*"/"),plan.codebooks)
            suffix+=1; namespace=Symbol("atomic_",suffix)
        end
    end
    any(b->b.variable===namespace || startswith(string(b.variable),string(namespace)*"/"),plan.codebooks) &&
        throw(ArgumentError("reserved atomic namespace collides with source"))
    factors=QUBOComponent{BigInt}[]; labels=String[]; total=big(0); allbits=BitID[b for book in plan.codebooks for b in book.bits]
    function allocate(owner,n)
        big(length(allbits))+n<=max_bits || throw(CompilationLimit(:atomic_bits,big(length(allbits))+n,big(max_bits)))
        bits=[BitID(Symbol(namespace,"/",owner),k;role=:semantic_auxiliary) for k in 1:n]
        append!(allbits,bits); bits
    end
    function square(label,offset,terms)
        # Each retained component has its own small bit order; the final sum shares identities.
        bits=sort!(unique(BitID[first(t) for t in terms])); ids=Dict(b=>i for (i,b) in enumerate(bits))
        builder=_ArithmeticBuilder(bits,Int(max_terms))
        _arithmetic_square!(builder,big(offset),[(ids[b],big(w)) for (b,w) in terms])
        total+=length(builder.linear)+length(builder.quadratic)
        total<=max_terms || throw(CompilationLimit(:atomic_terms,total,big(max_terms)))
        push!(factors,QUBOComponent(bits;linear=builder.linear,quadratic=builder.quadratic,offset=builder.offset,provenance=label))
        push!(labels,label)
    end
    exactly(label,bits)=square(label,-1,[(b,1) for b in bits])
    wires=[allocate("wire/$(i)",length(n.domain)) for (i,n) in enumerate(plan.nodes)]
    selectors=Dict{Int,Vector{BitID}}(); input_codes=Dict{Int,Vector{BitVector}}()
    for (i,n) in enumerate(plan.nodes)
        exactly("node $(i) $(n.operator): value validity",wires[i])
        n.operator===:constant && continue
        if n.operator===:variable
            book=plan.codebooks[n.data]
            mode=book isa StructuredCodebook ? book.encoding :
                _compatible(book,codebook(book.variable,book.values;encoding=:one_hot)) ? :one_hot :
                _compatible(book,codebook(book.variable,book.values;encoding=:domain_wall)) ? :domain_wall : :explicit
            if mode in (:one_hot,:zero_one_hot,:domain_wall,:native)
                for (j,bit) in enumerate(book.bits)
                    ks=mode===:one_hot ? [j] : mode===:zero_one_hot ? [j+1] : mode===:domain_wall ? collect(j+1:length(n.domain)) :
                        length(n.domain)==2 ? [2] : Int[]
                    square("input $(book.variable): $(mode) bit $(j)",0,[(bit,1);[(wires[i][k],-1) for k in ks]])
                end
            else
                count=big(2)^length(book.bits)
                count<=max_code_states || throw(CompilationLimit(:atomic_source_code_states,count,big(max_code_states)))
                codes=BitVector[]; indices=Int[]
                for mask in 0:Int(count)-1
                    code=BitVector(!iszero(mask & (1<<(j-1))) for j in eachindex(book.bits))
                    k=_value_index(book,code)
                    k==0 && continue
                    push!(codes,code); push!(indices,k)
                end
                selected=allocate("codes/$(i)",length(codes)); selectors[i]=selected; input_codes[i]=codes
                exactly("input $(book.variable): valid-code selection",selected)
                for (j,bit) in enumerate(book.bits)
                    square("input $(book.variable): code bit $(j)",0,[(bit,1);[(selected[k],-1) for k in eachindex(codes) if codes[k][j]]])
                end
                for k in eachindex(n.domain)
                    square("input $(book.variable): decoded value $(k)",0,[(wires[i][k],1);[(selected[j],-1) for j in eachindex(codes) if indices[j]==k]])
                end
            end
        else
            selected=allocate("tuples/$(i)",length(n.rows)); selectors[i]=selected
            exactly("node $(i) $(n.operator): local tuple",selected)
            for (port,wire) in enumerate([n.inputs;i]), k in eachindex(wires[wire])
                square("node $(i) $(n.operator): port $(port) value $(k)",0,
                    [(wires[wire][k],1);[(selected[r],-1) for r in eachindex(n.rows) if n.rows[r][port]==k]])
            end
        end
    end
    accepted=findfirst(==(1),plan.nodes[plan.root].domain)
    square("root acceptance",1,accepted===nothing ? Tuple{BitID,Int}[] : [(wires[plan.root][accepted],-1)])
    ids=Dict(b=>i for (i,b) in enumerate(allbits))
    q=QUBOComponent(allbits;
        linear=[(ids[f.bits[i]],v) for f in factors for (i,v) in _linear_terms(f)],
        quadratic=[(ids[f.bits[i]],ids[f.bits[j]],v) for f in factors for (i,j,v) in _quadratic_terms(f)],
        offset=sum(f.offset for f in factors),codebooks=plan.codebooks,
        auxiliary_meanings=Dict(b=>"atomic wire/code/tuple $(b.owner)" for b in allbits if b.role!==:primary),
        applicability="finite integer atomic graph; exact zero set, integer gap >= 1; total retained operations",
        provenance="atomic-square-v1; $(plan.semantic_id); $(length(plan.nodes)) nodes; $(length(factors)) factors")
    AtomicCompilation(plan,q,factors,labels,wires,selectors,input_codes)
end

"""Construct the unique semantic execution and a matching tuple witness for a source code.

For a valid source assignment its energy is 0 iff accepted, otherwise 1. The
compiler proof, not this witness alone, rules out a different false zero.
"""
function atomic_witness(c::AtomicCompilation, source_codes::AbstractDict)
    assignment=Dict{Symbol,BigInt}(); bitvalues=Dict{BitID,Bool}()
    for book in c.plan.codebooks
        code=source_codes[book.variable]; result=decode_code(book,code)
        result.valid || throw(ArgumentError("invalid source code"))
        assignment[book.variable]=BigInt(result.value)
        for (b,v) in zip(book.bits,code); bitvalues[b]=Bool(v); end
    end
    values=atomic_values(c.plan,assignment)
    for (i,n) in enumerate(c.plan.nodes)
        k=findfirst(==(values[i]),n.domain)
        for (j,b) in enumerate(c.wires[i]); bitvalues[b]=(j==k); end
        if haskey(c.selectors,i)
            chosen=if n.operator===:variable
                findfirst(==(BitVector(source_codes[c.plan.codebooks[n.data].variable])),c.input_codes[i])
            else
                row=[Int[findfirst(==(values[v]),c.plan.nodes[v].domain) for v in n.inputs];k]
                findfirst(==(row),n.rows)
            end
            chosen===nothing && error("missing atomic witness row")
            for (j,b) in enumerate(c.selectors[i]); bitvalues[b]=(j==chosen); end
        end
    end
    BitVector(bitvalues[b] for b in c.component.bits)
end
