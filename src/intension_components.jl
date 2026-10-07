"""Explicit pure expression node; integers/Symbols are constants/variable references.

No Julia expression evaluation or XML parsing. Repeated subexpressions are interned
into a bounded topological DAG before finite truth/Rosenberg quadratization.
"""
struct IntensionNode
    operator::Symbol
    arguments::Tuple
end
IntensionNode(op::Symbol,args...)=IntensionNode(op,args)

function _intension_apply(op,a,max_value_bits)
    boolean(x)=x in (0,1) ? x==1 : throw(ArgumentError("logical operand outside 0:1"))
    result=if op===:add
        sum(a)
    elseif op===:mul
        foldl(a;init=big(1)) do x,y
            ndigits(abs(x);base=2)+ndigits(abs(y);base=2)<=max_value_bits+1 ||
                throw(CompilationLimit(:intension_value_bits,big(ndigits(abs(x);base=2)+ndigits(abs(y);base=2)-1),big(max_value_bits)))
            x*y
        end
    elseif op===:sub
        a[1]-a[2]
    elseif op===:neg
        -a[1]
    elseif op===:abs
        abs(a[1])
    elseif op===:sqr
        _intension_apply(:mul,[a[1],a[1]],max_value_bits)
    elseif op===:dist
        abs(a[1]-a[2])
    elseif op===:min
        minimum(a)
    elseif op===:max
        maximum(a)
    elseif op===:eq
        all(==(first(a)),a)
    elseif op===:ne
        a[1]!=a[2]
    elseif op===:lt
        a[1]<a[2]
    elseif op===:le
        a[1]<=a[2]
    elseif op===:gt
        a[1]>a[2]
    elseif op===:ge
        a[1]>=a[2]
    elseif op===:not
        !boolean(a[1])
    elseif op in (:and,:or,:xor,:iff,:imp)
        b=boolean.(a)
        op===:and ? all(b) : op===:or ? any(b) : op===:xor ? isodd(count(identity,b)) :
            op===:iff ? all(==(first(b)),b) : (!b[1] || b[2])
    elseif op===Symbol("if")
        boolean(a[1]) ? a[2] : a[3]
    else
        throw(ArgumentError("unsupported intension operator $(op)"))
    end
    bits=ndigits(abs(big(result));base=2)
    bits<=max_value_bits || throw(CompilationLimit(:intension_value_bits,big(bits),big(max_value_bits)))
    return big(result)
end

"""Finite exact compilation of a pure bounded integer/Boolean expression DAG.

Supports add/mul/sub/neg/abs/sqr/dist/min/max, comparisons and Boolean operations/if.
div/mod/pow/sets are explicitly unsupported; no implicit convention on undefined
arithmetic. All supported operations are total on their checked domains. This is
truth-table + AND quadratization, not a scalable arithmetic gate library.
"""
function intension_component(books,expression;max_nodes::Integer=1000,max_edges::Integer=10000,max_depth::Integer=64,
        max_value_bits::Integer=256,kwargs...)
    max_nodes>0 && max_edges>=0 && max_depth>0 && max_value_bits>0 || throw(ArgumentError("invalid expression budgets"))
    books,_,positions=_finite_scope(books,Symbol[])
    nodes=Tuple{Symbol,Any}[]; intern=Dict{Any,Int}(); heights=Int[]
    seen=IdDict{IntensionNode,Int}(); edges=big(0)
    unary=(:neg,:abs,:sqr,:not); binary=(:sub,:dist,:ne,:lt,:le,:gt,:ge,:imp)
    variadic=(:add,:mul,:min,:max,:eq,:and,:or,:xor,:iff)
    function visit(e,depth)
        depth<=max_depth || throw(CompilationLimit(:intension_depth,big(depth),big(max_depth)))
        e isa IntensionNode && haskey(seen,e) && return seen[e]
        key=if e isa Integer
            ndigits(abs(big(e));base=2)<=max_value_bits || throw(CompilationLimit(:intension_value_bits,big(ndigits(abs(big(e));base=2)),big(max_value_bits)))
            (:constant,big(e))
        elseif e isa Symbol
            haskey(positions,e) || throw(ArgumentError("unknown expression variable $(e)"))
            bits=maximum(v->ndigits(abs(big(v));base=2),_book_values(books[positions[e]]))
            bits<=max_value_bits || throw(CompilationLimit(:intension_value_bits,big(bits),big(max_value_bits)))
            (:variable,positions[e])
        elseif e isa IntensionNode
            op=e.operator; n=length(e.arguments)
            valid=op in unary ? n==1 : op in binary ? n==2 : op in variadic ? n>=2 : op===Symbol("if") ? n==3 : false
            valid || throw(ArgumentError("unsupported operator or arity: $(op)/$(n)"))
            edges+=n
            edges<=max_edges || throw(CompilationLimit(:intension_edges,edges,big(max_edges)))
            (op,Tuple(visit(a,depth+1) for a in e.arguments))
        else
            throw(ArgumentError("expression requires IntensionNode, integer or Symbol"))
        end
        if haskey(intern,key)
            e isa IntensionNode && (seen[e]=intern[key])
            return intern[key]
        end
        length(nodes)<max_nodes || throw(CompilationLimit(:intension_nodes,big(length(nodes)+1),big(max_nodes)))
        push!(nodes,key); intern[key]=length(nodes)
        e isa IntensionNode && (seen[e]=length(nodes))
        return length(nodes)
    end
    root=visit(expression,1)
    for (op,data) in nodes
        h=op in (:constant,:variable) ? 1 : 1+maximum(heights[i] for i in data)
        h<=max_depth || throw(CompilationLimit(:intension_depth,big(h),big(max_depth)))
        push!(heights,h)
    end
    function predicate(v)
        values=BigInt[]
        for (op,data) in nodes
            push!(values,op===:constant ? data : op===:variable ? big(v[data]) :
                _intension_apply(op,[values[i] for i in data],max_value_bits))
        end
        values[root] in (0,1) || throw(ArgumentError("intension root must be Boolean"))
        return values[root]==1
    end
    return truth_component(books,predicate;oracle_id="Q7/intension/finite-DAG-v1/$(nodes)/root=$(root)",kwargs...)
end
