"""Lower partial arithmetic to value/definedness wires, without a global truth table.

Invalid arithmetic is NOT given an arbitrary satisfying value: the dummy value 0
is accompanied by a false definedness bit. Only the selected branch of `if` needs
to be defined. Other operators are strict. A root is accepted only if defined and
true. On every mathematically defined integer expression this preserves semantics.
"""
function _defined_expression(expression; max_depth=1000, max_nodes=100000)
    memo=IdDict{IntensionNode,Tuple{Any,Any}}(); visited=Ref(0)
    N=IntensionNode
    both(xs)=begin
        ys=Any[x for x in xs if x!=1]
        isempty(ys) ? 1 : length(ys)==1 ? only(ys) : N(:and,ys...)
    end
    function visit(e,depth)
        depth<=max_depth || throw(CompilationLimit(:definedness_depth,big(depth),big(max_depth)))
        e isa IntensionNode || return (e,1)
        haskey(memo,e) && return memo[e]
        visited[]+=1
        visited[]<=max_nodes || throw(CompilationLimit(:definedness_nodes,big(visited[]),big(max_nodes)))
        op=e.operator
        op===:set && return (e,1)
        parts=[visit(a,depth+1) for a in e.arguments]
        values=first.(parts); defined=last.(parts)
        valid=both(defined)
        value=if op in (:div,:mod,:pow)
            length(values)==2 || throw(ArgumentError("partial arithmetic arity"))
            guard=op===:pow ? N(:pow_defined,values...) : N(:ne,values[2],0)
            valid=both([valid,guard])
            N(op===:div ? :div_total : op===:mod ? :mod_total : :pow_total,values...)
        elseif op===Symbol("if")
            length(values)==3 || throw(ArgumentError("if arity"))
            valid=both([defined[1],defined[2]==defined[3] ? defined[2] : N(Symbol("if"),values[1],defined[2],defined[3])])
            N(op,values...)
        else
            N(op,values...)
        end
        memo[e]=(value,valid)
    end
    value,defined=visit(expression,1)
    defined==1 ? value : N(:and,defined,value)
end
