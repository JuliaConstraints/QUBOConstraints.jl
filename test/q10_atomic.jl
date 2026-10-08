@testitem "Q10 local atoms, projected minima and source-code bridges" tags=[:q10] begin
    using LinearAlgebra
    N=IntensionNode
    for mode in (:native,:one_hot,:domain_wall), op in (:and,:or,:xor,:iff,:imp,:eq,:ne,:lt,:le)
        books=[structured_codebook(s,0:1;encoding=mode) for s in (:x,:y)]
        p=atomic_plan(books,N(op,:x,:y);semantic_id="q10/Boolean-atom/$(op)")
        c=compile_atomic(p)
        expected(v)=op===:and ? v[1]==1 && v[2]==1 : op===:or ? v[1]==1 || v[2]==1 :
            op in (:xor,:ne) ? v[1]!=v[2] : op in (:iff,:eq) ? v[1]==v[2] : op===:imp ? v[1]==0 || v[2]==1 :
            op===:lt ? v[1]<v[2] : v[1]<=v[2]
        @test exhaustive_check(c.component,expected;oracle_id="q10/independent-Boolean",profile=:indicator_exact,max_states=1<<20).status===:pass
        for x in 0:1,y in 0:1
            z=atomic_witness(c,Dict(:x=>encode_code(books[1],x),:y=>encode_code(books[2],y)))
            @test energy(c.component,z)==Int(!expected([x,y]))
            m=dense_qubo(c.component)
            @test istriu(m.matrix)
            @test dot(z,m.matrix*z)+m.offset==energy(c.component,z)
        end
    end
    # Preserve redundant codes and expose Fibonacci/binary differences beyond 0:2.
    for mode in (:native,:one_hot,:zero_one_hot,:domain_wall,:unary,:binary,:bounded_binary,:gray,:arithmetic,:bounded_coefficient,:fibonacci)
        domain=mode===:native ? [-3,7] : [-3,0,2,7,11,20]
        book=structured_codebook(:x,domain;encoding=mode)
        plan=atomic_plan([book],N(:in,:x,N(:set,-3,7));semantic_id="q10/encoding/$(mode)")
        c=compile_atomic(plan)
        for mask in 0:(1<<length(book.bits))-1
            code=BitVector(!iszero(mask & (1<<(j-1))) for j in eachindex(book.bits))
            d=decode_code(book,code)
            if d.valid
                z=atomic_witness(c,Dict(:x=>code))
                @test energy(c.component,z)==Int(!(d.value in (-3,7)))
            else
                @test_throws ArgumentError atomic_witness(c,Dict(:x=>code))
            end
        end
    end
    # Source-only plans make invalid-code minimization small enough to enumerate.
    for mode in (:one_hot,:zero_one_hot,:domain_wall,:unary,:binary,:bounded_binary,:gray,:arithmetic,:bounded_coefficient,:fibonacci)
        b=structured_codebook(:x,[0,1];encoding=mode)
        c=compile_atomic(atomic_plan([b],:x;semantic_id="q10/source-root"))
        @test exhaustive_check(c.component,v->only(v)==1;oracle_id="q10/source-validity",profile=:indicator_exact).status===:pass
    end
    books=[structured_codebook(:x,0:1)]
    p=atomic_plan(books,N(:not,:x);semantic_id="q10/mutation")
    @test validate_atomic_plan(p)
    p.nodes[end].rows[1][end]=p.nodes[end].rows[2][end]
    @test_throws ArgumentError compile_atomic(p)
    @test_throws ArgumentError atomic_plan(books,N(:div,1,:x);semantic_id="q10/undefined")
    @test_throws ArgumentError atomic_plan(books,N(:pow,2,-1);semantic_id="q10/negative-power")
    @test_throws CompilationLimit atomic_plan(books,N(:eq,N(:pow,2,1000000),0);semantic_id="q10/power-budget")
    @test_throws ArgumentError core_atomic_plan(books,:sum;condition=AtomicCondition(:eq,1),ignored_parameter=42)
    @test_throws ArgumentError core_atomic_plan(books,:unknown)
    @test_throws ArgumentError core_atomic_plan(books,:cumulative;origins=[:x],lengths=[1],heights=[1],condition=AtomicCondition(:le,1),time_profile=:naturals)
    @test_throws CompilationLimit compile_atomic(atomic_plan(books,:x;semantic_id="q10/budget");max_bits=1)
    shared=N(:not,N(:not,:x))
    @test_throws CompilationLimit atomic_plan(books,N(:and,shared,N(:not,shared));semantic_id="q10/shared-depth",max_depth=4)
    @test_throws ArgumentError compile_atomic(atomic_plan([structured_codebook(Symbol("atomic/x"),0:1)],1;semantic_id="q10/namespace");namespace=:atomic)
    @test atomic_values(atomic_plan([],N(:eq,N(:pow,0,0),1);semantic_id="q10/power-zero"),Dict())[end]==1
end

@testitem "Atomic validation fixed-arity truth rows and rejection guards" tags=[:q10] begin
    N=IntensionNode
    unary=Dict(:neg=>(x->-x),:abs=>abs,:sqr=>(x->x*x),:not=>(x->Int(x==0)))
    binary=Dict(
        :add=>+, :sub=>-, :mul=>*, :div=>((x,y)->div(x,y,RoundToZero)),
        :mod=>rem, :pow=>^, :dist=>((x,y)->abs(x-y)), :min=>min, :max=>max,
        :eq=>==, :ne=>!=, :lt=><, :le=><=, :ge=>>=, :gt=>>,
        :and=>((x,y)->x==1 && y==1), :or=>((x,y)->x==1 || y==1),
        :xor=>((x,y)->(x==1)!=(y==1)), :iff=>((x,y)->x==y),
        :imp=>((x,y)->x==0 || y==1),
        :div_total=>((x,y)->y==0 ? 0 : div(x,y,RoundToZero)),
        :mod_total=>((x,y)->y==0 ? 0 : rem(x,y)),
        :pow_total=>((x,y)->y>=0 ? x^y : x==1 ? 1 : x==-1 ? (isodd(y) ? -1 : 1) : 0),
        :pow_defined=>((x,y)->y>=0 || abs(x)==1))
    # Expected results use small native integers and independent primitive
    # operations, rather than either atomic evaluator or generated QUBO energy.
    for (op,oracle) in unary
        domain=op===:not ? [0,1] : [-2,0,3]
        p=atomic_plan([structured_codebook(:x,domain)],N(:eq,N(op,:x),0);semantic_id="validator/unary/$(op)")
        @test validate_atomic_plan(p)
        node=p.nodes[findfirst(n->n.operator===op,p.nodes)]
        for row in node.rows
            @test node.domain[last(row)]==oracle(domain[first(row)])
        end
    end
    for (op,oracle) in binary
        logical=op in (:and,:or,:xor,:iff,:imp)
        x=logical ? [0,1] : [-2,0,3]
        y=logical ? [0,1] : op in (:div,:mod) ? [-2,1,3] :
            op===:pow ? [0,1,3] : op in (:pow_total,:pow_defined) ? [-3,-1,0,2] : [-1,0,2]
        p=atomic_plan([structured_codebook(:x,x),structured_codebook(:y,y)],
            N(:eq,N(op,:x,:y),0);semantic_id="validator/binary/$(op)")
        @test validate_atomic_plan(p)
        node=p.nodes[findfirst(n->n.operator===op,p.nodes)]
        for row in node.rows
            @test node.domain[last(row)]==oracle(x[row[1]],y[row[2]])
        end
    end
    books=[structured_codebook(:c,0:1),structured_codebook(:x,[-2,3]),structured_codebook(:y,[-1,4])]
    p=atomic_plan(books,N(:eq,N(Symbol("if"),:c,:x,:y),0);semantic_id="validator/ternary")
    @test validate_atomic_plan(p)
    node=p.nodes[findfirst(n->n.operator===Symbol("if"),p.nodes)]
    for row in node.rows
        @test node.domain[last(row)]==(row[1]==2 ? books[2].values[row[2]] : books[3].values[row[3]])
    end
    repeated=atomic_plan([structured_codebook(:x,[-2,3])],N(:eq,N(:add,:x,:x),0);semantic_id="validator/repeated-input")
    @test validate_atomic_plan(repeated)
    add=repeated.nodes[2]
    @test add.inputs==[1,1] && length(add.rows)==4
    for row in add.rows
        @test add.domain[last(row)]==sum(repeated.nodes[1].domain[row[k]] for k in 1:2)
    end

    # Mutate valid plans at every accepted arity. The complete row count must
    # never allow duplicated, out-of-range, misshaped or incorrect rows through.
    for expression in (N(:not,:x),N(:eq,:x,:y),N(Symbol("if"),:x,:y,:z))
        original=atomic_plan([structured_codebook(s,0:1) for s in (:x,:y,:z)],expression;semantic_id="validator/rejection")
        index=length(original.nodes);arity=length(original.nodes[index].inputs)
        for column in 1:arity+1, value in (-1,0,3)
            changed=deepcopy(original);changed.nodes[index].rows[1][column]=value
            @test_throws ArgumentError validate_atomic_plan(changed)
        end
        for delta in (-1,1)
            changed=deepcopy(original)
            delta==-1 ? pop!(changed.nodes[index].rows[1]) : push!(changed.nodes[index].rows[1],1)
            @test_throws ArgumentError validate_atomic_plan(changed)
        end
        changed=deepcopy(original);changed.nodes[index].rows[2]=copy(changed.nodes[index].rows[1])
        @test_throws ArgumentError validate_atomic_plan(changed)
        changed=deepcopy(original);pop!(changed.nodes[index].rows)
        @test_throws ArgumentError validate_atomic_plan(changed)
        changed=deepcopy(original);changed.nodes[index].rows[1][end]=3-changed.nodes[index].rows[1][end]
        @test_throws ArgumentError validate_atomic_plan(changed)
        @test_throws ArgumentError compile_atomic(changed)
        for input in (0,index,index+1)
            changed=deepcopy(original);changed.nodes[index].inputs[1]=input
            @test_throws ArgumentError validate_atomic_plan(changed)
        end
        changed=deepcopy(original);push!(changed.nodes[index].domain,changed.nodes[index].domain[1])
        @test_throws ArgumentError validate_atomic_plan(changed)
        changed=deepcopy(original);empty!(changed.nodes[index].domain)
        @test_throws ArgumentError validate_atomic_plan(changed)
        @test_throws CompilationLimit validate_atomic_plan(original;max_rows=length(original.nodes[index].rows)-1)
        @test validate_atomic_plan(original;max_rows=length(original.nodes[index].rows))
    end
    p=atomic_plan([structured_codebook(s,0:1) for s in (:x,:y)],N(:eq,:x,:y);semantic_id="validator/leaves")
    replace_node(plan,i;operator=plan.nodes[i].operator,inputs=plan.nodes[i].inputs,
        data=plan.nodes[i].data,domain=plan.nodes[i].domain,rows=plan.nodes[i].rows)=
        (plan.nodes[i]=AtomicNode(operator,inputs,data,domain,rows);plan)
    for data in (0,3,"x")
        @test_throws ArgumentError validate_atomic_plan(replace_node(deepcopy(p),1;data))
    end
    @test_throws ArgumentError validate_atomic_plan(replace_node(deepcopy(p),1;inputs=[1]))
    @test_throws ArgumentError validate_atomic_plan(replace_node(deepcopy(p),1;rows=[[1]]))
    @test_throws ArgumentError validate_atomic_plan(replace_node(deepcopy(p),1;domain=BigInt[0,2]))
    @test_throws ArgumentError validate_atomic_plan(replace_node(deepcopy(p),2;data=1))
    @test_throws ArgumentError validate_atomic_plan(replace_node(deepcopy(p),3;operator=:unknown))
    @test_throws ArgumentError validate_atomic_plan(replace_node(deepcopy(p),3;inputs=[1]))
    @test_throws ArgumentError validate_atomic_plan(AtomicPlan(p.codebooks,p.nodes,0,p.semantic_id))
    @test_throws ArgumentError validate_atomic_plan(AtomicPlan(p.codebooks,p.nodes,4,p.semantic_id))
    changed=deepcopy(p);changed.codebooks[2]=deepcopy(changed.codebooks[1])
    @test_throws ArgumentError validate_atomic_plan(changed)
    @test_throws ArgumentError validate_atomic_plan(atomic_plan([structured_codebook(:x,[-2,3])],N(:eq,:x,0);semantic_id="validator/nonbool") |> q->AtomicPlan(q.codebooks,q.nodes,1,q.semantic_id))
    budget=atomic_plan([structured_codebook(:x,[2,3])],N(:eq,N(:sqr,:x),4);semantic_id="validator/value-budget")
    @test_throws CompilationLimit validate_atomic_plan(budget;max_value_bits=2)
    @test validate_atomic_plan(budget;max_value_bits=4)
end

@testitem "Atomic validation tuples preserve large integer and partial-operation checks" tags=[:q10] begin
    N=IntensionNode
    large=big(2)^80+1
    for op in (:add,:sub,:mul,:dist,:min,:max)
        books=[structured_codebook(:x,[-large,large]),structured_codebook(:y,[big(2)^75,-big(2)^75])]
        p=atomic_plan(books,N(:eq,N(op,:x,:y),0);semantic_id="validator/large/$(op)",max_value_bits=512)
        @test validate_atomic_plan(p;max_value_bits=512)
        node=p.nodes[3]
        oracle=op===:add ? (+) : op===:sub ? (-) : op===:mul ? (*) :
            op===:dist ? ((x,y)->abs(x-y)) : op===:min ? min : max
        for row in node.rows
            @test node.domain[last(row)]==oracle(books[1].values[row[1]],books[2].values[row[2]])
        end
    end
    for op in (:pow,:pow_total), sign in (-1,1)
        exponents=sign==1 || op===:pow_total ? sign .* [big(2)^80,big(2)^80+1] : BigInt[0,1]
        books=[structured_codebook(:x,[-1,1]),structured_codebook(:y,exponents)]
        p=atomic_plan(books,N(:eq,N(op,:x,:y),1);semantic_id="validator/large-exponent/$(op)/$(sign)")
        @test validate_atomic_plan(p)
        for row in p.nodes[3].rows
            @test p.nodes[3].domain[last(row)]==(books[1].values[row[1]]==-1 && isodd(exponents[row[2]]) ? -1 : 1)
        end
    end
    # Constructed plans are untrusted: checked table sizes do not permit an
    # undefined primitive operation or a value-growth exception to be bypassed.
    p=atomic_plan([structured_codebook(:x,[0,1]),structured_codebook(:y,[0,1])],N(:eq,:x,:y);semantic_id="validator/partial")
    for op in (:div,:mod)
        changed=deepcopy(p);node=changed.nodes[3]
        changed.nodes[3]=AtomicNode(op,node.inputs,node.data,node.domain,node.rows)
        @test_throws ArgumentError validate_atomic_plan(changed)
    end
    p=atomic_plan([structured_codebook(:x,[2,3]),structured_codebook(:y,[-1,0])],N(:eq,:x,:y);semantic_id="validator/negative-exponent")
    node=p.nodes[3];p.nodes[3]=AtomicNode(:pow,node.inputs,node.data,node.domain,node.rows)
    @test_throws ArgumentError validate_atomic_plan(p)
end
