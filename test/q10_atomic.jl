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
