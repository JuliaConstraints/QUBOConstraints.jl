@testitem "Q6 nary supports, conflicts, wildcards and repeated identities" tags=[:q6] begin
    for mode in (:native,:one_hot,:zero_one_hot,:domain_wall), supports in (true,false)
        books=[structured_codebook(Symbol(:x,i),[0,1];encoding=mode) for i in 1:3]
        for rows in ([],[[0,1,0]],[[ANY_VALUE,1,0],[1,ANY_VALUE,1]],[[ANY_VALUE,ANY_VALUE,ANY_VALUE]])
            q=table_component(books,[:x1,:x2,:x3],rows;supports)
            oracle=v->any(r->all(r[i] isa AnyValue || r[i]==v[i] for i in 1:3),rows)==supports
            @test exhaustive_check(q,oracle;oracle_id="q6/table").status===:pass
        end
        rows=[[0,1,0],[1,1,ANY_VALUE]]
        q=table_component(books,[:x1,:x1,:x2],rows;supports)
        @test exhaustive_check(q,v->(v[1]==1)==supports;oracle_id="q6/repeated").status===:pass
        q=table_component(books,Symbol[],[[]];supports)
        @test exhaustive_check(q,v->supports;oracle_id="q6/nullary").status===:pass
    end
    for mode in (:one_hot,:zero_one_hot,:domain_wall)
        books=[structured_codebook(:x,[0,1];encoding=:native)]
        q=table_component(books,[:x],[[0],[1]];selector_encoding=mode)
        @test exhaustive_check(q,v->true;oracle_id="q6/selector").status===:pass
    end
end

@testitem "Q6 automaton and layered paths: independent language oracle" tags=[:q6] begin
    function accepts(word,edges,start,finals)
        states=Set([start])
        for (i,value) in enumerate(word)
            states=Set(e[3] for e in edges[i] if e[1] in states && e[2]==value)
        end
        return any(in(finals),states)
    end
    for mode in (:native,:one_hot,:zero_one_hot,:domain_wall), edge_mode in (:one_hot,:domain_wall)
        books=[structured_codebook(Symbol(:x,i),[0,1];encoding=mode) for i in 1:2]
        transitions=[(0,0,0),(0,1,1),(1,0,1),(1,1,0)]
        q=regular_component(books,[:x1,:x2],transitions;start=0,finals=[0],edge_encoding=edge_mode)
        @test exhaustive_check(q,v->iseven(sum(v));oracle_id="q6/parity").status===:pass
        q=regular_component(books,[:x1,:x1],transitions;start=0,finals=[0],edge_encoding=edge_mode)
        @test exhaustive_check(q,v->true;oracle_id="q6/repeated-path").status===:pass
        layers=[[(0,0,1),(0,0,2),(0,1,3)],[(1,1,9),(2,0,9),(3,1,9),(99,0,9)]]
        q=mdd_component(books,[:x1,:x2],layers;start=0,finals=[9],edge_encoding=edge_mode)
        @test exhaustive_check(q,v->accepts(v,layers,0,[9]);oracle_id="q6/nondeterminism").status===:pass
        q=regular_component(books,[:x1,:x2],transitions;start=0,finals=[99])
        @test exhaustive_check(q,v->false;oracle_id="q6/no-path").status===:pass
    end
    for finals in ([0],[1])
        q=mdd_component(StructuredCodebook[],Symbol[],[];start=0,finals)
        @test exhaustive_check(q,v->0 in finals;oracle_id="q6/empty-word").status===:pass
    end
end

@testitem "Q6 slide scopes and independent existential witnesses" tags=[:q6] begin
    for n in 0:6, width in 1:7, step in 1:8, circular in (true,false)
        scope=collect(1:n)
        q=slide_scopes(scope,width;offset=step,circular)
        starts=filter(i->circular || i+width<=n,collect(0:step:n-1))
        expected=[[mod(i+j,n)+1 for j in 0:width-1] for i in starts]
        @test q==expected
    end
    books=[structured_codebook(Symbol(:x,i),[0,1];encoding=:native) for i in 1:4]
    builder=(bs,scope,ns)->table_component(bs,scope,[[ANY_VALUE,ANY_VALUE,0]];namespace=ns)
    q=slide_component(books,[:x1,:x2,:x3,:x4],3,builder;circular=true)
    @test exhaustive_check(q,v->all(iszero,v);oracle_id="q6/slide-circular").status===:pass
    q=slide_component(books,[:x1,:x2,:x3,:x4],3,builder;circular=false)
    @test exhaustive_check(q,v->v[3]==v[4]==0;oracle_id="q6/slide-linear").status===:pass
    q=slide_component(books,[:x1,:x2,:x3,:x4],3,builder;circular=true,offset=3)
    @test exhaustive_check(q,v->v[3]==v[2]==0;oracle_id="q6/slide-offset").status===:pass
    # Different windows must not be forced to choose the same support tuple.
    equal=(bs,scope,ns)->table_component(bs,scope,[[0,0],[1,1]];namespace=ns)
    q=slide_component(books,[:x1,:x2,:x3,:x4],2,equal;offset=2)
    @test exhaustive_check(q,v->v[1]==v[2] && v[3]==v[4];oracle_id="q6/private").status===:pass
end

@testitem "Q6 structural input and constructive budget errors" tags=[:q6] begin
    b=[structured_codebook(:x,[0,1];encoding=:native)]
    @test_throws ArgumentError table_component(b,[:unknown],[[0]])
    @test_throws DimensionMismatch table_component(b,[:x],[[0,1]])
    @test_throws CompilationLimit table_component(b,[:x],[[0],[1]];max_rows=1)
    @test_throws CompilationLimit table_component(b,[:x],[[0]];max_auxiliaries=0)
    @test_throws ArgumentError table_component(b,[:x],[[0]];namespace=:x)
    @test_throws CompilationLimit mdd_component(b,[:x],[[(0,0,1),(0,1,1)]];start=0,finals=[1],max_edges=1)
    @test_throws DimensionMismatch mdd_component(b,[:x],[];start=0,finals=[1])
    @test_throws ArgumentError slide_scopes([:x],0)
    @test_throws ArgumentError slide_scopes([:x],1;offset=0)
    @test_throws CompilationLimit slide_scopes([:x],big(10)^10;circular=true)
    @test_throws CompilationLimit slide_scopes([:x,:x],1;max_windows=1)
end
