@testitem "Q7 cumulative whole timeline and finite variants" tags=[:q7] begin
    for mode in (:domain_wall,:zero_one_hot,:one_hot)
        books=[structured_codebook(:a,[-2,0];encoding=mode),structured_codebook(:b,[-2,0];encoding=mode)]
        for lengths in ([1,1],[0,1],[2,3]), heights in ([1,1],[0,2]), op in (:le,:eq,:ne,:ge,:lt,:gt)
            condition=IntegerCondition(op,1)
            oracle=v->all(t->satisfies(condition,sum(heights[i] for i in 1:2 if v[i]<=t<v[i]+lengths[i];init=0)),-3:4)
            for route in (:events,:finite)
                q=cumulative_component(books,[:a,:b],lengths,heights,condition;route)
                @test exhaustive_check(q,oracle;oracle_id="q7/cumulative-independent-grid").status===:pass
            end
        end
        # Empty exterior must be checked even for positive lower-bound conditions.
        for condition in (IntegerCondition(:eq,0),IntegerCondition(:in,[0,2]),IntegerCondition(:notin,1:2),IntegerCondition(:in,1:0))
            q=cumulative_component(books,[:a,:a,1],[1,2,1],[1,1,1],condition)
            oracle=v->all(t->satisfies(condition,Int(v[1]<=t<v[1]+1)+Int(v[1]<=t<v[1]+2)+Int(t==1)),-3:3)
            @test exhaustive_check(q,oracle;oracle_id="q7/cumulative-alias-and-constant").status===:pass
        end
    end
    books=[structured_codebook(:o,[-1,0];encoding=:domain_wall),structured_codebook(:l,[0,1];encoding=:native),structured_codebook(:h,[0,1];encoding=:native)]
    q=cumulative_component(books,[:o,0],[:l,1],[:h,1],IntegerCondition(:le,1))
    @test exhaustive_check(q,v->!(v[1]==0 && v[2]==1 && v[3]==1);oracle_id="q7/variable-tasks").status===:pass
    q=cumulative_component(books,[:o],[:l],[:h],IntegerCondition(:le,1);ends=[0])
    @test exhaustive_check(q,v->v[1]+v[2]==0;oracle_id="q7/explicit-ends").status===:pass
    @test_throws ArgumentError cumulative_component(books,[:o],[-1],[1],IntegerCondition(:le,1))
    @test_throws ArgumentError cumulative_component(books,[:o],[:o],[1],IntegerCondition(:le,1))
    @test_throws ArgumentError cumulative_component(books,[:o],[1],[1],IntegerCondition(:le,1);time_domain=:naturals)
    @test_throws ArgumentError cumulative_component(books,[:o],[:l],[1],IntegerCondition(:le,1);route=:events)
    @test_throws CompilationLimit cumulative_component(books,[:o],[1],[1],IntegerCondition(:le,1);max_events=0)
    huge=big(10)^30
    b=[structured_codebook(:a,[-huge,huge];encoding=:domain_wall)]
    q=cumulative_component(b,[:a,0],[1,1],[1,1],IntegerCondition(:le,1);max_events=6)
    @test exhaustive_check(q,v->true;oracle_id="q7/no-horizon-scan").status===:pass
    q=cumulative_component(b,Int[],Int[],Int[],IntegerCondition(:gt,0))
    @test exhaustive_check(q,v->false;oracle_id="q7/empty-exterior").status===:pass
end

@testitem "Q7 pairwise noOverlap and zero extents" tags=[:q7] begin
    for mode in (:native,:domain_wall,:zero_one_hot,:one_hot)
        books=[structured_codebook(:a,[0,1];encoding=mode),structured_codebook(:b,[0,1];encoding=mode)]
        for lengths in ([1,1],[0,1],[2,1]), ignored in (false,true)
            q=nooverlap_component(books,[:a,:b],lengths;zero_ignored=ignored)
            p=v->(ignored && 0 in lengths) || v[1]+lengths[1]<=v[2] || v[2]+lengths[2]<=v[1]
            @test exhaustive_check(q,p;oracle_id="q7/noOverlap-1d").status===:pass
        end
        for ignored in (false,true)
            q=nooverlap_component(books,[:a 1;0 :b],[1 0;1 2];zero_ignored=ignored)
            @test exhaustive_check(q,v->ignored || v[1]==1 || v[2]==1;oracle_id="q7/noOverlap-kd-zero").status===:pass
            q=nooverlap_component(books,[0,0],[:a,:b];zero_ignored=ignored)
            @test exhaustive_check(q,v->v[1]==0 || v[2]==0;oracle_id="q7/noOverlap-variable-length").status===:pass
        end
        q=nooverlap_component(books,[:a,:b,:a],[1,1,0])
        @test exhaustive_check(q,v->v[1]!=v[2];oracle_id="q7/noOverlap-shared-pairs").status===:pass
        @test_throws CompilationLimit nooverlap_component(books,[:a,:b],[1,1];max_pairs=0)
    end
end

@testitem "Q7 circuit requires one nontrivial cycle" tags=[:q7] begin
    for values in ([1,0,3,2],[1,0,2,3],[0,1,2,3],[1,2,3,0]), start in (0,5)
        books=[structured_codebook(Symbol(:s,i),[x+start];encoding=:domain_wall) for (i,x) in enumerate(values)]
        expected=values in ([1,0,2,3],[1,2,3,0])
        q=circuit_component(books,[b.variable for b in books];start_index=start)
        @test exhaustive_check(q,v->expected;oracle_id="q7/circuit-permutation-not-enough").status===:pass
    end
    books=[structured_codebook(:a,[0,1];encoding=:native),structured_codebook(:b,[0,1];encoding=:native)]
    for size in (nothing,1,2,3,:a)
        q=circuit_component(books,[:a,:b];size)
        @test exhaustive_check(q,v->v==[1,0] && size in (nothing,2);oracle_id="q7/circuit-size").status===:pass
    end
    books=[structured_codebook(:a,[1,2];encoding=:domain_wall),structured_codebook(:b,[0,2];encoding=:domain_wall),structured_codebook(:c,[0,1];encoding=:domain_wall)]
    q=circuit_component(books,[:a,:b,:c])
    @test exhaustive_check(q,v->v in ([1,2,0],[2,0,1]);oracle_id="q7/circuit-3cycle").status===:pass
end

@testitem "Q7 finite intension DAG operations and budgets" tags=[:q7] begin
    for mode in (:domain_wall,:zero_one_hot,:one_hot)
        books=[structured_codebook(:a,[-1,1];encoding=mode),structured_codebook(:b,[-1,1];encoding=mode)]
        for (op,f) in ((:add,+),(:sub,-),(:mul,*),(:dist,(x,y)->abs(x-y)),(:min,min),(:max,max))
            e=IntensionNode(:eq,IntensionNode(op,:a,:b),0)
            @test exhaustive_check(intension_component(books,e),v->f(v...)==0;oracle_id="q7/arithmetic-dag").status===:pass
        end
        for (op,f) in ((:neg,-),(:abs,abs),(:sqr,x->x^2))
            e=IntensionNode(:eq,IntensionNode(op,:a),1)
            @test exhaustive_check(intension_component(books,e),v->f(v[1])==1;oracle_id="q7/unary-dag").status===:pass
        end
        for (op,f) in ((:eq,==),(:ne,!=),(:lt,<),(:le,<=),(:gt,>),(:ge,>=))
            @test exhaustive_check(intension_component(books,IntensionNode(op,:a,:b)),v->f(v...);oracle_id="q7/comparison-dag").status===:pass
        end
        shared=IntensionNode(:add,:a,:b)
        e=IntensionNode(:eq,shared,shared)
        @test exhaustive_check(intension_component(books,e;max_nodes=4),v->true;oracle_id="q7/shared-dag").status===:pass
    end
    b=[structured_codebook(:a,[0,1];encoding=:native),structured_codebook(:b,[0,1];encoding=:native)]
    for (op,f) in ((:and,(x,y)->x==1 && y==1),(:or,(x,y)->x==1 || y==1),(:xor,(x,y)->x!=y),(:iff,==),(:imp,(x,y)->x==0 || y==1))
        @test exhaustive_check(intension_component(b,IntensionNode(op,:a,:b)),v->f(v...);oracle_id="q7/logical-dag").status===:pass
    end
    @test exhaustive_check(intension_component(b,IntensionNode(:not,:a)),v->v[1]==0;oracle_id="q7/not").status===:pass
    @test exhaustive_check(intension_component(b,IntensionNode(Symbol("if"),:a,:b,1)),v->v[1]==0 || v[2]==1;oracle_id="q7/if").status===:pass
    @test_throws ArgumentError intension_component(b,IntensionNode(:div,:a,:b))
    @test_throws ArgumentError intension_component(b,IntensionNode(:add,:a))
    @test_throws ArgumentError intension_component(b,IntensionNode(:not,2))
    @test_throws ArgumentError intension_component(b,IntensionNode(:add,:a,2))
    @test_throws CompilationLimit intension_component(b,IntensionNode(:eq,:a,:b);max_nodes=1)
    @test_throws CompilationLimit intension_component(b,IntensionNode(:eq,:a,:b);max_depth=1)
    @test_throws CompilationLimit intension_component(b,IntensionNode(:eq,:a,:b);max_edges=1)
    repeated=foldl((e,_)->IntensionNode(:add,e,e),1:20;init=:a)
    q=intension_component(b,IntensionNode(:eq,repeated,0);max_nodes=23,max_edges=42)
    @test exhaustive_check(q,v->v[1]==0;oracle_id="q7/shared-dag-no-exponential-traversal").status===:pass
    @test_throws CompilationLimit intension_component(b,IntensionNode(:eq,IntensionNode(:mul,7,7),1);max_value_bits=3)
end
