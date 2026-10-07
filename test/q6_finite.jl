@testitem "Q6 explicit finite lex, element, precedence and historical stretch" tags=[:q6] begin
    for mode in (:native,:domain_wall,:zero_one_hot,:binary,:gray,:fibonacci)
        books=[structured_codebook(Symbol(:x,i),[0,1];encoding=mode) for i in 1:3]
        for op in (:lt,:le,:gt,:ge)
            q=lex_component(books,[:x1,:x2],[:x2,:x3];operator=op)
            oracle=v->begin
                a=(v[1],v[2]); b=(v[2],v[3])
                op===:lt ? a<b : op===:le ? a<=b : op===:gt ? a>b : a>=b
            end
            @test exhaustive_check(q,oracle;oracle_id="q6/lex-independent",profile=:indicator_exact).status===:pass
        end
        q=element_component(books,[:x1,:x2],IntegerCondition(:eq,1);index=:x3)
        @test exhaustive_check(q,v->v[v[3]+1]==1;oracle_id="q6/element",profile=:indicator_exact).status===:pass
        q=element_component(books,[:x1,:x2],IntegerCondition(:eq,1);index=:x1)
        @test exhaustive_check(q,v->v[v[1]+1]==1;oracle_id="q6/element-alias",profile=:indicator_exact).status===:pass
        q=element_component(books,[:x1,:x2],IntegerCondition(:eq,1))
        @test exhaustive_check(q,v->1 in v[1:2];oracle_id="q6/element-exists",profile=:indicator_exact).status===:pass
        for covered in (false,true)
            q=precedence_component(books,[:x1,:x2,:x3],[0,1];covered)
            @test exhaustive_check(q,v->v[1]==0 && (!covered || 1 in v);oracle_id="q6/precedence",profile=:indicator_exact).status===:pass
        end
        q=stretch_component(books,[:x1,:x2,:x3],[0,1],[1:1,1:1])
        @test exhaustive_check(q,v->v[1]!=v[2] && v[2]!=v[3];oracle_id="q6/stretch",profile=:indicator_exact).status===:pass
        q=stretch_component(books,[:x1,:x2,:x3],[0,1],[1:3,1:3];patterns=[(0,1)])
        @test exhaustive_check(q,v->issorted(v);oracle_id="q6/stretch-patterns",profile=:indicator_exact).status===:pass
    end
    # The fixed invalid-code extension of the generic route is expensive for one-hot.
    # Certify smaller one-hot instances; retain an explicit UNKNOWN regression at n=3.
    oh=[structured_codebook(Symbol(:x,i),[0,1];encoding=:one_hot) for i in 1:2]
    for (q,p) in ((lex_component(oh,[:x1],[:x2];operator=:lt),v->v[1]<v[2]),
            (element_component(oh,[:x1],IntegerCondition(:eq,1);index=:x2),v->v[2]==0 && v[1]==1),
            (precedence_component(oh,[:x1,:x2],[0,1]),v->v[1]==0),
            (stretch_component(oh,[:x1,:x2],[0,1],[1:1,1:1]),v->v[1]!=v[2]))
        @test exhaustive_check(q,p;oracle_id="q6/small-onehot",profile=:indicator_exact).status===:pass
    end
    oh3=[oh;[structured_codebook(:x3,[0,1];encoding=:one_hot)]]
    large=lex_component(oh3,[:x1,:x2],[:x2,:x3];operator=:le)
    @test exhaustive_check(large,v->(v[1],v[2])<=(v[2],v[3]);oracle_id="q6/budget-boundary").status===:unknown
    books=[structured_codebook(:i,[-1,0,1,2];encoding=:domain_wall),structured_codebook(:x,[1];encoding=:domain_wall)]
    q=element_component(books,[:x],IntegerCondition(:eq,1);index=:i,start_index=1)
    @test exhaustive_check(q,v->v[1]==1;oracle_id="q6/index-bounds").status===:pass
    @test_throws CompilationLimit lex_component(books,[:i],[:x];max_primary_states=1)
    @test_throws ArgumentError precedence_component(books,[:x],[1,1])
    @test_throws ArgumentError stretch_component(books,[:x],[1],[0:1])
    @test_throws DimensionMismatch lex_component(books,[:x],[:x,:x])
end
