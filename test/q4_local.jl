@testitem "Q4 local relations have global validity" tags=[:q4] begin
    cases = [
        (:instantiation,[[-2,3],[1,5]],(;values=[3,1]),v->v==[3,1]),
        (:instantiation,[[1,2]],(;values=[9]),v->false),
        (:all_equal,[[1,3],[3,7],[1,3]],(;),v->all(==(v[1]),v)),
        (:all_different,[[0,1],[0,1],[0,2]],(;except=[0]),v->all(v[i]!=v[j] || v[i]==0 for i in 1:3 for j in i+1:3)),
        (:ordered,[[-3,2],[0,5],[1,9]],(;operator=:lt,distances=[2,-1]),v->v[1]+2<v[2] && v[2]-1<v[3]),
        (:ordered,[[2,5],[1,7]],(;operator=:ge,distances=[3]),v->v[1]+3>=v[2]),
        (:channel,[[-1,0,1],[0,1,2]],(;),v->all(x->x in 0:1,v) && all(v[v[i]+1]==i-1 for i in 1:2)),
        (:channel,[[4,5],[4,5]],(;start_index=4),v->all(v[v[i]-3]==i+3 for i in 1:2)),
        (:no_overlap,[[0,2],[1,4]],(;lengths=[0,2]),v->true),
        (:no_overlap,[[0,2],[1,4]],(;lengths=[0,2],zero_ignored=false),v->v[1]<=v[2] || v[2]+2<=v[1]),
        (:no_overlap,[[0,2],[1,4]],(;lengths=[2,2]),v->v[1]+2<=v[2] || v[2]+2<=v[1])]
    for mode in (:one_hot,:domain_wall,:zero_one_hot), (family,domains,kwargs,predicate) in cases
        q = local_constraint(family,domains;encoding=mode,kwargs...)
        @test all(b->b.role===:primary,q.bits)
        @test exhaustive_check(q,predicate;oracle_id="q4/$(family)/$(mode)").status===:pass
    end
    for mode in (:one_hot,:domain_wall,:zero_one_hot,:native)
        q = local_constraint(:all_different,[[false,true],[false,true]];encoding=mode)
        @test exhaustive_check(q,v->v[1]!=v[2];oracle_id="q4/native").status===:pass
    end
    for mode in (:one_hot,:domain_wall,:zero_one_hot)
        q = local_constraint(:ordered,[[typemax(Int)],[typemax(Int)]];encoding=mode,distances=[1],operator=:gt)
        @test exhaustive_check(q,v->true;oracle_id="q4/overflow").status===:pass
        q = local_constraint(:all_equal,[];encoding=mode)
        @test exhaustive_check(q,v->true;oracle_id="q4/empty").status===:pass
    end
end

@testitem "Q4 tables and scope checks" tags=[:q4] begin
    for mode in (:one_hot,:domain_wall,:zero_one_hot), supports in (false,true)
        q = table_constraint([[1,3],[2,7]],[(1,ANY_VALUE),(3,7)];encoding=mode,supports)
        @test exhaustive_check(q,v->((v[1]==1 || v[2]==7)==supports);oracle_id="q4/table").status===:pass
        q = table_constraint([[1,3,5]],[(1,),(5,)];encoding=mode,supports)
        @test exhaustive_check(q,v->((v[1] in (1,5))==supports);oracle_id="q4/unary-table").status===:pass
        q = table_constraint([[1,2]],Tuple{Int}[];encoding=mode,supports)
        @test exhaustive_check(q,v->!supports;oracle_id="q4/empty-table").status===:pass
        book = structured_codebook(:x,0:2;encoding=mode)
        q = local_component([book];binary=[BinaryRelation(:x,:x,[i!=j for i in 0:2,j in 0:2])])
        @test exhaustive_check(q,v->false;oracle_id="q4/repeated-scope").status===:pass
    end
    q = local_constraint(:all_equal,[[1,3],[1,3]];encoding=[:one_hot,:domain_wall])
    @test exhaustive_check(q,v->v[1]==v[2];oracle_id="q4/mixed").status===:pass
    q = table_constraint([["*","a"]],[("*",)])
    @test exhaustive_check(q,v->v[1]=="*";oracle_id="q4/literal-star").status===:pass
    @test_throws ArgumentError local_constraint(:unknown,[[1]])
    @test_throws ArgumentError local_constraint(:all_equal,[[1],[2]];except=[0])
    @test_throws ArgumentError local_constraint(:ordered,[[1],[2]];operator=:eq)
    @test_throws ArgumentError local_constraint(:no_overlap,[[1],[2]];lengths=[-1,1])
    @test_throws DimensionMismatch table_constraint([[1],[2]],[(1,)])
    @test_throws ArgumentError table_constraint([[1],[2],[3]],[(1,2,3)])
    @test_throws ArgumentError local_component([structured_codebook(:x,0:4;encoding=:binary)])
    @test_throws DimensionMismatch local_component([structured_codebook(:x,0:2)];unary=[UnaryRelation(:x,[true])])
    @test_throws CompilationLimit local_constraint(:all_equal,[0:1000000,0:1000000];max_table_entries=100)
    @test_throws CompilationLimit table_constraint([0:1000000,0:1000000],[];max_table_entries=100)
    @test_throws CompilationLimit local_component([structured_codebook(:x,0:100)];max_terms=10)
end

@testitem "Q4 validity guard preserves independent auxiliary minima" tags=[:q4] begin
    for mode in (:one_hot,:domain_wall,:binary,:gray,:fibonacci,:unary)
        book = structured_codebook(:x,0:3;encoding=mode)
        y = BitID(:a,1;role=:semantic_auxiliary)
        # Existing auxiliary minimizes to zero independently of the primary code.
        raw = QUBOComponent([book.bits;y];linear=[(length(book.bits)+1,3)],codebooks=[book])
        guarded = guard_encodings(raw)
        @test exhaustive_check(guarded.component,v->true;oracle_id="q4/independent-guard").status===:pass
    end
    for mode in (:one_hot,:domain_wall)
        book = codebook(:x,0:3;encoding=mode)
        @test exhaustive_check(encoding_validity(book),v->true;oracle_id="q4/explicit-validity").status===:pass
    end
end
