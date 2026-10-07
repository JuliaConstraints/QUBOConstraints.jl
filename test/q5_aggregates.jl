@testitem "Q5 aggregate families: exhaustive witnesses and invalid codes" tags=[:q5] begin
    for mode in (:one_hot,:domain_wall,:zero_one_hot), n in (1,2,3)
        books = [structured_codebook(Symbol(:x,i),[-1,2];encoding=mode) for i in 1:n]
        for c in (IntegerCondition(:eq,1),IntegerCondition(:le,2),IntegerCondition(:ne,1))
            q = nvalues_component(books,c;slack_encoding=:domain_wall)
            @test exhaustive_check(q,v->satisfies(c,length(unique(v)));oracle_id="q5/nvalues").status===:pass
            q = nvalues_component(books,c;except=[-1],slack_encoding=:domain_wall)
            @test exhaustive_check(q,v->satisfies(c,length(unique(filter(!=(-1),v))));oracle_id="q5/nvalues-except").status===:pass
        end
        for kind in (:minimum,:maximum), op in (:eq,:ne,:lt,:le,:gt,:ge,:in,:notin)
            c = IntegerCondition(op,op in (:in,:notin) ? [-1,0] : 0)
            q = extremum_component(books,kind,c;slack_encoding=:domain_wall)
            f = kind===:minimum ? minimum : maximum
            @test exhaustive_check(q,v->satisfies(c,f(v));oracle_id="q5/extremum").status===:pass
        end
        for closed in (false,true)
            cs = [IntegerCondition(:le,1)]
            q = cardinality_component(books,[2],cs;closed,slack_encoding=:domain_wall)
            @test exhaustive_check(q,v->count(==(2),v)<=1 && (!closed || all(==(2),v));oracle_id="q5/cardinality").status===:pass
        end
        q = knapsack_component(books,collect(1:n),fill(-2,n),IntegerCondition(:le,2),IntegerCondition(:ge,-3);
            slack_encoding=:bounded_binary)
        @test exhaustive_check(q,v->sum(i*big(v[i]) for i in eachindex(v))<=2 && -2sum(v)>=-3;
            oracle_id="q5/knapsack").status===:pass
    end
end

@testitem "Q5 bin packing: empty bins, zero size and separated semantic profiles" tags=[:q5] begin
    for mode in (:one_hot,:domain_wall,:zero_one_hot), sizes in ([1,1],[0,1],[0,0])
        books = [structured_codebook(Symbol(:x,i),[0,1];encoding=mode) for i in 1:2]
        for op in (:eq,:ne,:le,:gt)
            c = IntegerCondition(op,1)
            q = binpacking_component(books,sizes,c)
            predicate = v->all(satisfies(c,sum((sizes[i] for i in eachindex(v) if v[i]==label);init=0)) for label in unique(v))
            @test exhaustive_check(q,predicate;oracle_id="q5/bin-common").status===:pass
            cs = [IntegerCondition(:eq,2),IntegerCondition(:eq,1)]
            for profile in (:validator,:pdf_occupied)
                q = binpacking_component(books,sizes,cs;bins=[0,1],profile)
                predicate = v->all((profile===:pdf_occupied && !(label in v)) ||
                    satisfies(cs[label+1],sum((sizes[i] for i in eachindex(v) if v[i]==label);init=0)) for label in 0:1)
                @test exhaustive_check(q,predicate;oracle_id="q5/bin-$(profile)").status===:pass
            end
        end
    end
    books = [structured_codebook(:x,[0,1];encoding=:native)]
    @test_throws ArgumentError binpacking_component(books,[-1],IntegerCondition(:eq,0))
    @test_throws ArgumentError binpacking_component(books,[1],[IntegerCondition(:eq,0)];bins=[0])
    @test_throws ArgumentError binpacking_component(books,[1],[IntegerCondition(:eq,0)])
    @test_throws CompilationLimit binpacking_component(books,[10],IntegerCondition(:le,10);max_auxiliaries=2)
    q = binpacking_component(StructuredCodebook[],Int[],[IntegerCondition(:eq,1)];bins=[0])
    @test exhaustive_check(q,v->false;oracle_id="q5/bin-empty-declared").status===:pass
    q = binpacking_component(StructuredCodebook[],Int[],IntegerCondition(:eq,1))
    @test exhaustive_check(q,v->true;oracle_id="q5/bin-empty-common").status===:pass
end

@testitem "Q5 witness hygiene and explicit finite variants" tags=[:q5] begin
    books = [structured_codebook(Symbol(:x,i),[0,1];encoding=:native) for i in 1:3]
    # x1 is both a summand and a variable coefficient: one primary identity, repeated use.
    predicate = v->v[1]*v[2]+v[1]*v[3]==v[2]
    q = finite_arithmetic_component(books,predicate;semantic_id="q5/sum/variable-coefficients/repeated")
    @test exhaustive_check(q,predicate;oracle_id="q5/finite",profile=:indicator_exact).status===:pass
    a = nvalues_component(books,IntegerCondition(:eq,1))
    b = nvalues_component(books,IntegerCondition(:eq,2))
    combined = compose(a,b)
    @test exhaustive_check(combined,v->false;oracle_id="q5/private-contradiction").status===:pass
    @test all(bit.role!==:primary || startswith(string(bit.owner),"x") for bit in a.bits)
    @test length(a.codebooks)==3
    @test_throws CompilationLimit nvalues_component(books,IntegerCondition(:eq,1);max_auxiliaries=1)
    @test_throws ArgumentError nvalues_component([structured_codebook(Symbol("nvalues/x"),[0,1])],IntegerCondition(:eq,1))
    @test_throws ArgumentError cardinality_component(books,[0,0],[IntegerCondition(:eq,0),IntegerCondition(:eq,1)])
    @test_throws ArgumentError extremum_component(StructuredCodebook[],:minimum,IntegerCondition(:eq,0))
    q = nvalues_component(StructuredCodebook[],IntegerCondition(:eq,0))
    @test exhaustive_check(q,v->true;oracle_id="q5/nvalues-empty").status===:pass
    singleton = [structured_codebook(:x,[big(typemax(Int))+1];encoding=:domain_wall)]
    q = extremum_component(singleton,:maximum,IntegerCondition(:eq,big(typemax(Int))+1))
    @test exhaustive_check(q,v->true;oracle_id="q5/extremum-bigint").status===:pass
    for mode in (:one_hot,:domain_wall,:zero_one_hot)
        nonuniform = [structured_codebook(Symbol(:x,i),[-3,0,5];encoding=mode) for i in 1:2]
        for kind in (:minimum,:maximum), c in (IntegerCondition(:in,[-3,5]),IntegerCondition(:notin,-1:4))
            q = extremum_component(nonuniform,kind,c;slack_encoding=:domain_wall)
            f = kind===:minimum ? minimum : maximum
            @test exhaustive_check(q,v->satisfies(c,f(v));oracle_id="q5/nonuniform").status===:pass
        end
        q = nvalues_component(nonuniform,IntegerCondition(:ne,1);slack_encoding=:domain_wall)
        @test exhaustive_check(q,v->length(unique(v))!=1;oracle_id="q5/nonuniform-count").status===:pass
    end
    @test_throws CompilationLimit binpacking_component(books,ones(Int,3),IntegerCondition(:le,2);max_auxiliaries=7)
end
