@testitem "Q4 all three-bit fixed truth tables" tags=[:q4] begin
    books = [structured_codebook(Symbol(:x,i),0:1;encoding=:native) for i in 1:3]
    for f in 0:255
        predicate = v->iszero(f & (1<<sum(v[i]<<(i-1) for i in 1:3)))
        q = truth_component(books,predicate;oracle_id="q4/truth3/$(f)")
        @test exhaustive_check(q,predicate;oracle_id="q4/truth3/$(f)",profile=:indicator_exact).status===:pass
        @test count(b->b.role!==:primary,q.bits)<=1
    end
end

@testitem "Q4 compact table compilation and exactness boundaries" tags=[:q4] begin
    for mode in (:one_hot,:zero_one_hot,:domain_wall,:unary,:binary,:bounded_binary,:gray,:arithmetic,:bounded_coefficient,:fibonacci)
        book = structured_codebook(:x,0:3;encoding=mode)
        predicate = v->v[1] in (1,3)
        q = truth_component([book],predicate;oracle_id="q4/compact/$(mode)")
        @test exhaustive_check(q,predicate;oracle_id="q4/compact/$(mode)",profile=:indicator_exact).status===:pass
    end
    for n in 4:5, predicate in (v->isodd(sum(v)),v->all(==(1),v),v->sum(v)<=2,v->v[1]+v[2]==v[3])
        books = [structured_codebook(Symbol(:x,i),0:1;encoding=:native) for i in 1:n]
        q = truth_component(books,predicate;oracle_id="q4/higher/$(n)")
        @test exhaustive_check(q,predicate;oracle_id="q4/higher/$(n)",profile=:indicator_exact).status===:pass
    end
    # Fixed invalid-code energy=1 is stricter than the ordinary invalid-code >=1 contract.
    books = [structured_codebook(Symbol(:x,i),0:1;encoding=:one_hot) for i in 1:2]
    predicate = v->v[1]!=v[2]
    direct = local_constraint(:all_different,[0:1,0:1])
    fixed = truth_component(books,predicate;oracle_id="q4/overconstrained-extension")
    @test all(b->b.role===:primary,direct.bits)
    @test any(b->b.role!==:primary,fixed.bits)
    @test exhaustive_check(direct,predicate;oracle_id="q4/direct",profile=:indicator_exact).status===:pass
    @test exhaustive_check(fixed,predicate;oracle_id="q4/fixed",profile=:indicator_exact).status===:pass
    @test exhaustive_check(truth_component([],v->isempty(v);oracle_id="empty"),v->true;oracle_id="empty").status===:pass
end

@testitem "Q4 budgets never imply general infeasibility" tags=[:q4] begin
    books = [structured_codebook(Symbol(:x,i),0:1;encoding=:native) for i in 1:3]
    calls = Ref(0)
    @test_throws CompilationLimit truth_component(books,v->(calls[]+=1;true);oracle_id="budget",max_primary_states=4)
    @test calls[]==0
    @test_throws CompilationLimit truth_component(books,v->!all(==(1),v);oracle_id="budget",max_auxiliaries=0)
    @test_throws CompilationLimit truth_component(books,v->!all(==(1),v);oracle_id="budget",zero_aux_only=true)
    @test_throws CompilationLimit truth_component(books,v->isodd(sum(v));oracle_id="budget",max_terms=1)
    @test_throws CompilationLimit truth_component(books,v->isodd(sum(v));oracle_id="budget",max_work=1)
    @test_throws ArgumentError truth_component(books,v->true;oracle_id="")
    @test_throws ArgumentError truth_component(books,v->1;oracle_id="nonbool")
    @test_throws ArgumentError truth_component(books,v->true;oracle_id="budget",max_auxiliaries=-1)
end
