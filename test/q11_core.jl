@testitem "Q11 core completeness regressions" tags=[:q11] begin
    using QUBOConstraints
    N=IntensionNode; C=AtomicCondition
    books=[structured_codebook(s,0:1) for s in (:x,:y)]
    cases=[
        (N(:eq,N(:div,:x,N(:add,N(:sub,:x,:x),1)),:x),(x,y)->true),
        (N(:eq,N(Symbol("if"),N(:eq,:x,0),7,N(:div,7,:x)),7),(x,y)->true),
        (N(:eq,N(:div,1,:x),1),(x,y)->x==1),
        (N(:eq,N(:mod,1,:x),0),(x,y)->x==1),
        (N(:eq,N(:pow,-1,N(:sub,0,:x)),1),(x,y)->x==0),
        (N(:eq,N(Symbol("if"),N(:eq,:x,0),1,N(:pow,2,-1)),1),(x,y)->x==0),
        (N(:eq,N(:pow,2,-1),0),(x,y)->false),
    ]
    for (e,oracle) in cases
        p=core_atomic_plan(books,:intension;expression=e)
        @test validate_atomic_plan(p)
        q=compile_atomic(p)
        for x in 0:1,y in 0:1
            @test energy(q.component,atomic_witness(q,Dict(:x=>encode_code(books[1],x),:y=>encode_code(books[2],y))))==!oracle(x,y)
        end
    end
    # Exhaust the auxiliary cube of the smallest guarded relation independently.
    b=[structured_codebook(:x,0:1;encoding=:native)]
    q=compile_atomic(atomic_plan(b,N(:div_total,1,:x);semantic_id="Q11/totalized-div-atom"))
    @test exhaustive_check(q.component,v->only(v)==1;oracle_id="Q11/defined-div",profile=:indicator_exact).status===:pass
    # Distinct runtime values are part of the all-variable distribute form.
    p=core_atomic_plan(books,:cardinality;list=[:x,:x],values=[:y,:y],occurs=[:x,:x])
    @test atomic_values(p,Dict(:x=>0,:y=>1))[p.root]==0
    unguarded=core_atomic_plan(books,:cardinality;list=[:x,:x],values=[:y,:y],occurs=[:x,:x],distinct_values=false)
    @test atomic_values(unguarded,Dict(:x=>0,:y=>1))[unguarded.root]==1
    legacy=core_atomic_plan(books,:cardinality;list=[:x,:y],values=[:x,:y],occurs=[2,2],distinct_values=false)
    @test atomic_values(legacy,Dict(:x=>1,:y=>1))[legacy.root]==1
    for family in (:noOverlap,:cumulative,:binPacking)
        kwargs=family===:noOverlap ? (;origins=[:x,:y],lengths=[-1,1]) : family===:cumulative ?
            (;origins=[:x,:y],lengths=[-1,1],heights=[1,1],condition=C(:le,1)) : (;sizes=[-1,1],condition=C(:le,1))
        p=core_atomic_plan(books,family;parameter_signs=:algebraic,kwargs...)
        @test validate_atomic_plan(p)
    end
    p=core_atomic_plan(books,:mdd;transitions=[(:s,0,:a),(:a,1,:t)])
    @test atomic_values(p,Dict(:x=>0,:y=>1))[p.root]==1
    @test_throws ArgumentError core_atomic_plan(books,:mdd;transitions=[(:s,0,:a),(:a,1,:s)])
    p=core_atomic_plan(books,:ordered;form=:lists,lists=[[:x,:y],[1,1]],operator=:le)
    @test atomic_values(p,Dict(:x=>1,:y=>0))[p.root]==1
    p=core_atomic_plan(books,:element;value=:y)
    @test atomic_values(p,Dict(:x=>0,:y=>1))[p.root]==1
    named=[structured_codebook(:atomic,0:1)]
    c=compile_atomic(core_atomic_plan(named,:intension;expression=:atomic))
    @test energy(c.component,atomic_witness(c,Dict(:atomic=>encode_code(named[1],1))))==0
end
