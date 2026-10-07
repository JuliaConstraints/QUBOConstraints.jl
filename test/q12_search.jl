@testitem "Q12 learned decision representation" tags=[:q12] begin
    using QUBOConstraints
    N=IntensionNode
    books=[structured_codebook(s,0:1) for s in (:x,:y)]
    space=AtomicSearchSpace(books;node_slots=16,constant_bound=8)
    @test fieldnames(AtomicSearchSpace)==(:codebooks,:node_slots,:constant_bound,:operators,
        :max_nodes,:max_local_rows,:max_total_rows,:max_value_bits,:max_depth)
    for expression in (N(:eq,:x,:y),N(:ne,:x,:y),N(:eq,N(:add,:x,:y),1),
            N(:eq,N(:div,1,:x),1),N(:eq,N(Symbol("if"),N(:eq,:x,0),7,N(:div,7,:x)),7))
        plan=core_atomic_plan(books,:intension;expression)
        weights=encode_atomic_weights(space,plan)
        domains=atomic_decision_domains(space)
        @test all(lo<=v<=hi for (v,(lo,hi)) in zip(weights,domains))
        # A fresh model and plain integers suffice, with no original recipe retained.
        fresh=AtomicSearchSpace(deepcopy(books);node_slots=16,constant_bound=8)
        restored=parse.(BigInt,string.(weights))
        @test canonical_equal(compile_atomic(plan).component,compile_atomic_weights(fresh,restored).component)
        @test restored==weights
    end
    # The SAME fixed structure selects genuinely different relations through weights.
    tiny=AtomicSearchSpace(books;node_slots=1,constant_bound=0,operators=(:eq,:ne))
    eqweights=BigInt[1,1,2,0,0,3]
    newweights=copy(eqweights); newweights[1]=2
    assignments=[Dict(:x=>x,:y=>y) for x in 0:1 for y in 0:1]
    eqplan=decode_atomic_weights(tiny,eqweights); neplan=decode_atomic_weights(tiny,newweights)
    @test [atomic_values(eqplan,d)[eqplan.root] for d in assignments]==[1,0,0,1]
    @test [atomic_values(neplan,d)[neplan.root] for d in assignments]==[0,1,1,0]
    aliased=copy(eqweights); aliased[3]=1
    aliasplan=decode_atomic_weights(tiny,aliased)
    @test all(atomic_values(aliasplan,d)[aliasplan.root]==1 for d in assignments)
    rooted=copy(eqweights); rooted[end]=1
    rootplan=decode_atomic_weights(tiny,rooted)
    @test [atomic_values(rootplan,d)[rootplan.root] for d in assignments]==[0,0,1,1]
    for changed in (BigInt[1,3,2,0,0,3],BigInt[1,1,2,1,0,3],BigInt[0,0,0,0,0,3],
            BigInt[1,1,2,0,1,3],BigInt[3,1,2,0,0,3])
        @test_throws ArgumentError decode_atomic_weights(tiny,changed)
    end
    @test_throws DimensionMismatch decode_atomic_weights(tiny,eqweights[1:5])
    @test_throws ArgumentError AtomicSearchSpace(books;node_slots=1,operators=(:knapsack,))
    @test_throws ArgumentError AtomicSearchSpace(books;node_slots=1,constant_bound=-1)
    @test_throws CompilationLimit encode_atomic_weights(AtomicSearchSpace(books;node_slots=0),eqplan)
    literalspace=AtomicSearchSpace(books;node_slots=2,constant_bound=1)
    literalplan=atomic_plan(books,N(:eq,:x,2);semantic_id="literal")
    @test_throws CompilationLimit encode_atomic_weights(literalspace,literalplan)
    foreign=[structured_codebook(s,0:1;encoding=:domain_wall) for s in (:x,:y)]
    @test_throws ArgumentError encode_atomic_weights(AtomicSearchSpace(foreign;node_slots=16),eqplan)
    # The public learner actually ranges over these decisions, without a witness seed.
    labels=Bool[d[:x]==d[:y] for d in assignments]
    learned=QUBOConstraints.train(tiny,assignments,labels;optimizer=AtomicEnumerativeOptimizer(max_candidates=300))
    @test learned.status===:training_fit && learned.mismatches==0
    @test learned.weights!==nothing && learned.candidates>1
    warm=QUBOConstraints.train(tiny,assignments,labels;initial_weights=eqweights,
        optimizer=AtomicEnumerativeOptimizer(max_candidates=1))
    @test warm.status===:training_fit && warm.weights==eqweights && warm.candidates==1
    @test_throws ArgumentError QUBOConstraints.train(tiny,assignments,labels;initial_weights=BigInt[100])
    learnedmatrix=compile_atomic_weights(tiny,learned.weights)
    @test exhaustive_check(learnedmatrix.component,v->v[1]==v[2];oracle_id="Q12/learned-equality",
        profile=:indicator_exact,max_states=1_048_576).status===:pass
    short=QUBOConstraints.train(tiny,assignments,labels;optimizer=AtomicEnumerativeOptimizer(max_candidates=1))
    @test short.status===:budget_exhausted
    impossible=QUBOConstraints.train(tiny,[assignments[1],assignments[1]],Bool[true,false];
        optimizer=AtomicEnumerativeOptimizer(max_candidates=300))
    @test impossible.status===:space_exhausted && impossible.mismatches==1
    @test_throws ArgumentError QUBOConstraints.train(tiny,[Dict(:x=>3,:y=>0)],Bool[false])
end
