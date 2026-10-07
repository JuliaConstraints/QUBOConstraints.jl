@testitem "Q5 separable sums and all primary representations" tags=[:q5] begin
    modes = (:native,:one_hot,:zero_one_hot,:domain_wall,:unary,:binary,:bounded_binary,:arithmetic,:bounded_coefficient,:fibonacci)
    for mode in modes
        domains = mode===:native ? [[-2,3],[0,2]] : [[-2,0,2],[1,3,5]]
        books = [structured_codebook(Symbol(:x,i),d;encoding=mode) for (i,d) in enumerate(domains)]
        for op in (:eq,:le,:lt,:ge,:gt), rhs in (-3,0,4), slack in (:unary,:binary,:arithmetic,:bounded_coefficient,:fibonacci)
            condition = IntegerCondition(op,rhs)
            q = sum_component(books,condition;coefficients=[2,-1],slack_encoding=slack)
            @test exhaustive_check(q,v->satisfies(condition,2v[1]-v[2]);oracle_id="q5/sum/$(mode)/$(op)/$(rhs)/$(slack)").status===:pass
        end
    end
end

@testitem "Q5 condition sets and membership counts" tags=[:q5] begin
    books = [structured_codebook(Symbol(:x,i),[0,2,5];encoding=:domain_wall) for i in 1:2]
    for condition in (IntegerCondition(:ne,2),IntegerCondition(:in,[1,2,5,7]),IntegerCondition(:notin,2:5),IntegerCondition(:in,Int[]))
        q = sum_component(books,condition;slack_encoding=:one_hot)
        @test exhaustive_check(q,v->satisfies(condition,sum(v));oracle_id="q5/set").status===:pass
        q = count_component(books,[0,5],condition;slack_encoding=:domain_wall)
        @test exhaustive_check(q,v->satisfies(condition,count(in((0,5)),v));oracle_id="q5/count").status===:pass
    end
    @test_throws ArgumentError IntegerCondition(:eq,[1])
    @test_throws ArgumentError IntegerCondition(:in,1)
    @test_throws CompilationLimit sum_component(books,IntegerCondition(:le,100);max_witness_values=2)
    @test_throws ArgumentError sum_component(books,IntegerCondition(:ne,2);slack_encoding=:binary)
    q = sum_component([structured_codebook(:x,[typemax(Int)])],IntegerCondition(:eq,2big(typemax(Int)));coefficients=[2])
    @test exhaustive_check(q,v->true;oracle_id="q5/overflow").status===:pass
end
