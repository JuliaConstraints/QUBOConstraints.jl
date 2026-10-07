@testitem "Q3 compact codecs and exact validity" tags=[:q3] begin
    using TOML
    modes = (:native,:one_hot,:zero_one_hot,:domain_wall,:unary,:binary,
        :bounded_binary,:gray,:arithmetic,:bounded_coefficient,:fibonacci)
    for mode in modes, n in 1:(mode===:native ? 2 : 8)
        book = structured_codebook(:x,collect(-3:n-4);encoding=mode)
        @test all(decode_code(book,encode_code(book,v)).valid &&
            decode_code(book,encode_code(book,v)).value==v for v in book.values)
        counts = zeros(Int,n)
        for mask in 0:2^length(book.bits)-1
            z = [!iszero(mask & (1<<(i-1))) for i in eachindex(book.bits)]
            d = decode_code(book,z)
            d.valid && (counts[d.value+4]+=1)
        end
        @test all(>(0),counts)
        if mode in (:native,:one_hot,:zero_one_hot,:domain_wall,:binary,:gray,:fibonacci)
            @test all(==(1),counts)
        elseif mode===:unary
            @test counts==[binomial(n-1,k) for k in 0:n-1]
        end
        @test all(size(encoded_codes(book,v),2)==counts[v+4] for v in book.values)
        q = encoding_validity(book)
        r = exhaustive_check(q,v->true;oracle_id="q3/validity/$(mode)/$(n)",max_states=1_048_576)
        @test r.status===:pass
        @test r.evaluations==2^length(q.bits)
        @test component_artifact(q;report=r)["schema_version"]=="qubo-component/2"
        @test TOML.parse(sprint(io->write_component(io,q;report=r)))["sha256"]==r.component_digest
    end
end

@testitem "Q3 boundaries and non-exponential representation" tags=[:q3] begin
    # Neighbours of powers of two, Fibonacci and triangular numbers, with overlap removed.
    boundaries = sort!(unique(vcat([max(1,k+d) for k in (2,4,8,16,32,64,128,256,512,1024,
        3,5,13,21,34,55,89,144,233,377,610,987,6,10,15,28,36,45,66,78,91) for d in -1:1])))
    for mode in (:binary,:bounded_binary,:gray,:arithmetic,:bounded_coefficient,:fibonacci), K in boundaries
        book = structured_codebook(:x,0:K;encoding=mode,coefficient_bound=3)
        @test all(decode_code(book,encode_code(book,v)).value==v for v in (0,K÷2,K))
        mode===:bounded_coefficient && @test maximum(book.weights)<=3
    end
    book = structured_codebook(:x,0:1000;encoding=:unary)
    @test Base.summarysize(book)<200_000
    @test count(encode_code(book,700))==700
    @test decode_code(book,encode_code(book,700)).value==700
    @test_throws ArgumentError encoded_codes(book,0)
    @test isempty(encoding_validity(book).quadratic.nzval)
    @test exhaustive_check(encoding_validity(book),v->true;oracle_id="budget").status===:unknown
    @test_throws ArgumentError structured_codebook(:x,[])
    @test_throws ArgumentError structured_codebook(:x,[1,1])
    @test_throws ArgumentError structured_codebook(:x,0:2;encoding=:native)
    @test_throws ArgumentError structured_codebook(:x,0:2;encoding=:mixed_radix)
    @test_throws ArgumentError structured_codebook(:x,0:2;coefficient_bound=0)
    @test_throws ArgumentError encoding_validity(structured_codebook(:x,1:5);max_terms=2)
    book = structured_codebook(:x,[nothing,"present"];encoding=:zero_one_hot)
    @test decode_code(book,[false]).valid && decode_code(book,[false]).value===nothing
    @test !decode_code(book,[2]).valid
    @test !decode_code(book,[false,false]).valid
    @test length(encoding_registry().entries)==18
end

@testitem "Q3 explicit codebook compatibility and proof binding" tags=[:q3] begin
    for mode in (:one_hot,:domain_wall), n in 1:7
        explicit = codebook(:x,1:n;encoding=mode)
        compact = structured_codebook(:x,1:n;encoding=mode)
        @test all(encode_code(explicit,v)==encode_code(compact,v) for v in 1:n)
        q = encoding_validity(compact)
        @test exhaustive_check(compose(q,q),v->true;oracle_id="compose").status===:pass
        mode===:one_hot && @test canonical_equal(one_hot_validity([compact]),q)
    end
    q = encoding_validity(structured_codebook(:x,0:3;encoding=:unary))
    bad = QUBOComponent(q.bits;linear=[(1,-1),(2,1)],codebooks=q.codebooks)
    # One zero energy representative is not sufficient in a redundant encoding.
    @test exhaustive_check(bad,v->true;oracle_id="all-representations").status===:fail
    r = exhaustive_check(q,v->true;oracle_id="valid")
    @test_throws ArgumentError component_artifact(bad;report=r)
    @test_throws ArgumentError compose(q,encoding_validity(structured_codebook(:x,0:3;encoding=:domain_wall)))
    mixed = compose(one_hot_constraint(:all_different,[0:1];variables=[:a]),
        encoding_validity(structured_codebook(:b,0:2;encoding=:fibonacci)))
    @test length(mixed.codebooks)==2
    @test exhaustive_check(mixed,v->true;oracle_id="mixed-codebooks").status===:pass
end

@testitem "Q3 independent small-code fixtures and rank method" tags=[:q3] begin
    gray = structured_codebook(:g,0:3;encoding=:gray)
    @test [collect(encode_code(gray,v)) for v in 0:3]==[[false,false],[true,false],[true,true],[false,true]]
    fib = structured_codebook(:f,0:8;encoding=:fibonacci)
    @test fib.weights==[1,2,3,5,8]
    @test encode_code(fib,7)==[false,true,false,true,false]
    @test !decode_code(fib,[true,true,false,false,false]).valid
    @test !decode_code(fib,[false,true,false,false,true]).valid
    for K in (1,3,7,15,31)
        b = structured_codebook(:x,0:K;encoding=:gray)
        q = encoding_validity(b)
        @test q.bits==sort(b.bits) && iszero(q.offset) && isempty(q.linear.nzval) && isempty(q.quadratic.nzval)
    end
    for K in (1,2,4,7,12,20,33)
        b = structured_codebook(:x,0:K;encoding=:fibonacci)
        @test encoding_validity(b).bits==sort(b.bits)
    end
    for (mode,weights) in ((:bounded_binary,[1,2,4,4]),(:arithmetic,[1,1,2,3,4]),(:bounded_coefficient,[1,2,2,3,3]))
        b = structured_codebook(:x,0:11;encoding=mode,coefficient_bound=3)
        @test b.weights==weights
    end
    include(joinpath(@__DIR__,"..","perf","q3","campaign.jl"))
    @test exact_rank([1 2 3;2 4 6])==1
    @test exact_rank([0 1 0;0 0 1;0 0 0])==2
    @test exact_rank(zeros(Int,3,4))==0
    @test exact_rank([1 0;0 1])==2
    for mode in (:one_hot,:domain_wall)
        data = valid_rows(structured_codebook(:x,0:3;encoding=mode))
        @test capabilities(data.rows,data.labels,4)["all_semantic_functions_on_valid_codes"]
    end
    data = valid_rows(structured_codebook(:x,0:3;encoding=:unary))
    @test !capabilities(data.rows,data.labels,4)["all_semantic_functions_on_valid_codes"]
    graph = graph_metrics([0,1,3,7],3)
    @test graph["components"]==1 && graph["largest_component_diameter"]==3
end
