@testitem "Q2 globally exact one-hot families" tags=[:q2] begin
    using SparseArrays
    channel = v->all(i->v[v[i]]==i,eachindex(v))
    @testset "Channel no spurious ground states" begin
        for n in 1:4
            q = one_hot_constraint(:channel,[collect(1:n) for _ in 1:n])
            r = exhaustive_check(q,channel;oracle_id="channel/involution-v1")
            @test r.status === :pass && r.evaluations == 2^(n*n)
            @test count(b->b.role!==:primary,q.bits) == 0
            @test maximum(abs,nonzeros(q.quadratic);init=big(0)) <= 2
        end
        q = one_hot_constraint(:channel,[1:3 for _ in 1:3])
        @test energy(q,[1,1,0,1,1,0,0,0,1]) == 2
        for domains in ([[1,2],[1,2],[3]], [[2],[1],[1,3]], [[2],[3],[1]])
            q = one_hot_constraint(:channel,domains)
            @test exhaustive_check(q,channel;oracle_id="channel/subdomains-v1").status === :pass
        end
    end
    @testset "Other initial families" begin
        cases = [
            (:all_different,[[1,3],[3,5],[1,5]],(;),v->allunique(v)),
            (:all_different,[["a","b"],["b","c"]],(;),v->allunique(v)),
            (:ordered,[[5,-2],[0,9],[1,10]],(;),v->issorted(v)),
            (:linear_sum,[[-2,3],[1,5],[-1,2]],(;rhs=4,coefficients=[2,-1,3]),v->2v[1]-v[2]+3v[3]==4),
            (:linear_sum,[[typemax(Int)],[1]],(;rhs=big(typemax(Int))+1),v->sum(big.(v))==big(typemax(Int))+1),
            (:linear_sum,[[0,1],[0,1]],(;rhs=9),v->sum(v)==9),
            (:no_overlap,[[0,2],[1,4],[3,6]],(;lengths=[2,1,3]),
                v->all(v[i]+[2,1,3][i]<=v[j] || v[j]+[2,1,3][j]<=v[i] for i in 1:3 for j in i+1:3)),
            (:no_overlap,[[0,1],[0,2]],(;lengths=[0,2]),v->v[1]<=v[2] || v[2]+2<=v[1])]
        for (family,domains,kwargs,predicate) in cases
            q = one_hot_constraint(family,domains;kwargs...)
            r = exhaustive_check(q,predicate;oracle_id="q2/$(family)-independent")
            @test r.status === :pass
            @test r.evaluations == 2^sum(length,domains)
        end
    end
    @testset "Validation and composition" begin
        @test_throws ArgumentError one_hot_constraint(:unknown,[[1]])
        @test_throws ArgumentError one_hot_constraint(:channel,[[0,1]])
        @test_throws ArgumentError one_hot_constraint(:no_overlap,[[1]];lengths=[-1])
        @test_throws ArgumentError one_hot_constraint(:linear_sum,[[1]];coefficients=[1.0])
        @test_throws ArgumentError one_hot_constraint(:channel,[[1],[2]];variables=[:x,:x])
        a = one_hot_constraint(:all_different,[1:3,1:3,1:3])
        b = one_hot_constraint(:ordered,[1:3,1:3,1:3])
        @test exhaustive_check(compose(a,b),v->allunique(v)&&issorted(v);oracle_id="q2/composition").status === :pass
        raw = QUBOComponent(a.bits;linear=[(i,-1) for i in 1:length(a.bits)],offset=3,codebooks=a.codebooks)
        guarded = guard_one_hot(raw)
        @test guarded.lower_bound == -6 && guarded.weight == 7
        @test exhaustive_check(guarded.component,v->true;oracle_id="q2/guard").status === :pass
        @test_throws ArgumentError one_hot_validity([codebook(:v,0:2;encoding=:domain_wall)])
        @test_throws ArgumentError guard_one_hot(raw;gap=0)
    end
end

@testitem "Q2 exact auxiliary calibration" tags=[:q2] begin
    include(joinpath(@__DIR__,"..","perf","q2","calibration.jl"))
    for n in 1:10, kind in (:positive_monomial,:negative_monomial_shifted,:even_parity)
        q = calibration_qubo(kind,n)
        report = exhaustive_check(q,v->calibration_truth(kind,v);oracle_id="q2/$(kind)",profile=:indicator_exact)
        @test report.status === :pass
        expected = kind===:negative_monomial_shifted ? Int(n>=3) : max(0,ndigits(n-1;base=2)-1)
        @test count(b->b.role!==:primary,q.bits)==expected
    end
end

@testitem "Q2 dense sum has global validity" tags=[:q2] begin
    using LinearAlgebra
    for n in 1:3, rhs in (0,n,2n)
        q = QUBO_linear_sum(n,rhs)
        @test eltype(q)===BigInt && istriu(q)
        correct = true
        for mask in 0:2^(n*n)-1
            z = [!iszero(mask & (1<<(i-1))) for i in 1:n*n]
            decoded_sum = sum(mod(i-1,n)*z[i] for i in eachindex(z))
            validity = sum((sum(z[(i-1)*n+1:i*n])-1)^2 for i in 1:n)
            correct &= dot(z,q*z)+big(rhs)^2+n == (decoded_sum-rhs)^2+validity
        end
        @test correct
    end
    @test_throws ArgumentError QUBO_linear_sum(0,0)
    @test QUBO_linear_sum(1,big(typemax(Int))+1)==reshape(BigInt[-1],1,1)
end
