@testitem "Q2 candidate motifs and validity" tags=[:q2] begin
    include(joinpath(@__DIR__,"..","perf","q2","patterns.jl"))
    using LinearAlgebra
    square(k) = [square_pattern(k,a,b,2,4,5,3) for a in 1:5,b in 1:5]
    @test square(1)==Matrix{Int}(I,5,5)
    @test square(2)==ones(Int,5,5)-Matrix{Int}(I,5,5)
    @test square(3)==[Int(a>b) for a in 1:5,b in 1:5]
    @test square(4)==[Int(a>=b) for a in 1:5,b in 1:5]
    @test square(5)==[Int(a<b) for a in 1:5,b in 1:5]
    @test square(6)==[Int(a<=b) for a in 1:5,b in 1:5]
    @test square(7)==[0 0 0 -1 0; -1 -1 -1 -1 -1; 0 0 0 -1 0; 0 0 0 -1 0; 0 0 0 -1 0]
    @test square(8)==-square(7)
    @test square(9)==[-1 0 -1 0 -1; 0 0 0 0 0; -1 0 -1 0 -1; 0 0 0 0 0; -1 0 -1 0 -1]
    @test square(10)==-square(9)
    @test square(11)==[0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;1 -1 1 1 1;0 1 0 0 0]
    @test square(12)==[3 2 1 0 0;2 3 2 1 0;1 2 3 2 1;0 1 2 3 2;0 0 1 2 3]
    @test square(13)==[0 0 1 2 3;0 0 0 1 2;1 0 0 0 1;2 1 0 0 0;3 2 1 0 0]
    @test square(14)==[2a*b for a in 1:5,b in 1:5]
    @test length(selections())==756
    for triangle in 1:3
        q = raw_pattern_component(Int[],triangle,1,3,4)
        for mask in 0:7
            z = [!iszero(mask & (1<<(i-1))) for i in 1:3]
            base = 2sum(z[i]*z[j] for i in 1:3 for j in i+1:3)
            expected = triangle==1 ? base : triangle==2 ? base-sum(z) : base+(sum(i*z[i] for i in 1:3)-4)^2-16
            @test energy(q,z)==expected
        end
    end
    raw = raw_pattern_component([7,9,11],2,3,3,0;offset=6)
    @test energy(raw,[1,1,0,1,1,0,0,0,1])==-1
    safe = guard_one_hot(raw)
    @test exhaustive_check(safe.component,v->all(i->v[v[i]]==i,1:3);oracle_id="q2/guarded-atoms").status===:pass
end
