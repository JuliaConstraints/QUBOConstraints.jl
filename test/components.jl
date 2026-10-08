@testitem "Q1 canonical polynomial and codebooks" tags=[:q1] begin
    using LinearAlgebra
    x, y = BitID(:x, 1), BitID(:y, 1)
    q = QUBOComponent([y, x]; linear = [(1, 2), (2, -1), (1, -2)],
        quadratic = [(1, 2, 4), (2, 1, -1), (2, 2, 2)], offset = -3)
    expected = QUBOComponent([x, y]; linear = [(1, 1)], quadratic = [(1, 2, 3)], offset = -3)
    @test canonical_equal(q, expected)
    d = dense_qubo(q)
    for z in ([false,false], [false,true], [true,false], [true,true])
        @test energy(q, z) == -3 + z[1] + 3z[1]*z[2]
        @test energy(q, z) == dot(z, d.matrix*z) + d.offset
    end
    @test_throws ArgumentError energy(q, [2, 0])
    @test_throws DimensionMismatch energy(q, [true])
    @test_throws ArgumentError QUBOComponent([x, x])
    @test_throws ArgumentError QUBOComponent([x]; linear = [(2, 1)])
    @test_throws ArgumentError QUBOComponent([x]; offset = Inf, coefficient_type = Float64)
    huge = QUBOComponent([x]; linear = [(1, typemax(Int)), (1, 1)])
    @test energy(huge, [true]) == big(typemax(Int)) + 1
    rational = QUBOComponent([x]; linear = [(1, 1//3)], coefficient_type = Rational{BigInt})
    @test energy(rational, [true]) == 1//3
    @test_throws ArgumentError BitID(:x, 0)
    @test_throws ArgumentError BitID(:x, 1; role = :embedding)
    for enc in (:one_hot, :domain_wall), values in ([3, 9, 17], [typemax(Int)], ["a", "b"])
        book = codebook(:v, values; encoding = enc)
        for value in values
            decoded = decode_code(book, encode_code(book, value))
            @test decoded.valid && decoded.value == value
        end
        @test !decode_code(book, fill(2, length(book.bits) + 1)).valid
    end
    redundant = Codebook(:x, [x], [7], reshape([false, true], 1, 2), [1, 1])
    @test size(encoded_codes(redundant, 7)) == (1, 2)
    @test decode_code(redundant, [true]).value == 7
    @test_throws ArgumentError encode_code(redundant, 8)
    @test_throws ArgumentError Codebook(:x, [x], [0, 1], [false false], [1,2])
    @test_throws ArgumentError Codebook(:x, [x], [0, 1], reshape([false], 1, 1), [1])
    for n in 1:5
        bits = reverse([BitID(:v,i) for i in 1:n])
        terms = [(i,j,mod(3i+j,7)-3) for i in 1:n for j in 1:n]
        candidate = QUBOComponent(bits; quadratic=terms,offset=2)
        for mask in 0:(2^n-1)
            z = [!iszero(mask & (1 << (i-1))) for i in 1:n]
            original = reverse(z)
            expected_polynomial = 2 + sum(v*original[i]*original[j] for (i,j,v) in terms)
            @test energy(candidate,z) == expected_polynomial
        end
    end
    floatq = QUBOComponent([x,y]; linear=[(1,2.0)],coefficient_type=Float64)
    @test (@inferred QUBOComponent{Float64}([x,y]; linear=[(1,2.0)])) isa QUBOComponent{Float64}
    @test (@inferred energy(floatq,[true,false])) == 2.0
    @test (@inferred energy(floatq,[1,0])) == 2.0
    @test (@inferred energy(floatq,[1.0,0.0])) == 2.0
    @testset "Binary numeric input preserves canonical energy" begin
        for T in (BigInt, Rational{BigInt}, Float64)
            numeric = QUBOComponent{T}([y,x]; linear=[(1,-2),(2,3)],
                quadratic=[(1,2,4),(2,1,-1),(1,1,2)], offset=-3)
            for mask in 0:3
                z = [!iszero(mask & (1 << (i-1))) for i in 1:2]
                expected = energy(numeric,z)
                for R in (Int, BigInt, Float64, Rational{Int})
                    input = R.(z)
                    @test energy(numeric,input) == expected
                    @test input == R.(z)
                    @test energy(numeric,view(input,:)) == expected
                end
            end
        end
        @test_throws ArgumentError energy(floatq,[0.5,1.0])
        @test_throws ArgumentError energy(floatq,[NaN,0.0])
        @test_throws ArgumentError energy(floatq,[Inf,0.0])
        @test_throws DimensionMismatch energy(floatq,[1,0,1])
        emptyq = QUBOComponent{Float64}(BitID[];offset=2)
        @test energy(emptyq,Int[]) == 2.0
    end
end

@testitem "Sorted sparse views retain canonical reductions and owned storage" begin
    using SparseArrays
    for T in (Float64,BigInt,Rational{BigInt}),width in (0,1,4,32)
        bits=reverse([BitID(:view_fixture,i) for i in 1:width])
        large=T===Float64 ? T(1e16) : T(big(2)^80)
        half=T===BigInt ? T(3) : T(3//2)
        small=T===BigInt ? T(1) : T(1//10)
        linear=Tuple{Int,T}[(i,large) for i in 1:width]
        append!(linear,[(i,-large) for i in 1:width])
        append!(linear,[(i,T(1)) for i in 1:width])
        quadratic=Tuple{Int,Int,T}[(i,i,half) for i in 1:width]
        for i in 1:width-1
            append!(quadratic,[(i,i+1,large),(i+1,i,-large),(i,i+1,small)])
        end
        # Independently remap, order and reduce the original coefficients. In
        # Float64 this includes cancellation that depends on canonical order.
        remap=Dict(bit=>i for (i,bit) in enumerate(sort(bits)))
        terms=[(remap[bits[i]],remap[bits[i]],v) for (i,v) in linear]
        append!(terms,[(minmax(remap[bits[i]],remap[bits[j]])...,v) for (i,j,v) in quadratic])
        sort!(terms;by=t->(t[2],t[1],t[3]))
        coefficients=Dict{Tuple{Int,Int},T}()
        for (i,j,v) in terms
            key=(i,j);coefficients[key]=haskey(coefficients,key) ? coefficients[key]+v : v
        end
        q=QUBOComponent{T}(bits;linear,quadratic,offset=2)
        @test q.bits==sort(bits)
        @test q.offset==T(2)
        for i in 1:width,j in i:width
            @test (i==j ? q.linear[i] : q.quadratic[i,j])==get(coefficients,(i,j),zero(T))
        end
        @test all(i<j for j in 1:width for i in rowvals(q.quadratic)[nzrange(q.quadratic,j)])
        @test all(!iszero,nonzeros(q.linear)) && all(!iszero,nonzeros(q.quadratic))
        before_linear,before_quadratic=copy(q.linear),copy(q.quadratic)
        empty!(linear);empty!(quadratic);empty!(bits)
        @test q.linear==before_linear && q.quadratic==before_quadratic && length(q.bits)==width
    end
    @test_throws ArgumentError QUBOComponent{Float64}([BitID(:overflow,1)];linear=[(1,1e308),(1,1e308)])
    @test_throws ArgumentError QUBOComponent{Float64}([BitID(:nonfinite,1)];quadratic=[(1,1,NaN)])
end

@testitem "Q1 owned heterogeneous codebook coverage" tags=[:q1] begin
    import QUBOConstraints as QC
    struct TupleBitsCodebook <: QC.AbstractCodebook
        variable::Symbol
        bits::Tuple{BitID,BitID}
        values::Vector{Int}
    end
    for T in (BigInt, Rational{BigInt}, Float64)
        explicit = codebook(:a, [-3, 5])
        structured = structured_codebook(:b, [-3, 5]; encoding=:one_hot)
        books = QC.AbstractCodebook[explicit, structured]
        bits = reverse([explicit.bits; structured.bits])
        linear = [(1, 2), (2, -3), (4, 1)]
        quadratic = [(1, 2, 4), (2, 1, -1), (3, 3, 2)]
        q = QUBOComponent{T}(bits; linear, quadratic, offset=2, codebooks=books)
        @test q.codebooks !== books
        for (owned, input) in zip(q.codebooks, books)
            @test typeof(owned) === typeof(input)
            @test owned !== input
            @test owned.bits == input.bits && owned.bits !== input.bits
            @test owned.values == input.values && owned.values !== input.values
        end
        for mask in 0:15
            original = [!iszero(mask & (1 << (i-1))) for i in 1:4]
            z = [original[findfirst(==(bit), bits)] for bit in q.bits]
            expected = T(2) + sum(T(v)*original[i] for (i,v) in linear) +
                sum(T(v)*original[i]*original[j] for (i,j,v) in quadratic)
            @test energy(q,z) == expected
        end
        explicit.bits[1] = BitID(:changed,1)
        explicit.values[1] = -9
        explicit.codes[1,1] = false
        @test q.codebooks[1].bits == [BitID(:a,1),BitID(:a,2)]
        @test q.codebooks[1].values == [-3,5]
        @test q.codebooks[1].codes[1,1]
    end
    book = codebook(:a, [-3,5])
    alias = Codebook(:alias, book.bits, book.values, book.codes, book.value_indices)
    @test_throws ArgumentError QUBOComponent(book.bits;codebooks=[book,book])
    @test_throws ArgumentError QUBOComponent(book.bits;codebooks=[book,alias])
    @test_throws ArgumentError QUBOComponent([book.bits[1]];codebooks=[book])
    singleton = codebook(:constant,[17];encoding=:domain_wall)
    q = QUBOComponent(BitID[];codebooks=[singleton],offset=3)
    @test energy(q,Bool[]) == 3
    @test q.codebooks[1].values == [17]
    subset = QUBOComponent([book.bits;BitID(:unbound,1)];codebooks=[book])
    @test subset.codebooks[1].bits == book.bits
    custom = TupleBitsCodebook(:custom,(BitID(:custom,1),BitID(:custom,2)),[0,1])
    qcustom = QUBOComponent(collect(custom.bits);codebooks=[custom])
    @test qcustom.codebooks[1].bits == custom.bits
    @test qcustom.codebooks[1].values !== custom.values
end

@testitem "Q1 exhaustive oracle, auxiliaries and invalid encodings" tags=[:q1] begin
    book = codebook(:x, 0:1; encoding = :domain_wall)
    x = only(book.bits)
    aux = BitID(:witness, 1; role = :semantic_auxiliary)
    # E = (x-a)^2 + x; minimization over a gives x exactly.
    q = QUBOComponent([x, aux]; linear = [(1,2), (2,1)], quadratic = [(1,2,-2)],
        codebooks = [book], auxiliary_meanings = Dict(aux => "copy of x"))
    r = exhaustive_check(q, v -> v[1] == 0; oracle_id = "x-eq-zero", profile = :indicator_exact)
    @test r.status === :pass && r.proof === :exhaustive
    @test r.evaluations == 4 && r.valid_codes == 2 && r.invalid_codes == 0
    @test r.min_violation == 1
    @test exhaustive_check(q, _ -> true; oracle_id = "always").status === :fail
    @test exhaustive_check(q, _ -> 1; oracle_id = "bad-oracle").status === :failed
    @test exhaustive_check(q, _ -> error("oracle failed"); oracle_id = "bad-oracle").status === :failed
    @test exhaustive_check(q, _ -> true; oracle_id = "budget", max_states = 2).status === :unknown
    @test exhaustive_check(q, _ -> true; oracle_id = "timeout", time_limit = 1e-12).status === :unknown
    @test exhaustive_check(q, _ -> true; oracle_id = "surrogate", profile = :surrogate).status === :surrogate
    @test_throws ArgumentError exhaustive_check(q, _ -> true; oracle_id = "gap", gap=0)
    @test_throws ArgumentError exhaustive_check(q, _ -> true; oracle_id = "gap", gap=1, atol=0.5)
    near = QUBOComponent(book.bits; linear=[(1, (big(2)^54-1)//big(2)^54)],
        coefficient_type=Rational{BigInt},codebooks=[book])
    @test exhaustive_check(near, v -> only(v)==0; oracle_id="strict-gap",gap=1.0).status === :fail
    oh = codebook(:v, [5, 9])
    # (1-u-v)^2 protects both invalid one-hot codes.
    valid = QUBOComponent(oh.bits; offset=1, linear=[(1,-1),(2,-1)],
        quadratic=[(1,2,2)], codebooks=[oh])
    vr = exhaustive_check(valid, _ -> true; oracle_id="one-hot-validity")
    @test vr.status === :pass && vr.invalid_codes == 2
    invalid = QUBOComponent(oh.bits; codebooks=[oh])
    ir = exhaustive_check(invalid, _ -> true; oracle_id="missing-validity")
    @test ir.status === :fail && !ir.counterexample.valid_code
    @test ir.counterexample.primary_mask == 0
    # A redundant value must have correct energy for BOTH codes to compose safely.
    redundant = Codebook(:x, [x], [7], [false true], [1,1])
    mismatch = QUBOComponent([x]; linear=[(1,1)], codebooks=[redundant])
    @test exhaustive_check(mismatch, _ -> true; oracle_id="redundancy").status === :fail
    constant = QUBOComponent(BitID[]; codebooks=[codebook(:c, [42]; encoding=:domain_wall)])
    @test exhaustive_check(constant, v -> only(v)==42; oracle_id="singleton").evaluations == 1
    floatq = QUBOComponent(book.bits; linear=[(1,1.0)], coefficient_type=Float64, codebooks=[book])
    @test exhaustive_check(floatq, v -> only(v)==0; oracle_id="float").proof === :validated
    @test_throws ArgumentError exhaustive_check(QUBOComponent([x]), _ -> true; oracle_id="missing-book")
end

@testitem "Q1 composition and certificate artifacts" tags=[:q1] begin
    using TOML
    book = codebook(:x, 0:1; encoding=:domain_wall)
    x = only(book.bits)
    a = BitID(:a, 1; role=:quadratization_auxiliary)
    q = QUBOComponent([x,a]; linear=[(1,2),(2,1)], quadratic=[(1,2,-2)], codebooks=[book])
    combined = compose(q,q)
    @test length(combined.bits) == 3
    @test length(combined.auxiliary_meanings) == 2
    r = exhaustive_check(combined, v -> only(v)==0; oracle_id="double-zero")
    @test r.status === :pass && r.min_violation == 2
    @test exhaustive_check(combined, v -> only(v)==0; oracle_id="double-zero", profile=:indicator_exact).status === :fail
    different = codebook(:x, [1,0]; encoding=:domain_wall)
    @test_throws ArgumentError compose(q, QUBOComponent(different.bits; codebooks=[different]))
    @test_throws ArgumentError compose(rename_auxiliaries(q,:right), q)
    data = component_artifact(combined; report=r)
    @test data["schema_version"] == "qubo-component/1"
    @test data["certificate"]["proof"] == "exhaustive"
    @test length(data["sha256"]) == 64
    text = sprint(io -> write_component(io, combined; report=r))
    @test TOML.parse(text)["sha256"] == data["sha256"]
    @test_throws ArgumentError component_artifact(q; report=r)
    bad = exhaustive_check(q, _ -> true; oracle_id="bad")
    @test haskey(component_artifact(q; report=bad)["certificate"], "counterexample")
    # Altered codebooks must invalidate old evidence even with an identical polynomial.
    combined.codebooks[1].values[1] = 100
    @test_throws ArgumentError component_artifact(combined; report=r)
end
@testitem "Composition streams exact owned sparse terms" tags=[:q1] begin
    rename(bit)=bit.role===:primary ? bit :
        BitID(Symbol(:right,"/",bit.owner),bit.index;role=bit.role)
    for T in (BigInt,Rational{BigInt},Float64),S in (BigInt,Rational{BigInt},Float64),width in (0,1,2,8)
        books=[codebook(Symbol(:x,i),0:1;encoding=:domain_wall) for i in 1:width]
        primary=BitID[b for book in books for b in book.bits]
        bits=vcat(primary,[BitID(:aux,i;role=:quadratization_auxiliary) for i in 1:width])
        left_linear=[(i,T(mod(i,5)-2)) for i in eachindex(bits)]
        right_linear=[(i,S(mod(i,7)-3)) for i in eachindex(bits)]
        left_quadratic=[(i,i+1,T(2)) for i in 1:length(bits)-1]
        right_quadratic=[(i,i+1,S(-1)) for i in 1:length(bits)-1]
        a=QUBOComponent{T}(bits;linear=left_linear,quadratic=left_quadratic,offset=1,
            codebooks=books,auxiliary_meanings=Dict(bit=>"left" for bit in bits if bit.role!==:primary),
            applicability="left",provenance="left")
        b=QUBOComponent{S}(reverse(bits);linear=right_linear,quadratic=right_quadratic,offset=-2,
            codebooks=books,auxiliary_meanings=Dict(bit=>"right" for bit in bits if bit.role!==:primary),
            applicability="right",provenance="right")
        combined=compose(a,b)
        R=T===Float64 || S===Float64 ? Float64 : promote_type(T,S)
        @test combined isa QUBOComponent{R}
        @test length(combined.bits)==3width
        @test combined.offset==R(-1)
        @test combined.applicability=="(left) AND (right)"
        @test combined.provenance=="add(left, right)"
        @test all(combined.auxiliary_meanings[bit]=="left" for bit in a.bits if bit.role!==:primary)
        @test all(combined.auxiliary_meanings[rename(bit)]=="right" for bit in b.bits if bit.role!==:primary)
        @test length(combined.codebooks)==width
        @test combined.bits!==a.bits && combined.bits!==b.bits
        @test combined.linear!==a.linear && combined.quadratic!==b.quadratic
        @test all(zip(books,combined.codebooks)) do pair
            input,owned=pair
            input!==owned && input.values==owned.values && input.values!==owned.values &&
                input.codes==owned.codes && input.codes!==owned.codes
        end
        positions=Dict(bit=>i for (i,bit) in enumerate(combined.bits))
        patterns=width<=2 ? [[isodd(mask>>(i-1)) for i in eachindex(combined.bits)]
            for mask in 0:(1<<length(combined.bits))-1] :
            [[pattern(i) for i in eachindex(combined.bits)] for pattern in
                (i->false,i->true,isodd,i->mod(i,3)==0)]
        for z in patterns
            left=[z[positions[bit]] for bit in bits]
            right=[z[positions[rename(bit)]] for bit in reverse(bits)]
            expected=T(1)+sum(v*left[i] for (i,v) in left_linear;init=zero(T))+
                sum(v*left[i]*left[j] for (i,j,v) in left_quadratic;init=zero(T))+
                S(-2)+sum(v*right[i] for (i,v) in right_linear;init=zero(S))+
                sum(v*right[i]*right[j] for (i,j,v) in right_quadratic;init=zero(S))
            @test energy(combined,z)==expected
        end
        artifact=component_artifact(combined)
        if width>0
            a.codebooks[1].values[1]=99
            a.auxiliary_meanings[only(filter(bit->bit.index==1 && bit.role!==:primary,a.bits))]="changed"
            @test component_artifact(combined)==artifact
        end
    end
    bits=[BitID(:fraction,1)]
    rational=QUBOComponent{Rational{BigInt}}(bits;linear=[(1,1//3)])
    float=QUBOComponent{Float64}(bits)
    @test energy(compose(float,rational),[true])==Float64(1//3)
    @test energy(compose(rational,float),[true])==Float64(1//3)
end
