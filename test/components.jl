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
@testitem "Canonical diagonal storage is exact and owned" tags=[:q1] begin
    using SparseArrays
    for T in (Float64,BigInt,Rational{BigInt}),n in (0,1,2,8,32),mode in (:none,:full,:holes,:cancelled)
        coefficient(k)=T===Rational{BigInt} ? T(k//3) : T(k)
        bits=reverse([BitID(:diagonal,i) for i in 1:n])
        linear=Tuple{Int,T}[]
        quadratic=[(i,j,coefficient(mod(i+j,5)-2)) for j in 1:n for i in 1:j-1]
        for i in 1:n
            if mode===:full || (mode===:holes && isodd(i))
                push!(linear,(i,coefficient(i)))
                push!(quadratic,(i,i,coefficient(-2)))
            elseif mode===:cancelled
                push!(linear,(i,coefficient(3)))
                push!(quadratic,(i,i,coefficient(-3)))
            end
        end
        expected=zeros(T,n,n)
        for (i,v) in linear;expected[n-i+1,n-i+1]+=v;end
        for (i,j,v) in quadratic
            row,column=minmax(n-i+1,n-j+1)
            expected[row,column]+=v
        end
        q=QUBOComponent{T}(bits;linear,quadratic,offset=2)
        diagonal=[expected[i,i] for i in 1:n]
        indices=findall(!iszero,diagonal)
        @test q.linear.nzind==indices
        @test nonzeros(q.linear)==diagonal[indices]
        @test length(q.linear)==n && nnz(q.linear)==length(indices)
        @test all(iszero(q.quadratic[i,i]) for i in 1:n)
        @test dense_qubo(q).matrix==expected
        @test issorted(q.linear.nzind) && allunique(q.linear.nzind)
        patterns=unique([[f(i) for i in 1:n] for f in (i->false,i->true,isodd,i->mod(i,3)==0)])
        for original in patterns
            value=T(2)+sum(v*original[i] for (i,v) in linear;init=zero(T))+
                sum(v*original[i]*original[j] for (i,j,v) in quadratic;init=zero(T))
            @test energy(q,reverse(original))==value
        end
        before=component_artifact(q)
        isempty(linear) || (linear[1]=(1,coefficient(99)))
        isempty(quadratic) || (quadratic[1]=(1,1,coefficient(99)))
        @test component_artifact(q)==before
    end
    for T in (BigInt,Rational{BigInt})
        huge=T(big(2)^256)
        q=QUBOComponent{T}([BitID(:huge_diagonal,1)];linear=[(1,huge)],offset=-huge)
        @test energy(q,[false])==-huge && energy(q,[true])==0
        @test nonzeros(q.linear)==[huge]
    end
end
@testitem "Canonical bit sorting owns vector and view inputs" tags=[:q1] begin
    for T in (Float64,BigInt,Rational{BigInt}),n in (0,1,2,8,32),shape in (:vector,:view,:strided)
        coefficient(v)=T===Rational{BigInt} ? T(v//3) : T(v)
        original=reverse([BitID(Symbol(:owner,mod(i,3)),i;
            role=isodd(i) ? :primary : :semantic_auxiliary) for i in 1:n])
        storage=shape===:vector ? copy(original) :
            [BitID(:padding,i) for i in 1:2n+2]
        bits=if shape===:vector
            storage
        elseif shape===:view
            storage[2:n+1]=original
            view(storage,2:n+1)
        else
            storage[2:2:2n]=original
            view(storage,2:2:2n)
        end
        before=copy(storage)
        linear=[(i,coefficient(mod(i,7)-3)) for i in 1:n]
        quadratic=[(i,i,coefficient(1)) for i in 1:n]
        append!(quadratic,[(i,i+1,coefficient(3)) for i in 1:n-1])
        append!(quadratic,[(i+1,i,coefficient(-1)) for i in 1:n-1])
        q=QUBOComponent{T}(bits;linear,quadratic,offset=2)
        @test storage==before && bits==original
        @test q.bits==sort(original) && q.bits!==bits
        @test allunique(q.bits) && issorted(q.bits)
        positions=Dict(bit=>i for (i,bit) in enumerate(q.bits))
        for pattern in (i->false,i->true,isodd,i->mod(i,3)==0)
            input=[pattern(i) for i in 1:n]
            ordered=falses(n)
            for i in 1:n
                ordered[positions[original[i]]]=input[i]
            end
            value=T(2)+sum(v*input[i] for (i,v) in linear;init=zero(T))+
                sum(v*input[i]*input[j] for (i,j,v) in quadratic;init=zero(T))
            @test energy(q,ordered)==value
        end
        artifact=component_artifact(q)
        if n>0
            bits[1]=BitID(:changed_input,1)
            @test component_artifact(q)==artifact
            input_snapshot=copy(storage)
            q.bits[1]=BitID(:changed_result,1)
            @test storage==input_snapshot
        end
    end
    duplicate=[BitID(:same,1),BitID(:same,1)]
    original=copy(duplicate)
    @test_throws ArgumentError QUBOComponent(duplicate)
    @test duplicate==original
end
@testitem "Immutable bit copies preserve owned buffer relationships" tags=[:q1] begin
    @test all(T->T===Symbol || isbitstype(T),fieldtypes(BitID))
    for owner in (:x,Symbol("space/λ"),Symbol("")),role in (:primary,:semantic_auxiliary,:quadratization_auxiliary),index in (1,17,typemax(Int))
        bit=BitID(owner,index;role)
        @test deepcopy(bit)===bit
        @test deepcopy((bit,bit))===(bit,bit)
    end
    for n in (0,1,2,16,128)
        bits=[BitID(:buffer,i) for i in 1:n]
        copied=deepcopy(bits)
        @test copied==bits && copied!==bits
        graph=deepcopy((bits,bits))
        @test graph[1]===graph[2] && graph[1]!==bits
        if n>0
            copied[1]=BitID(:changed,1)
            @test bits[1]==BitID(:buffer,1)
        end
    end
    bits=[BitID(:alias,i) for i in 1:8]
    graph=(bits,bits,reshape(bits,2,4),view(bits,2:5))
    copied=deepcopy(graph)
    @test copied[1]===copied[2] && copied[1]!==bits
    @test copied[3]==reshape(bits,2,4) && copied[4]==view(bits,2:5)
    copied[1][2]=BitID(:changed,2)
    @test copied[3][2,1]==copied[1][2] && copied[4][1]==copied[1][2]
    @test bits[2]==BitID(:alias,2)
    # Unassigned slots and shared direct memory remain valid after copying.
    memory=Memory{BitID}(undef,4);memory[1]=BitID(:memory,1);memory[3]=BitID(:memory,3)
    duplicated=deepcopy((memory,memory))
    @test duplicated[1]===duplicated[2] && duplicated[1]!==memory
    @test [isassigned(duplicated[1],i) for i in 1:4]==[isassigned(memory,i) for i in 1:4]
    @test duplicated[1][1]==memory[1] && duplicated[1][3]==memory[3]
    duplicated[1][1]=BitID(:changed_memory,1)
    @test memory[1]==BitID(:memory,1)
    # Default copying of mutable/non-BitID semantic values is still recursive.
    values=[[1],[2]];book=codebook(:nested_values,values)
    owned=deepcopy(book)
    @test owned.bits==book.bits && owned.bits!==book.bits
    @test owned.values==book.values && owned.values!==book.values && owned.values[1]!==book.values[1]
    @test owned.codes==book.codes && owned.codes!==book.codes
    owned.values[1][1]=99
    @test book.values[1]==[1]
    owned.bits[1]=BitID(:changed_book,1)
    @test book.bits[1]!=owned.bits[1]
    cycle=Any[bits];push!(cycle,cycle)
    owned_cycle=deepcopy(cycle)
    @test owned_cycle[2]===owned_cycle && owned_cycle[1]!==bits && owned_cycle[1]==bits
end

@testitem "Canonical identity guard preserves rejection and caller ownership" tags=[:q1] begin
    for T in (Float64,BigInt,Rational{BigInt}), n in (2,3,8,64),
            placement in (:start,:middle,:end), shape in (:vector,:strided)
        unique_bits=[BitID(Symbol(:guard_owner,mod(i,3)),i;
            role=isodd(i) ? :primary : :semantic_auxiliary) for i in 1:n]
        duplicate=placement===:start ? unique_bits[end] :
            placement===:middle ? unique_bits[1] : unique_bits[div(n,2)]
        input=reverse(unique_bits)
        position=placement===:start ? 1 : placement===:middle ? div(n,2)+1 : n+1
        insert!(input,position,duplicate)
        storage=shape===:vector ? copy(input) : [BitID(:guard_padding,i) for i in 1:2length(input)+2]
        bits=if shape===:vector
            storage
        else
            storage[2:2:2length(input)]=input
            view(storage,2:2:2length(input))
        end
        before=copy(storage)
        @test length(Set(bits))==length(bits)-1
        visited=Ref(false)
        linear=((visited[]=true; (0,1)) for _ in 1:1)
        @test_throws ArgumentError QUBOComponent{T}(bits;linear)
        @test !visited[] && storage==before
    end
    # Owner and role are part of the identity even when indices agree.
    identities=[BitID(:a,1),BitID(:b,1),
        BitID(:a,1;role=:semantic_auxiliary),
        BitID(:a,1;role=:quadratization_auxiliary)]
    for T in (Float64,BigInt,Rational{BigInt}), permutation in
            ([1,2,3,4],[4,3,2,1],[2,4,1,3])
        input=identities[permutation]
        snapshot=copy(input)
        q=QUBOComponent{T}(input;linear=[(i,i) for i in eachindex(input)])
        @test q.bits==sort(identities) && input==snapshot && q.bits!==input
        @test length(Set(q.bits))==length(input)
        for pattern in (i->false,i->true,isodd)
            assignment=[pattern(i) for i in eachindex(input)]
            ordered=[assignment[findfirst(==(bit),input)] for bit in q.bits]
            @test energy(q,ordered)==sum(i*assignment[i] for i in eachindex(input))
        end
    end
end

@testitem "Indexed codebook coverage preserves ownership and validation priority" tags=[:q1] begin
    function coverage_book(variable,bits)
        n=length(bits)
        n==0 && return Codebook(variable,bits,[0],falses(0,1),[1])
        return Codebook(variable,bits,[-1,1],hcat(falses(n),trues(n)),[1,2])
    end
    for T in (Float64,BigInt,Rational{BigInt}), n in (0,1,2,16,127,128,129,512)
        coefficient(v)=T===Rational{BigInt} ? T(v//3) : T(v)
        identities=[BitID(:covered,i) for i in 1:n]
        bits=reverse(identities)
        split=cld(n,2)
        books=AbstractCodebook[coverage_book(:left,identities[1:split]),
            coverage_book(:right,identities[split+1:n])]
        snapshot=copy(bits)
        linear=[(i,coefficient(mod(i,7)-3)) for i in 1:n]
        quadratic=[(i,i+1,coefficient(2)) for i in 1:n-1]
        q=QUBOComponent{T}(bits;linear,quadratic,offset=2,codebooks=books)
        @test q.bits==identities && q.bits!==bits && bits==snapshot
        @test all(zip(q.codebooks,books)) do (owned,input)
            owned.variable==input.variable && owned.bits==input.bits &&
                owned.bits!==input.bits && owned.values==input.values &&
                owned.values!==input.values && owned.codes==input.codes &&
                owned.codes!==input.codes
        end
        source_positions=Dict(bit=>i for (i,bit) in enumerate(bits))
        for pattern in (i->false,i->true,isodd,i->mod(i,3)==0)
            original=[pattern(i) for i in 1:n]
            ordered=[original[source_positions[bit]] for bit in q.bits]
            value=T(2)+sum(v*original[i] for (i,v) in linear;init=zero(T))+
                sum(v*original[i]*original[j] for (i,j,v) in quadratic;init=zero(T))
            @test energy(q,ordered)==value
        end
        if n>0
            artifact=component_artifact(q)
            books[1].bits[1]=BitID(:changed_coverage_input,1)
            @test component_artifact(q)==artifact
        end
        n>=2 || continue
        outside=BitID(:absent_coverage,1)
        fixtures=[
            ("overlapping codebooks",AbstractCodebook[
                coverage_book(:first,identities[1:split]),
                coverage_book(:second,identities[split:n])]),
            ("codebook bit missing from component",AbstractCodebook[
                coverage_book(:first,identities[1:n-1]),
                coverage_book(:second,[outside])]),
            ("overlapping codebooks",AbstractCodebook[
                coverage_book(:first,identities),
                coverage_book(:second,[identities[1],outside])]),
            ("overlapping codebooks",AbstractCodebook[
                coverage_book(:first,identities),
                coverage_book(:second,[outside]),coverage_book(:third,[outside])]),
            ("duplicate codebook variable",AbstractCodebook[
                coverage_book(:same,identities[1:split]),
                coverage_book(:same,identities[split+1:n])])]
        for (message,invalid_books) in fixtures
            original_books=deepcopy([book.bits for book in invalid_books])
            error=try
                QUBOComponent{T}(bits;codebooks=invalid_books)
                nothing
            catch exception
                exception
            end
            @test error isa ArgumentError
            @test sprint(showerror,error)=="ArgumentError: "*message
            @test bits==snapshot && [book.bits for book in invalid_books]==original_books
        end
    end
end

@testitem "Direct codebook variable checks preserve heterogeneous key semantics" tags=[:q1] begin
    struct VariableKeyFixture{V} <: AbstractCodebook
        variable::V
        bits::Vector{BitID}
        values::Vector{Int}
    end
    for kind in (:symbol,:integer,:string,:mixed), n in (0,1,16,63,64,65,128,1024)
        key(i)=kind===:symbol ? Symbol(:key,i) : kind===:integer ? i :
            kind===:string ? "key_"*string(i) : isodd(i) ? Float64(i) : i
        books=AbstractCodebook[VariableKeyFixture(key(i),[BitID(:variable_key,i)],[-1,1]) for i in 1:n]
        bits=reverse([BitID(:variable_key,i) for i in 1:n])
        snapshot=copy(bits)
        keys=[book.variable for book in books]
        @test length(Set(keys))==n
        linear=[(i,Float64(mod(i,5)-2)) for i in 1:n]
        q=QUBOComponent{Float64}(bits;linear,codebooks=books)
        @test q.bits==reverse(bits) && q.bits!==bits && bits==snapshot &&
            all(zip(q.codebooks,books)) do (owned,input)
                isequal(owned.variable,input.variable) &&
                    owned.bits==input.bits && owned.bits!==input.bits &&
                    owned.values==input.values && owned.values!==input.values
            end
        positions=Dict(bit=>i for (i,bit) in enumerate(bits))
        for pattern in (i->false,i->true,isodd,i->mod(i,3)==0)
            assignment=[pattern(i) for i in 1:n]
            ordered=[assignment[positions[bit]] for bit in q.bits]
            @test energy(q,ordered)==sum(v*assignment[i] for (i,v) in linear;init=0.0)
        end
        if n>0
            owned_values=copy(q.codebooks[1].values)
            books[1].values[1]=99
            @test q.codebooks[1].values==owned_values &&
                isequal([book.variable for book in books],keys)
        end
        n>=2 || continue
        invalid=copy(books)
        invalid[end]=VariableKeyFixture(first(keys),copy(books[end].bits),[-1,1])
        input_bits=deepcopy([book.bits for book in invalid])
        error=try
            QUBOComponent{Float64}(bits;codebooks=invalid)
            nothing
        catch exception
            exception
        end
        @test error isa ArgumentError
        @test sprint(showerror,error)=="ArgumentError: duplicate codebook variable"
        @test bits==snapshot && [book.bits for book in invalid]==input_bits
    end
    for (first_key,second_key,duplicate) in
            ((1,1.0,true),(big(1),1//1,true),(NaN,NaN,true),
             (-0.0,0.0,false),(:x,"x",false))
        keys=Any[first_key,second_key,collect(3:64)...]
        books=AbstractCodebook[VariableKeyFixture(key,[BitID(:special_key,i)],[-1,1])
            for (i,key) in enumerate(keys)]
        bits=[BitID(:special_key,i) for i in eachindex(keys)]
        @test (length(Set(keys))<length(keys))==duplicate
        result=try
            QUBOComponent{Float64}(bits;codebooks=books)
        catch exception
            exception
        end
        if duplicate
            @test result isa ArgumentError
            @test sprint(showerror,result)=="ArgumentError: duplicate codebook variable"
        else
            @test length(result.codebooks)==length(books)
            @test isequal([book.variable for book in result.codebooks],keys)
        end
    end
    # Duplicate variables retain priority even when coverage is also invalid.
    for missing in (false,true)
        books=AbstractCodebook[VariableKeyFixture(Symbol(:priority,i),[BitID(:priority,i)],[-1,1]) for i in 1:64]
        books[end]=VariableKeyFixture(:priority1,[missing ? BitID(:absent_priority,1) : BitID(:priority,1)],[-1,1])
        bits=[BitID(:priority,i) for i in 1:64]
        error=try
            QUBOComponent{Float64}(bits;codebooks=books)
            nothing
        catch exception
            exception
        end
        @test error isa ArgumentError
        @test sprint(showerror,error)=="ArgumentError: duplicate codebook variable"
    end
end

@testitem "Validated offsets preserve values ownership and error priority" tags=[:q1] begin
    for T in (Float64,BigInt,Rational{BigInt}),n in (0,1,2,16),
            raw_offset in (0,2,2.0,big(2),2//1,big(2)//big(1))
        bits=reverse([BitID(:offset_regression,i) for i in 1:n])
        linear=[(i,T(mod(i,5)-2)) for i in 1:n]
        quadratic=[(i,i+1,T(1)) for i in 1:n-1]
        before=deepcopy((bits,linear,quadratic,raw_offset))
        q=QUBOComponent{T}(bits;linear,quadratic,offset=raw_offset)
        @test q.offset==T(raw_offset)
        @test (bits,linear,quadratic,raw_offset)==before
        @test q.bits==sort(bits) && q.bits!==bits
        for pattern in (i->false,i->true,isodd)
            input=[pattern(i) for i in 1:n]
            expected=T(raw_offset)+sum(v*input[i] for (i,v) in linear;init=zero(T))+
                sum(v*input[i]*input[j] for (i,j,v) in quadratic;init=zero(T))
            @test energy(q,reverse(input))==expected
        end
        if raw_offset isa T
            @test q.offset===raw_offset
        end
    end
    for T in (Float64,BigInt,Rational{BigInt})
        empty=QUBOComponent{T}(BitID[])
        @test empty.offset==zero(T) && energy(empty,Bool[])==zero(T)
    end
    for T in (Float64,Rational{BigInt})
        @test QUBOComponent{T}(BitID[];offset=2.5).offset==T(5//2)
        for offset in (NaN,Inf,-Inf,1//0)
            @test_throws ArgumentError QUBOComponent{T}(BitID[];offset)
        end
    end
    for offset in (2.5,NaN,Inf,-Inf,1//0)
        @test_throws InexactError QUBOComponent{BigInt}(BitID[];offset)
    end
    huge=big(2)^4096
    for T in (BigInt,Rational{BigInt})
        q=QUBOComponent{T}([BitID(:huge_offset,1)];linear=[(1,3)],offset=huge)
        @test energy(q,[0])==huge && energy(q,[1])==huge+3
        @test huge==big(2)^4096
    end
    function offset_error(bits;kwargs...)
        try
            QUBOComponent{Float64}(bits;kwargs...)
            error("expected rejection")
        catch e
            sprint(showerror,e)
        end
    end
    bits=[BitID(:offset_guard,1)]
    @test offset_error(bits;linear=[(1,NaN)],offset="unsupported")=="ArgumentError: nonfinite coefficient"
    @test offset_error(bits;linear=[(0,NaN)],offset=Inf)=="ArgumentError: linear index out of bounds"
    @test offset_error([bits;bits];offset="unsupported")=="ArgumentError: duplicate bit identity"
    @test offset_error(bits;linear=[(1,1e308),(1,1e308)],offset=2)=="ArgumentError: coefficient overflow"
    @test offset_error(bits;linear=[(1,1e308),(1,1e308)],offset=Inf)=="ArgumentError: nonfinite coefficient"
end

@testitem "Canonical remaps preserve arbitrary input order and codebook validation" tags=[:q1] begin
    function remap_book(variable,bits)
        n=length(bits)
        n==0 && return Codebook(variable,bits,[0],falses(0,1),[1])
        Codebook(variable,bits,[-1,1],hcat(falses(n),trues(n)),[1,2])
    end
    for T in (Float64,BigInt,Rational{BigInt}),n in (0,1,2,31,32,33,128,129),
            order in (:canonical,:reversed),shape in (:vector,:strided),with_books in (false,true)
        identities=sort([BitID(Symbol(:remap_owner,mod(i,3)),i;
            role=!with_books&&mod(i,3)==0 ? :semantic_auxiliary : :primary) for i in 1:n])
        input=order===:canonical ? copy(identities) : reverse(identities)
        storage=shape===:vector ? copy(input) : [BitID(:remap_padding,i) for i in 1:2n+2]
        bits=if shape===:vector
            storage
        else
            storage[2:2:2n]=input
            view(storage,2:2:2n)
        end
        before=copy(storage)
        coefficient(v)=T===Rational{BigInt} ? T(v//3) : T(v)
        linear=[(i,coefficient(mod(i,7)-3)) for i in 1:n]
        quadratic=[(i,i+1,coefficient(2)) for i in 1:n-1]
        append!(quadratic,[(i+1,i,coefficient(-1)) for i in 1:n-1])
        append!(quadratic,[(i,i,coefficient(1)) for i in 1:n])
        split=cld(n,2)
        books=with_books ? AbstractCodebook[remap_book(:left,identities[1:split]),
            remap_book(:right,identities[split+1:n])] : AbstractCodebook[]
        meanings=Dict(bit=>"meaning-$i" for (i,bit) in enumerate(input) if bit.role!==:primary)
        q=QUBOComponent{T}(bits;linear,quadratic,offset=2,codebooks=books,
            auxiliary_meanings=meanings,applicability="remap fixture",provenance="original order")
        @test q.bits==identities && q.bits!==bits && storage==before
        @test q.auxiliary_meanings==meanings
        @test q.applicability=="remap fixture" && q.provenance=="original order"
        @test length(q.codebooks)==length(books)
        @test all(zip(q.codebooks,books)) do (owned,original)
            owned!==original && owned.variable==original.variable && owned.bits==original.bits &&
                owned.bits!==original.bits && owned.values==original.values && owned.values!==original.values
        end
        positions=Dict(bit=>i for (i,bit) in enumerate(input))
        for pattern in (i->false,i->true,isodd,i->mod(i,3)==0)
            original=[pattern(i) for i in 1:n]
            z=[original[positions[bit]] for bit in q.bits]
            expected=T(2)+sum(v*original[i] for (i,v) in linear;init=zero(T))+
                sum(v*original[i]*original[j] for (i,j,v) in quadratic;init=zero(T))
            @test energy(q,z)==expected
        end
        @test storage==before && all(zip(q.codebooks,books)) do (owned,original)
            owned.bits==original.bits && owned.values==original.values
        end
        if n>0
            bits[1]=BitID(:changed_input,1)
            @test q.bits==identities
            q.bits[1]=BitID(:changed_output,1)
            @test bits[1]==BitID(:changed_input,1)
        end
    end
    function constructor_error(bits;kwargs...)
        try
            QUBOComponent(bits;kwargs...)
            error("expected constructor rejection")
        catch e
            sprint(showerror,e)
        end
    end
    for n in (1,16,128),order in (:canonical,:reversed)
        identities=[BitID(:guard_remap,i) for i in 1:n]
        bits=order===:canonical ? identities : reverse(identities)
        overlap=AbstractCodebook[remap_book(:a,identities),remap_book(:b,identities)]
        @test constructor_error(bits;codebooks=overlap)=="ArgumentError: overlapping codebooks"
        absent=BitID(:absent_remap,1)
        missing=AbstractCodebook[remap_book(:a,[absent])]
        @test constructor_error(bits;codebooks=missing)=="ArgumentError: codebook bit missing from component"
        repeated_missing=AbstractCodebook[remap_book(:a,[absent]),remap_book(:b,[absent])]
        @test constructor_error(bits;codebooks=repeated_missing)=="ArgumentError: overlapping codebooks"
        duplicate_variables=AbstractCodebook[remap_book(:a,identities),remap_book(:a,[absent])]
        @test constructor_error(bits;codebooks=duplicate_variables)=="ArgumentError: duplicate codebook variable"
        @test constructor_error(bits;linear=[(0,1)],codebooks=duplicate_variables)=="ArgumentError: linear index out of bounds"
        visited=Ref(false)
        linear=((visited[]=true;(0,1)) for _ in 1:1)
        @test constructor_error([bits;bits[1]];linear)=="ArgumentError: duplicate bit identity"
        @test !visited[]
    end
end

@testitem "Owned exact energy accumulators preserve polynomials and retained capacity" tags=[:q1] begin
    using SparseArrays
    function ordinary_exact_value(q,z)
        value=q.offset
        for p in eachindex(nonzeros(q.linear))
            z[q.linear.nzind[p]]==1 && (value+=nonzeros(q.linear)[p])
        end
        for j in axes(q.quadratic,2),p in nzrange(q.quadratic,j)
            z[j]==1 && z[rowvals(q.quadratic)[p]]==1 &&
                (value+=nonzeros(q.quadratic)[p])
        end
        value
    end
    for n in (0,1,2,4,8,9,10,16,64), precision in (0,4096), shape in (:linear,:mixed)
        huge=big(2)^precision
        bits=reverse([BitID(:exact_energy,i;role=mod(i,3)==0 ? :semantic_auxiliary : :primary) for i in 1:n])
        linear=precision==0 ? [(i,BigInt(mod(i,7)-3)) for i in 1:n] :
            [(i,i==9 ? huge : i==10 ? -huge : BigInt(1)) for i in 1:n]
        quadratic=shape===:linear ? Tuple{Int,Int,BigInt}[] :
            [(i,i+1,BigInt(2)) for i in 1:n-1]
        if shape===:mixed
            append!(quadratic,[(i+1,i,BigInt(-1)) for i in 1:n-1])
            append!(quadratic,[(i,i,BigInt(1)) for i in 1:n])
        end
        q=QUBOComponent(bits;linear,quadratic,offset=2)
        snapshot=deepcopy((q.offset,nonzeros(q.linear),nonzeros(q.quadratic),q.bits))
        positions=Dict(bit=>i for (i,bit) in enumerate(bits))
        for pattern in (i->false,i->true,isodd,i->mod(i,3)==0)
            original=[pattern(i) for i in 1:n]
            ordered=[original[positions[bit]] for bit in q.bits]
            expected=BigInt(2)+sum(v*BigInt(original[i]) for (i,v) in linear;init=BigInt(0))+
                sum(v*BigInt(original[i])*BigInt(original[j]) for (i,j,v) in quadratic;init=BigInt(0))
            for kind in (:bool,:integer,:float,:big,:mixed,:view)
                input=kind===:bool ? copy(ordered) : kind===:integer ? Int.(ordered) :
                    kind===:float ? Float64.(ordered) : kind===:big ? BigInt.(ordered) :
                    kind===:mixed ? Any[isodd(i) ? Int(ordered[i]) : Float64(ordered[i]) for i in 1:n] :
                    view(UInt8.(ordered),1:n)
                input_snapshot=deepcopy(input)
                ordinary=ordinary_exact_value(q,input)
                result=energy(q,input)
                @test result isa BigInt && result==expected
                @test (q.offset,nonzeros(q.linear),nonzeros(q.quadratic),q.bits)==snapshot
                @test input==input_snapshot
                @test result.alloc<=ordinary.alloc && Base.summarysize(result)<=Base.summarysize(ordinary)
                active=any(p->input[q.linear.nzind[p]]==1,eachindex(nonzeros(q.linear))) ||
                    any(j->input[j]==1 && any(p->input[rowvals(q.quadratic)[p]]==1,nzrange(q.quadratic,j)),axes(q.quadratic,2))
                if active
                    @test result!==q.offset && all(v->result!==v,[nonzeros(q.linear);nonzeros(q.quadratic)])
                    Base.GMP.MPZ.add!(result,result,BigInt(1))
                    @test (q.offset,nonzeros(q.linear),nonzeros(q.quadratic),q.bits)==snapshot && input==input_snapshot
                else
                    @test result===q.offset
                end
            end
        end
        @test_throws DimensionMismatch energy(q,falses(n+1))
        if n>0
            for value in (0.5,NaN,Inf,-Inf)
                invalid=zeros(Float64,n);invalid[1]=value
                @test_throws ArgumentError energy(q,invalid)
            end
            @test (q.offset,nonzeros(q.linear),nonzeros(q.quadratic),q.bits)==snapshot
        end
    end
    # Explicit small final value after a large peak beyond the ordinary prefix.
    for precision in (256,4096,16384)
        huge=big(2)^precision
        q=QUBOComponent([BitID(:peak,i) for i in 1:16];
            linear=[(i,i==9 ? huge : i==10 ? -huge : BigInt(1)) for i in 1:16])
        result=energy(q,trues(16));ordinary=ordinary_exact_value(q,trues(16))
        @test result==14
        @test result.alloc<=ordinary.alloc
        @test Base.summarysize(result)<=Base.summarysize(ordinary)
        @test q.linear.nzval[9]==huge && q.linear.nzval[10]==-huge
    end
end
