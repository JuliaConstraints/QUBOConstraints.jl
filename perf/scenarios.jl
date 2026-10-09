import QUBOConstraints as QC
using Random
using SparseArrays
const VG = QC.ValuePairGuidance

function _fixture_bit_order(identities, parameters)
    order = get(parameters, "input_order", "reversed")
    order == "reversed" && return reverse(identities)
    order == "canonical" && return sort!(copy(identities))
    order == "shuffled" && return shuffle(MersenneTwister(41), identities)
    throw(ArgumentError("unknown fixture bit order"))
end

"Deep-copy bit storage while preserving repeated-array, reshape and view relationships."
function bit_storage_copy_case(parameters)
    n=get(parameters,"bits",128)
    n>=4 || throw(ArgumentError("at least four alias fixture bits required"))
    bits=[QC.BitID(Symbol(:owner,mod(i,7)),i;role=isodd(i) ? :primary : :semantic_auxiliary) for i in 1:n]
    input=(bits,bits,reshape(bits,1,n),view(bits,2:n-1))
    prepare=()->input
    operation=deepcopy
    verify=(state,owned)->begin
        owned==state && owned[1]!==state[1] && owned[1]===owned[2] || return false
        expected=owned[1][2]
        changed=QC.BitID(:changed,2)
        owned[1][2]=changed
        valid=owned[3][1,2]==changed && owned[4][1]==changed && state[1][2]==expected
        owned[1][2]=expected
        valid && owned==state
    end
    (;prepare,operation,verify)
end

function rational_component_case(parameters)
    n=get(parameters,"bits",64)
    n>=0 || throw(ArgumentError("nonnegative rational fixture width required"))
    dense=get(parameters,"dense",false)
    T=Rational{BigInt}
    bits=_fixture_bit_order([QC.BitID(:rational,i) for i in 1:n],parameters)
    linear=[(i,T((mod(i,7)-3)//3)) for i in 1:n]
    quadratic=dense ? [(i,j,T((mod(i+j,7)-3)//5)) for j in 1:n for i in 1:j] :
        [(i,i+1,T(2//5)) for i in 1:n-1]
    append!(quadratic,[(i,i,T(-1//3)) for i in 1:n])
    positions=Dict(bit=>i for (i,bit) in enumerate(bits))
    prepare=()->(;bits,linear,quadratic)
    operation=s->QC.QUBOComponent{Rational{BigInt}}(s.bits;linear=s.linear,quadratic=s.quadratic,offset=2)
    verify=(s,q)->begin
        for pattern in (i->false,i->true,isodd,i->mod(i,3)==0)
            original=[pattern(i) for i in 1:n]
            expected=T(2)+sum(v*original[i] for (i,v) in s.linear;init=zero(T))+
                sum(v*original[i]*original[j] for (i,j,v) in s.quadratic;init=zero(T))
            QC.energy(q,[original[positions[bit]] for bit in q.bits])==expected || return false
        end
        true
    end
    (;prepare,operation,verify)
end


"Complete composition with shared primary books, renamed auxiliaries and promoted coefficients."
function component_composition_case(parameters)
    width=get(parameters,"width",64)
    width>=0 || throw(ArgumentError("nonnegative composition width required"))
    shape=get(parameters,"shape","mixed_float")
    T,S,with_books=shape=="bare_float" ? (Float64,Float64,false) :
        shape=="bare_exact" ? (BigInt,BigInt,false) :
        shape=="mixed_float" ? (Float64,BigInt,true) :
        shape=="mixed_rational" ? (BigInt,Rational{BigInt},true) :
        throw(ArgumentError("unknown composition shape"))
    books=with_books ? [QC.codebook(Symbol(:x,i),0:1;encoding=:domain_wall) for i in 1:width] :
        QC.AbstractCodebook[]
    primary=with_books ? QC.BitID[b for book in books for b in book.bits] :
        [QC.BitID(:x,i) for i in 1:width]
    bits=vcat(primary,[QC.BitID(:aux,i;role=:quadratization_auxiliary) for i in 1:width])
    right_bits=reverse(bits)
    left_linear=[(i,T(mod(i,5)-2)) for i in eachindex(bits)]
    right_linear=[(i,S===Rational{BigInt} ? S((mod(i,7)-3)//3) : S(mod(i,7)-3)) for i in eachindex(bits)]
    left_quadratic=[(i,i+1,T(2)) for i in 1:length(bits)-1]
    right_quadratic=[(i,i+1,S(-1)) for i in 1:length(bits)-1]
    meanings=Dict(bit=>"original auxiliary" for bit in bits if bit.role!==:primary)
    a=QC.QUBOComponent{T}(bits;linear=left_linear,quadratic=left_quadratic,offset=1,
        codebooks=books,auxiliary_meanings=meanings,applicability="left",provenance="left")
    b=QC.QUBOComponent{S}(right_bits;linear=right_linear,quadratic=right_quadratic,offset=-2,
        codebooks=books,auxiliary_meanings=meanings,applicability="right",provenance="right")
    rename(bit)=bit.role===:primary ? bit :
        QC.BitID(Symbol(:right,"/",bit.owner),bit.index;role=bit.role)
    expected_bits=sort!(unique(vcat(bits,rename.(right_bits))))
    R=T===Float64 || S===Float64 ? Float64 : promote_type(T,S)
    state=(;a,b,bits,right_bits,left_linear,right_linear,left_quadratic,right_quadratic,expected_bits)
    prepare=()->state
    operation=s->QC.compose(s.a,s.b)
    verify=(s,q)->begin
        q isa QC.QUBOComponent{R} && q.bits==s.expected_bits &&
            q.applicability=="(left) AND (right)" && q.provenance=="add(left, right)" || return false
        length(q.codebooks)==length(s.a.codebooks) || return false
        all(zip(s.a.codebooks,q.codebooks)) do pair
            input,owned=pair
            input!==owned && input.values==owned.values && input.values!==owned.values &&
                input.codes==owned.codes && input.codes!==owned.codes
        end || return false
        q.auxiliary_meanings==merge(s.a.auxiliary_meanings,
            Dict(rename(bit)=>text for (bit,text) in s.b.auxiliary_meanings)) || return false
        positions=Dict(bit=>i for (i,bit) in enumerate(q.bits))
        for pattern in (i->false,i->true,isodd,i->mod(i,3)==0)
            z=[pattern(i) for i in eachindex(q.bits)]
            left=[z[positions[bit]] for bit in s.bits]
            right=[z[positions[rename(bit)]] for bit in s.right_bits]
            expected=s.a.offset+sum(v*left[i] for (i,v) in s.left_linear;init=zero(s.a.offset))+
                sum(v*left[i]*left[j] for (i,j,v) in s.left_quadratic;init=zero(s.a.offset))+
                s.b.offset+sum(v*right[i] for (i,v) in s.right_linear;init=zero(s.b.offset))+
                sum(v*right[i]*right[j] for (i,j,v) in s.right_quadratic;init=zero(s.b.offset))
            QC.energy(q,z)==expected || return false
        end
        true
    end
    (;prepare,operation,verify)
end


function atomic_fixture(parameters)
    width = get(parameters, "width", 16)
    arity = get(parameters, "arity", 2)
    width > 1 && arity in (1, 2, 3) || throw(ArgumentError("invalid atomic fixture shape"))
    N = QC.IntensionNode
    books = arity == 1 ? [QC.structured_codebook(:x, 0:width-1)] :
        arity == 2 ? [QC.structured_codebook(s, 0:width-1) for s in (:x, :y)] :
        [QC.structured_codebook(:c, 0:1), QC.structured_codebook(:x, 0:width-1),
            QC.structured_codebook(:y, 0:width-1)]
    expression = arity == 1 ? N(:eq, N(:neg, :x), 0) : arity == 2 ? N(:eq, :x, :y) :
        N(:eq, N(Symbol("if"), :c, :x, :y), 0)
    (; books, expression, semantic_id="perf/atomic-fixture/$(arity)/$(width)")
end

function atomic_fixture_oracle(plan)
    for node in plan.nodes
        node.operator in (:variable, :constant) && continue
        expected_rows = prod(length(plan.nodes[i].domain) for i in node.inputs)
        length(node.rows) == expected_rows || return false
        for row in node.rows
            args = [Int(plan.nodes[i].domain[k]) for (i, k) in zip(node.inputs, row)]
            expected = node.operator === :eq ? Int(args[1] == args[2]) :
                node.operator === :neg ? -args[1] : args[1] == 1 ? args[2] : args[3]
            node.domain[last(row)] == expected || return false
        end
    end
    true
end

"Revalidate every exact truth row of a prepared unary, binary or ternary plan."
function atomic_validation_case(parameters)
    fixture = atomic_fixture(parameters)
    plan = QC.atomic_plan(fixture.books, fixture.expression; semantic_id=fixture.semantic_id)
    prepare = () -> plan
    operation = state -> QC.validate_atomic_plan(state)
    verify = (state, result) -> result === true && atomic_fixture_oracle(state)
    (; prepare, operation, verify)
end

"Construct the complete ordered plan from already prepared source books and expression."
function atomic_construction_case(parameters)
    fixture = atomic_fixture(parameters)
    prepare = () -> fixture
    operation = state -> QC.atomic_plan(state.books, state.expression; semantic_id=state.semantic_id)
    verify = (state, plan) -> QC.validate_atomic_plan(plan) && atomic_fixture_oracle(plan) &&
        all(zip(state.books, plan.codebooks)) do (input, owned)
            input !== owned && input.values == owned.values && input.values !== owned.values &&
                input.bits == owned.bits && input.bits !== owned.bits
        end
    (; prepare, operation, verify)
end

"Canonical sparse construction with duplicate, reversed and diagonal input terms."
function component_construction_case(parameters)
    n = get(parameters, "bits", 128)
    n > 1 || throw(ArgumentError("at least two fixture bits required"))
    T = get(parameters, "coefficient_type", "float") == "exact" ? BigInt : Float64
    bits = _fixture_bit_order([QC.BitID(:x, i) for i in 1:n], parameters)
    linear = [(i, T(mod(i, 7) - 3)) for i in 1:n]
    quadratic = [(i, i + 1, T(2)) for i in 1:n-1]
    append!(quadratic, [(i + 1, i, T(-1)) for i in 1:n-1])
    append!(quadratic, [(i, i, T(1)) for i in 1:n])
    positions = Dict(bit => i for (i, bit) in enumerate(bits))
    prepare = () -> (; bits, linear, quadratic)
    operation = if T === BigInt
        state -> QC.QUBOComponent{BigInt}(state.bits; linear=state.linear, quadratic=state.quadratic, offset=2)
    else
        state -> QC.QUBOComponent{Float64}(state.bits; linear=state.linear, quadratic=state.quadratic, offset=2)
    end
    verify = (state, q) -> begin
        for pattern in (i -> false, i -> true, isodd, i -> mod(i, 3) == 0)
            original = [pattern(i) for i in 1:n]
            z = [original[positions[bit]] for bit in q.bits]
            expected = T(2) + sum(v * original[i] for (i, v) in state.linear) +
                sum(v * original[i] * original[j] for (i, j, v) in state.quadratic)
            QC.energy(q, z) == expected || return false
        end
        true
    end
    (; prepare, operation, verify)
end

"Construction with owned heterogeneous codebooks covering every supplied bit."
function component_codebook_case(parameters)
    count = get(parameters, "codebooks", 64)
    count > 0 || throw(ArgumentError("at least one fixture codebook required"))
    T = get(parameters, "coefficient_type", "float") == "exact" ? BigInt : Float64
    books = QC.AbstractCodebook[isodd(i) ? QC.codebook(Symbol(:x, i), [-3, 5]) :
        QC.structured_codebook(Symbol(:x, i), [-3, 5]; encoding=:one_hot) for i in 1:count]
    bits = _fixture_bit_order(QC.BitID[b for book in books for b in book.bits], parameters)
    n = length(bits)
    linear = [(i, T(mod(i, 7) - 3)) for i in 1:n]
    quadratic = [(i, i + 1, T(2)) for i in 1:n-1]
    append!(quadratic, [(i + 1, i, T(-1)) for i in 1:n-1])
    append!(quadratic, [(i, i, T(1)) for i in 1:n])
    positions = Dict(bit => i for (i, bit) in enumerate(bits))
    prepare = () -> (; bits, linear, quadratic, books)
    operation = if T === BigInt
        state -> QC.QUBOComponent{BigInt}(state.bits; linear=state.linear,
            quadratic=state.quadratic, offset=2, codebooks=state.books)
    else
        state -> QC.QUBOComponent{Float64}(state.bits; linear=state.linear,
            quadratic=state.quadratic, offset=2, codebooks=state.books)
    end
    verify = (state, q) -> begin
        for pattern in (i -> false, i -> true, isodd, i -> mod(i, 3) == 0)
            original = [pattern(i) for i in 1:n]
            z = [original[positions[bit]] for bit in q.bits]
            expected = T(2) + sum(v * original[i] for (i, v) in state.linear) +
                sum(v * original[i] * original[j] for (i, j, v) in state.quadratic)
            QC.energy(q, z) == expected || return false
        end
        length(q.codebooks) == count || return false
        all(zip(q.codebooks, state.books)) do (owned, input)
            owned !== input && owned.variable == input.variable &&
                owned.bits == input.bits && owned.bits !== input.bits &&
                owned.values == input.values && owned.values !== input.values
        end
    end
    (; prepare, operation, verify)
end

# Retained whole-interaction reference for matched-work comparisons. It has the
# same canonical polynomial and owned buffers as the sparse operation below.
function scan_delta_reference!(w, g, ids, replacements)
    fill!(w.change, 0.0)
    for (i, v) in zip(ids, replacements), k in g.by_variable[i]
        w.change[k] = Float64(g.atoms[k][2] == v) - w.active[k]
    end
    total = 0.0
    @inbounds @simd for i in eachindex(g.linear)
        total += g.linear[i] * w.change[i]
    end
    @inbounds @simd for k in eachindex(g.coefficient)
        i, j = g.left[k], g.right[k]
        total += g.coefficient[k] * (w.change[i] * w.active[j] +
            w.active[i] * w.change[j] + w.change[i] * w.change[j])
    end
    total
end

"One-hot ring guidance; sparse and reference delta perform the same fixed number of moves."
function value_pair_case(parameters)
    n = get(parameters, "variables", 128)
    depth = get(parameters, "depth", 4)
    repetitions = get(parameters, "repetitions", 1024)
    operation_name = get(parameters, "operation", "delta")
    mode = get(parameters, "mode", "conditional")
    moved_variables = get(parameters, "moved_variables", 2)
    1 <= moved_variables <= n || throw(ArgumentError("invalid fixture move width"))
    atoms = [(i, v) for i in 1:n for v in 1:8]
    pairs = Dict{Tuple{Int,Int},Float64}()
    for i in 1:n, a in 1:8, b in 1:8
        j = mod1(i + 1, n)
        i == j && continue
        key = minmax(8(i - 1) + a, 8(j - 1) + b)
        pairs[key] = get(pairs, key, 0.0) + (a == b ? 1.0 : -0.25)
    end
    g = VG.Guide(atoms, zeros(length(atoms)), [(a, b, c) for ((a, b), c) in sort!(collect(pairs); by=first)])
    prepare = () -> begin
        w = VG.Workspace(g)
        values = [mod1(i, 8) for i in 1:n]
        VG.refresh!(w, g, values)
        ids = collect(1:moved_variables)
        replacements = [mod1(values[i] + 2, 8) for i in ids]
        (; g, w, values, rng=Xoshiro(41), ids, replacements)
    end
    operation = if operation_name == "delta"
        state -> begin
            result = 0.0
            for _ in 1:repetitions
                result = VG.delta!(state.w, state.g, state.ids, state.replacements)
            end
            result
        end
    elseif operation_name == "scan_delta"
        state -> begin
            result = 0.0
            for _ in 1:repetitions
                result = scan_delta_reference!(state.w, state.g, state.ids, state.replacements)
            end
            result
        end
    elseif operation_name == "scope"
        state -> begin
            result = 0
            for _ in 1:repetitions
                result = length(VG.scope!(state.w, state.g, state.values, depth, state.rng; mode))
            end
            result
        end
    elseif operation_name == "proposal"
        state -> begin
            result = 0.0
            for _ in 1:repetitions
                result = VG.proposal!(state.w, state.g, state.values, depth, state.rng; mode).delta
            end
            result
        end
    else
        throw(ArgumentError("unknown guidance workload"))
    end
    polynomial(values) = begin
        z = [Float64(values[i] == v) for (i, v) in atoms]
        sum(c * z[a] * z[b] for ((a, b), c) in pairs)
    end
    verify = (state, result) -> begin
        if operation_name in ("delta", "scan_delta")
            candidate = copy(state.values)
            candidate[collect(state.ids)] = collect(state.replacements)
            return result == polynomial(candidate) - polynomial(state.values)
        elseif operation_name == "scope"
            return result == min(depth, n) && allunique(state.w.chosen)
        else
            candidate = copy(state.values)
            candidate[state.w.chosen] = state.w.replacements
            return result <= 0.0 && result == polynomial(candidate) - polynomial(state.values)
        end
    end
    (; prepare, operation, verify)
end

"Numeric binary energy traversal; the direct polynomial is an independent oracle."
function numeric_energy_case(parameters)
    n = get(parameters, "bits", 1024)
    repetitions = get(parameters, "repetitions", 1024)
    kind = get(parameters, "input_type", "integer")
    R = kind == "float" ? Float64 : kind == "boolean" ? Bool : Int
    q = QC.QUBOComponent{Float64}([QC.BitID(:x, i) for i in 1:n];
        linear=[(i, Float64(mod(i, 7) - 3)) for i in 1:n],
        quadratic=[(i, i + 1, 0.25) for i in 1:n-1])
    prepare = () -> (; q, z=R.(isodd.(1:n)))
    operation = state -> begin
        result = 0.0
        for _ in 1:repetitions
            result = QC.energy(state.q, state.z)
        end
        result
    end
    verify = (state, result) -> result == sum(Float64(mod(i, 7) - 3) * state.z[i] for i in 1:n) +
        sum(0.25 * state.z[i] * state.z[i + 1] for i in 1:n-1)
    (; prepare, operation, verify)
end

"Construction controls for validated offset conversion, excluding prepared inputs."
function component_offset_case(parameters)
    n = get(parameters, "bits", 0)
    kind = get(parameters, "coefficient_type", "exact")
    offset_kind = get(parameters, "offset_kind", "integer")
    n >= 0 || throw(ArgumentError("nonnegative offset fixture width required"))
    kind in ("float", "exact", "rational") || throw(ArgumentError("unknown offset coefficient type"))
    offset_kind in ("integer", "integral_float", "already_typed") ||
        throw(ArgumentError("unknown offset input type"))
    T = kind == "float" ? Float64 : kind == "exact" ? BigInt : Rational{BigInt}
    bits = [QC.BitID(:offset_fixture, i) for i in 1:n]
    linear = [(i, T(mod(i, 7) - 3)) for i in 1:n]
    quadratic = [(i, i + 1, T(1)) for i in 1:n-1]
    offset = offset_kind == "integer" ? 2 : offset_kind == "integral_float" ? 2.0 : T(2)
    before = deepcopy((bits, linear, quadratic, offset))
    prepare = () -> (; bits, linear, quadratic, offset)
    operation = if T === Float64
        s -> QC.QUBOComponent{Float64}(s.bits; linear=s.linear, quadratic=s.quadratic, offset=s.offset)
    elseif T === BigInt
        s -> QC.QUBOComponent{BigInt}(s.bits; linear=s.linear, quadratic=s.quadratic, offset=s.offset)
    else
        s -> QC.QUBOComponent{Rational{BigInt}}(s.bits; linear=s.linear, quadratic=s.quadratic, offset=s.offset)
    end
    verify = (s, q) -> begin
        (s.bits, s.linear, s.quadratic, s.offset) == before || return false
        q.bits == s.bits && q.bits !== s.bits && q.offset == T(s.offset) || return false
        for pattern in (i -> false, i -> true, isodd)
            z = [pattern(i) for i in 1:n]
            expected = T(s.offset) + sum(v*z[i] for (i,v) in s.linear; init=zero(T)) +
                sum(v*z[i]*z[j] for (i,j,v) in s.quadratic; init=zero(T))
            QC.energy(q, z) == expected || return false
        end
        offset_kind == "already_typed" && q.offset !== s.offset && return false
        true
    end
    (; prepare, operation, verify)
end

"Exact scalar energy, with independent input polynomials and GMP limb-capacity controls."
function exact_energy_case(parameters)
    n = get(parameters, "bits", 128)
    precision = get(parameters, "precision_bits", 0)
    growth_precision = get(parameters, "growth_precision_bits", 0)
    kind = get(parameters, "input_type", "boolean")
    pattern = get(parameters, "pattern", "one")
    n >= 0 && precision >= 0 && growth_precision >= 0 || throw(ArgumentError("nonnegative fixture sizes required"))
    kind in ("boolean", "integer") || throw(ArgumentError("unknown binary input type"))
    pattern in ("zero", "one", "odd") || throw(ArgumentError("unknown binary pattern"))
    huge = big(2)^precision
    second_huge = growth_precision > 0 ? big(2)^growth_precision : -huge
    linear = precision == 0 ? [(i, BigInt(mod(i, 7) - 3)) for i in 1:n] :
        [(i, i == 9 ? huge : i == 10 ? second_huge : BigInt(1)) for i in 1:n]
    quadratic = precision == 0 ? [(i, i + 1, BigInt(1)) for i in 1:n-1] :
        Tuple{Int,Int,BigInt}[]
    q = QC.QUBOComponent([QC.BitID(:exact_energy_fixture, i) for i in 1:n];
        linear, quadratic, offset=2)
    flags = pattern == "zero" ? falses(n) : pattern == "one" ? trues(n) : isodd.(1:n)
    z = kind == "boolean" ? flags : Int.(flags)
    before = deepcopy((q.offset, nonzeros(q.linear), nonzeros(q.quadratic)))
    z_before = copy(z)
    prepare = () -> (; q, z)
    operation = state -> QC.energy(state.q, state.z)
    verify = (state, result) -> begin
        expected = BigInt(2) + sum(v * BigInt(state.z[i]) for (i, v) in linear; init=BigInt(0)) +
            sum(v * BigInt(state.z[i]) * BigInt(state.z[j]) for (i, j, v) in quadratic; init=BigInt(0))
        result == expected || return false
        (state.q.offset, nonzeros(state.q.linear), nonzeros(state.q.quadratic)) == before || return false
        state.z == z_before || return false
        # Ordinary non-mutating arithmetic supplies the retention baseline;
        # its order follows canonical sparse storage, while value checks above
        # evaluate the independent original input polynomial.
        ordinary = state.q.offset
        for p in eachindex(nonzeros(state.q.linear))
            state.z[state.q.linear.nzind[p]] == 1 && (ordinary += nonzeros(state.q.linear)[p])
        end
        for j in axes(state.q.quadratic, 2), p in nzrange(state.q.quadratic, j)
            state.z[j] == 1 && state.z[rowvals(state.q.quadratic)[p]] == 1 &&
                (ordinary += nonzeros(state.q.quadratic)[p])
        end
        result.alloc <= ordinary.alloc && Base.summarysize(result) <= Base.summarysize(ordinary) || return false
        pattern == "zero" && result !== state.q.offset && return false
        true
    end
    (; prepare, operation, verify)
end

"Complete construction for dense and long repeated, cancelling or zero raw coordinate lists."
function raw_coordinate_component_case(parameters)
    n=get(parameters,"bits",2)
    n>0 || throw(ArgumentError("positive raw coordinate fixture width required"))
    shape=get(parameters,"shape","repeated")
    shape in ("dense","repeated","cancelled","zero") || throw(ArgumentError("unknown raw coordinate shape"))
    kind=get(parameters,"coefficient_type","float")
    T=kind=="float" ? Float64 : kind=="exact" ? BigInt : kind=="rational" ? Rational{BigInt} :
        throw(ArgumentError("unknown raw coordinate coefficient type"))
    coefficient(v)=T===Rational{BigInt} ? T(v//5) : T(v)
    bits=_fixture_bit_order([QC.BitID(:raw_coordinate,i) for i in 1:n],parameters)
    count=get(parameters,"raw_terms",512)
    count>=0 || throw(ArgumentError("nonnegative raw term count required"))
    linear=shape=="dense" ? [(i,coefficient(mod(i,7)-3)) for i in 1:n] :
        [(1,coefficient(shape=="zero" ? 0 : shape=="cancelled" && iseven(i) ? -1 : 1)) for i in 1:count]
    quadratic=shape=="dense" ? [(i,j,coefficient(mod(i+j,7)-3)) for j in 1:n for i in 1:j] :
        [(1,1,coefficient(shape=="zero" ? 0 : shape=="cancelled" && iseven(i) ? -1 : 1)) for i in 1:count]
    positions=Dict(bit=>i for (i,bit) in enumerate(bits))
    snapshot=deepcopy((bits,linear,quadratic))
    state=(;bits,linear,quadratic)
    prepare=()->state
    operation=kind=="float" ? s->QC.QUBOComponent{Float64}(s.bits;linear=s.linear,quadratic=s.quadratic,offset=2) :
        kind=="exact" ? s->QC.QUBOComponent{BigInt}(s.bits;linear=s.linear,quadratic=s.quadratic,offset=2) :
        s->QC.QUBOComponent{Rational{BigInt}}(s.bits;linear=s.linear,quadratic=s.quadratic,offset=2)
    verify=(s,q)->begin
        (s.bits,s.linear,s.quadratic)==snapshot && q.bits!==s.bits && q.bits==sort(s.bits) || return false
        for pattern in (i->false,i->true,isodd,i->mod(i,3)==0)
            original=[pattern(i) for i in 1:n]
            expected=T(2)+sum(v*original[i] for (i,v) in s.linear;init=zero(T))+
                sum(v*original[i]*original[j] for (i,j,v) in s.quadratic;init=zero(T))
            QC.energy(q,[original[positions[bit]] for bit in q.bits])==expected || return false
        end
        true
    end
    (;prepare,operation,verify)
end
