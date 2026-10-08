import QUBOConstraints as QC
using Random
const VG = QC.ValuePairGuidance

"Revalidate every exact truth row of a prepared unary, binary or ternary plan."
function atomic_validation_case(parameters)
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
    plan = QC.atomic_plan(books, expression; semantic_id="perf/atomic-validation/$(arity)/$(width)")
    prepare = () -> plan
    operation = state -> QC.validate_atomic_plan(state)
    verify = (state, result) -> begin
        result === true || return false
        for node in state.nodes
            node.operator in (:variable, :constant) && continue
            expected_rows = prod(length(state.nodes[i].domain) for i in node.inputs)
            length(node.rows) == expected_rows || return false
            for row in node.rows
                args = [Int(state.nodes[i].domain[k]) for (i, k) in zip(node.inputs, row)]
                expected = node.operator === :eq ? Int(args[1] == args[2]) :
                    node.operator === :neg ? -args[1] : args[1] == 1 ? args[2] : args[3]
                node.domain[last(row)] == expected || return false
            end
        end
        true
    end
    (; prepare, operation, verify)
end

"Canonical sparse construction with duplicate, reversed and diagonal input terms."
function component_construction_case(parameters)
    n = get(parameters, "bits", 128)
    n > 1 || throw(ArgumentError("at least two fixture bits required"))
    T = get(parameters, "coefficient_type", "float") == "exact" ? BigInt : Float64
    bits = reverse([QC.BitID(:x, i) for i in 1:n])
    linear = [(i, T(mod(i, 7) - 3)) for i in 1:n]
    quadratic = [(i, i + 1, T(2)) for i in 1:n-1]
    append!(quadratic, [(i + 1, i, T(-1)) for i in 1:n-1])
    append!(quadratic, [(i, i, T(1)) for i in 1:n])
    prepare = () -> (; bits, linear, quadratic)
    operation = if T === BigInt
        state -> QC.QUBOComponent{BigInt}(state.bits; linear=state.linear, quadratic=state.quadratic, offset=2)
    else
        state -> QC.QUBOComponent{Float64}(state.bits; linear=state.linear, quadratic=state.quadratic, offset=2)
    end
    verify = (state, q) -> begin
        for pattern in (i -> false, i -> true, isodd, i -> mod(i, 3) == 0)
            original = [pattern(i) for i in 1:n]
            z = reverse(original)
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
    bits = reverse(QC.BitID[b for book in books for b in book.bits])
    n = length(bits)
    linear = [(i, T(mod(i, 7) - 3)) for i in 1:n]
    quadratic = [(i, i + 1, T(2)) for i in 1:n-1]
    append!(quadratic, [(i + 1, i, T(-1)) for i in 1:n-1])
    append!(quadratic, [(i, i, T(1)) for i in 1:n])
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
            z = [original[findfirst(==(bit), state.bits)] for bit in q.bits]
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
