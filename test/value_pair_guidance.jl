@testitem "Sparse value-pair simultaneous delta and ownership" tags=[:guidance] begin
    using Random
    const VG = QUBOConstraints.ValuePairGuidance
    atoms = [(3, 11), (1, -7), (3, -2), (1, 5), (2, 9), (2, -4)]
    linear = [0.5, -1.0, 0.25, 2.0, -0.5, 0.0]
    terms = [(i, j, Float64(mod(3i + j, 11) - 5) / 4) for i in 1:6 for j in i:6]
    g = VG.Guide(atoms, linear, terms)
    w, other = VG.Workspace(g), VG.Workspace(g)
    # Direct input polynomial includes diagonal terms and refers to the original
    # non-sorted semantic bindings, independently of all adjacency structures.
    polynomial(values) = begin
        z = [Float64(values[i] == v) for (i, v) in atoms]
        sum(linear .* z) + sum(c * z[i] * z[j] for (i, j, c) in terms)
    end
    domains = ((-7, 5), (-4, 9), (-2, 11))
    for assignment in Iterators.product(domains...), replacement in Iterators.product(domains...)
        values, candidate = collect(assignment), collect(replacement)
        VG.refresh!(w, g, values)
        original_active = copy(w.active)
        @test VG.energy(g, w) == polynomial(values)
        for ids in ((1,), (2, 1), (3, 1), (1, 2, 3), (3, 2, 1))
            moved = copy(values)
            changes = [candidate[i] for i in ids]
            for (i, v) in zip(ids, changes)
                moved[i] = v
            end
            @test VG.delta!(w, g, ids, changes) == polynomial(moved) - polynomial(values)
            @test allunique(w.touched_terms)
            @test w.active == original_active
        end
        @test VG.delta!(w, g, (), ()) == 0.0
        @test isempty(w.changed_atoms) && isempty(w.touched_terms)
    end
    VG.refresh!(w, g, [-7, -4, -2])
    expected = VG.delta!(w, g, (3, 1), (11, 5))
    w.epoch = typemax(UInt)
    @test VG.delta!(w, g, (1, 3), (5, 11)) == expected
    @test w.epoch == one(UInt)
    allocated_delta(w, g) = @allocated VG.delta!(w, g, (1, 3), (5, 11))
    allocated_delta(w, g)
    @test allocated_delta(w, g) == 0
    @test_throws ArgumentError VG.delta!(w, g, (1, 1), (5, 5))
    @test_throws ArgumentError VG.delta!(w, g, (1,), ())
    @test_throws BoundsError VG.delta!(w, g, (4,), (0,))
    @test_throws DimensionMismatch VG.refresh!(w, g, [1])
    @test_throws ArgumentError VG.energy(VG.Guide(atoms, linear, terms), w)
    @test w.change !== other.change && w.term_epoch !== other.term_epoch
    @test w.chosen !== other.chosen && w.replacements !== other.replacements
    before = VG.energy(g, w)
    owned_atoms, owned_linear, owned_terms = copy(atoms), copy(linear), copy(terms)
    owned = VG.Guide(owned_atoms, owned_linear, owned_terms)
    owned_atoms[1] = (99, 99); owned_linear[1] = 99.0; empty!(owned_terms)
    @test owned.atoms == atoms && owned.linear == g.linear && owned.coefficient == g.coefficient
    @test VG.energy(g, w) == before
    @test VG.metadata(g)["authority"] == "guidance_only"
    @test_throws ArgumentError VG.Guide([(1, 2), (1, 2)], zeros(2), [])
    @test_throws ArgumentError VG.Guide(atoms, linear, [(1, 2, 1.0), (2, 1, 1.0)])
    @test_throws ArgumentError VG.Guide(atoms, linear, [(1, 7, 1.0)])
    @test_throws ArgumentError VG.Guide(atoms, fill(Inf, 6), [])
    @test_throws ArgumentError VG.Guide(atoms, linear, [(1, 2, Inf)])
end

@testitem "Sparse value-pair floating-point and dense-move oracle" tags=[:guidance] begin
    using Random
    const VG = QUBOConstraints.ValuePairGuidance
    rng = Xoshiro(20261008)
    atoms = [(i, v) for i in 1:12 for v in (-7, 5, 11)]
    # Each Float64 input is represented exactly by BigFloat. The oracle uses
    # original semantic bindings, including diagonal and same-variable terms.
    linear = [ldexp(randn(rng), rand(rng, -15:15)) for _ in atoms]
    terms = [(a, b, ldexp(randn(rng), rand(rng, -15:15)))
        for a in eachindex(atoms) for b in a:length(atoms) if rand(rng) < 0.2]
    g = VG.Guide(atoms, linear, terms)
    w = VG.Workspace(g)
    polynomial(values) = sum(BigFloat(c) * (values[i] == v)
        for ((i, v), c) in zip(atoms, linear)) +
        sum(BigFloat(c) * (values[atoms[a][1]] == atoms[a][2]) *
            (values[atoms[b][1]] == atoms[b][2]) for (a, b, c) in terms)
    # A conservative bound covers folding diagonal coefficients and summation,
    # including cancellation. No bitwise equivalence with SIMD reductions is
    # required for arbitrary floating-point coefficients.
    tolerance = 8eps(Float64) * (length(linear) + length(terms)) *
        (sum(abs, linear) + sum(abs(c) for (_, _, c) in terms))
    for trial in 1:100
        values = rand(rng, (-7, 5, 11), 12)
        candidate = rand(rng, (-7, 5, 11), 12)
        ids = randperm(rng, 12)[1:mod1(trial, 12)]
        moved = copy(values)
        moved[ids] = candidate[ids]
        VG.refresh!(w, g, values)
        delta = VG.delta!(w, g, ids, candidate[ids])
        @test abs(BigFloat(delta) - (polynomial(moved) - polynomial(values))) <= tolerance
        @test abs(BigFloat(VG.energy(g, w)) - polynomial(values)) <= tolerance
        @test VG.delta!(w, g, reverse(ids), reverse(candidate[ids])) == delta
        @test allunique(w.touched_terms)
    end
end

@testitem "Sparse value-pair scope and proposal polynomial oracle" tags=[:guidance] begin
    using Random
    const VG = QUBOConstraints.ValuePairGuidance
    atoms = [(i, v) for i in 1:6 for v in (-3, 8, 12)]
    # A connected chain plus same-variable interactions exercises the exclusion
    # in field updates without dropping those terms from the full-move delta.
    terms = [(a, b, Float64(mod(a + 2b, 9) - 4) / 4)
        for a in eachindex(atoms) for b in a+1:length(atoms)
        if abs(atoms[a][1] - atoms[b][1]) <= 1]
    g = VG.Guide(atoms, zeros(length(atoms)), terms)
    w = VG.Workspace(g)
    polynomial(values) = sum(c * (values[atoms[a][1]] == atoms[a][2]) *
        (values[atoms[b][1]] == atoms[b][2]) for (a, b, c) in terms)
    for mode in ("absolute", "conditional"), depth in (1, 2, 4, 8), cap in (1, 7, 32), seed in 1:10
        values = [(-3, 8, 12)[mod1(seed + i, 3)] for i in 1:6]
        original = copy(values)
        move = VG.proposal!(w, g, values, depth, Xoshiro(seed); mode, max_candidates=cap)
        @test length(move.ids) == min(depth, 6)
        @test allunique(move.ids) && issorted(move.ids)
        @test move.examined <= cap
        candidate = copy(values)
        candidate[move.ids] = move.values
        @test move.delta == polynomial(candidate) - polynomial(values)
        @test move.delta <= 0.0
        @test values == original
        @test move.ids === w.chosen && move.values === w.replacements
        @test all(w.active[k] == Float64(values[i] == v) for (k, (i, v)) in enumerate(g.atoms))
    end
    disconnected = VG.Guide([(1, -3), (3, 8)], zeros(2), [])
    dw = VG.Workspace(disconnected)
    @test length(VG.scope!(dw, disconnected, [-3, 0, 8], 4, Xoshiro(1))) == 3
    empty = VG.Guide(Tuple{Int,Int}[], Float64[], [])
    ew = VG.Workspace(empty)
    @test isempty(VG.scope!(ew, empty, Int[], 2, Xoshiro(1)))
    @test VG.proposal!(ew, empty, Int[], 2, Xoshiro(1)).delta == 0.0
    @test_throws ArgumentError VG.scope!(w, g, fill(-3, 6), 0, Xoshiro(1))
    @test_throws ArgumentError VG.scope!(w, g, fill(-3, 6), 2, Xoshiro(1); exploration=2)
    @test_throws ArgumentError VG.proposal!(w, g, fill(-3, 6), 2, Xoshiro(1); max_candidates=0)
end

@testitem "Sparse and full-scan guidance branches and reserved buffers" tags=[:guidance] begin
    using Random
    const VG = QUBOConstraints.ValuePairGuidance
    n = 128
    atoms = [(i, v) for i in 1:n for v in (0, 1)]
    terms = [(2(i - 1) + a, 2i + b, a == b ? 0.5 : -0.25)
        for i in 1:n-1 for a in 1:2 for b in 1:2]
    g = VG.Guide(atoms, zeros(length(atoms)), terms)
    w = VG.Workspace(g)
    values = Int.(isodd.(1:n))
    polynomial(values) = sum(c * (values[atoms[a][1]] == atoms[a][2]) *
        (values[atoms[b][1]] == atoms[b][2]) for (a, b, c) in terms)
    VG.refresh!(w, g, values)
    for ids in ([1, 2], collect(1:n), [n], collect(n:-1:1), Int[])
        replacements = [1 - values[i] for i in ids]
        moved = copy(values); moved[ids] = replacements
        @test VG.delta!(w, g, ids, replacements) == polynomial(moved) - polynomial(values)
        @test length(ids) == n ? w.pair_visits == length(g.coefficient) :
            w.pair_visits == length(w.touched_terms)
        @test w.active == Float64[values[i] == v for (i, v) in atoms]
    end
    allocations(w, g, values, depth, rng) = @allocated VG.proposal!(w, g, values, depth, rng; mode="conditional")
    allocations(w, g, values, n, Xoshiro(1))
    fresh = VG.Workspace(g)
    @test allocations(fresh, g, values, 8, Xoshiro(2)) == 0
    @test allocations(fresh, g, values, n, Xoshiro(3)) == 0
    @test fresh.chosen !== w.chosen && fresh.replacements !== w.replacements
end
