"""Sparse guidance on explicitly bound original `(variable, value)` atoms.

Guidance energy is a proxy and never replaces original constraint validation.
Guide arrays are owned and must be treated as read-only. Workspaces and their
returned scope/proposal buffers belong to one caller and may be reused.
"""
module ValuePairGuidance
using Random
export Guide, Workspace, refresh!, energy, delta!, scope!, proposal!, metadata

struct Guide
    atoms::Vector{Tuple{Int,Int}}
    linear::Vector{Float64}
    left::Vector{Int}
    right::Vector{Int}
    coefficient::Vector{Float64}
    by_variable::Vector{Vector{Int}}
    by_atom_terms::Vector{Vector{Int}}
    by_variable_terms::Vector{Vector{Int}}
    absolute_strength::Vector{Float64}
    provenance::String
    source_sha256::String
end

function Guide(atoms, linear, terms; provenance="unspecified guidance", source_sha256="",
        max_atoms=2048, max_terms=20000)
    a = Tuple{Int,Int}[(Int(i), Int(v)) for (i, v) in atoms]
    length(a) <= max_atoms && length(terms) <= max_terms || throw(ArgumentError("QUBO guide size cap exceeded"))
    allunique(a) && all(x -> x[1] > 0, a) || throw(ArgumentError("duplicate or invalid value atom"))
    length(linear) == length(a) && all(isfinite, linear) || throw(ArgumentError("invalid linear guide coefficients"))
    l = Float64.(linear)
    left, right, coefficient = Int[], Int[], Float64[]
    seen = Set{Tuple{Int,Int}}()
    for (i, j, c) in terms
        1 <= i <= length(a) && 1 <= j <= length(a) && isfinite(c) || throw(ArgumentError("invalid pair coefficient"))
        i, j = minmax(Int(i), Int(j))
        (i, j) in seen && throw(ArgumentError("duplicate canonical pair coefficient"))
        push!(seen, (i, j))
        if i == j
            l[i] += c
        elseif !iszero(c)
            push!(left, i); push!(right, j); push!(coefficient, Float64(c))
        end
    end
    all(isfinite, l) && all(isfinite, coefficient) || throw(ArgumentError("coefficient conversion overflow"))
    n = isempty(a) ? 0 : maximum(first, a)
    n <= max_atoms || throw(ArgumentError("variable index exceeds guide cap"))
    by = [Int[] for _ in 1:n]
    for (k, (i, _)) in enumerate(a)
        push!(by[i], k)
    end
    atom_terms = [Int[] for _ in a]
    variable_terms = [Int[] for _ in by]
    strength = zeros(n)
    for k in eachindex(coefficient)
        i, j = left[k], right[k]
        u, v = a[i][1], a[j][1]
        push!(atom_terms[i], k); push!(atom_terms[j], k)
        push!(variable_terms[u], k)
        u == v || push!(variable_terms[v], k)
        c = abs(coefficient[k])
        strength[u] += c; strength[v] += c
    end
    Guide(a, l, left, right, coefficient, by, atom_terms, variable_terms, strength,
        String(provenance), String(source_sha256))
end

mutable struct Workspace
    guide::Guide
    active::Vector{Float64}
    change::Vector{Float64}
    strength::Vector{Float64}
    affinity::Vector{Float64}
    chosen::Vector{Int}
    replacements::Vector{Int}
    selected::Vector{Bool}
    field::Vector{Float64}
    changed_atoms::Vector{Int}
    touched_terms::Vector{Int}
    term_epoch::Vector{UInt}
    epoch::UInt
    pair_visits::Int
end

Workspace(g::Guide) = Workspace(g, zeros(length(g.atoms)), zeros(length(g.atoms)),
    zeros(length(g.by_variable)), zeros(length(g.by_variable)),
    sizehint!(Int[], length(g.by_variable)), sizehint!(Int[], length(g.by_variable)),
    fill(false, length(g.by_variable)), zeros(length(g.atoms)),
    sizehint!(Int[], length(g.atoms)), sizehint!(Int[], length(g.coefficient)),
    zeros(UInt, length(g.coefficient)), zero(UInt), 0)

function _check_guide(w::Workspace, g::Guide)
    w.guide === g || throw(ArgumentError("workspace belongs to another guide"))
    nothing
end

function refresh!(w::Workspace, g::Guide, values)
    _check_guide(w, g)
    length(values) >= length(g.by_variable) || throw(DimensionMismatch("guide parent variables"))
    for k in eachindex(g.atoms)
        i, v = g.atoms[k]
        w.active[k] = Float64(values[i] == v)
    end
    w
end

function energy(g::Guide, w::Workspace)
    _check_guide(w, g)
    total = 0.0
    @inbounds @simd for i in eachindex(g.linear)
        total += g.linear[i] * w.active[i]
    end
    @inbounds @simd for k in eachindex(g.coefficient)
        total += g.coefficient[k] * w.active[g.left[k]] * w.active[g.right[k]]
    end
    total
end

"Simultaneous polynomial delta; shared changed/changed edges are visited once."
function delta!(w::Workspace, g::Guide, ids, replacements)
    _check_guide(w, g)
    length(ids) == length(replacements) || throw(ArgumentError("invalid move scope"))
    # Validate before touching scratch, including duplicate variables at small depths.
    for k in eachindex(ids)
        i = ids[k]
        1 <= i <= length(g.by_variable) || throw(BoundsError(g.by_variable, i))
        for j in firstindex(ids):k-1
            i != ids[j] || throw(ArgumentError("duplicate move variable"))
        end
    end
    for a in w.changed_atoms
        w.change[a] = 0.0
    end
    empty!(w.changed_atoms); empty!(w.touched_terms)
    w.pair_visits = 0
    w.epoch += one(UInt)
    if iszero(w.epoch)
        fill!(w.term_epoch, zero(UInt))
        w.epoch = one(UInt)
    end
    adjacency_work = 0
    for (i, v) in zip(ids, replacements), a in g.by_variable[i]
        change = Float64(g.atoms[a][2] == v) - w.active[a]
        iszero(change) && continue
        w.change[a] = change
        push!(w.changed_atoms, a)
        adjacency_work += length(g.by_atom_terms[a])
    end
    # Sorting a wide edge frontier costs more than scanning the full polynomial.
    # This independent reduction has no neighbor writes or fast-math. Both paths
    # use ordinary Float64 arithmetic and preserve the same simultaneous move.
    if adjacency_work > 0 && adjacency_work >= div(length(g.coefficient), 8)
        w.pair_visits = length(g.coefficient)
        total = 0.0
        @inbounds @simd for a in eachindex(g.linear)
            total += g.linear[a] * w.change[a]
        end
        @inbounds @simd for k in eachindex(g.coefficient)
            i, j = g.left[k], g.right[k]
            total += g.coefficient[k] * (w.change[i] * w.active[j] +
                w.active[i] * w.change[j] + w.change[i] * w.change[j])
        end
        return total
    end
    # Canonical linear traversal is independent of the caller's move ordering.
    sort!(w.changed_atoms; alg=Base.Sort.QuickSort)
    total = 0.0
    for a in w.changed_atoms
        total += g.linear[a] * w.change[a]
    end
    for a in w.changed_atoms
        for k in g.by_atom_terms[a]
            w.term_epoch[k] == w.epoch && continue
            w.term_epoch[k] = w.epoch
            push!(w.touched_terms, k)
        end
    end
    # Canonical pair traversal is independent of the caller's move ordering.
    # Julia's default integer sort may allocate radix/scratch buffers. These
    # in-place algorithms retain canonical order without per-move scratch.
    sort!(w.touched_terms; alg=Base.Sort.QuickSort)
    w.pair_visits = length(w.touched_terms)
    for k in w.touched_terms
        i, j = g.left[k], g.right[k]
        total += g.coefficient[k] * (w.change[i] * w.active[j] +
            w.active[i] * w.change[j] + w.change[i] * w.change[j])
    end
    total
end

function scope!(w::Workspace, g::Guide, values, depth, rng; mode="absolute", exploration=0.10)
    mode in ("absolute", "conditional") && 0 <= exploration <= 1 && depth > 0 || throw(ArgumentError("invalid guide policy"))
    refresh!(w, g, values)
    fill!(w.affinity, 0.0); fill!(w.selected, false); empty!(w.chosen)
    if mode == "absolute"
        copyto!(w.strength, g.absolute_strength)
    else
        fill!(w.strength, 0.0)
        for k in eachindex(g.coefficient)
            i, j = g.left[k], g.right[k]
            a, b = g.atoms[i][1], g.atoms[j][1]
            c = abs(g.coefficient[k]) * ((w.active[i] + w.active[j]) / 2)
            w.strength[a] += c; w.strength[b] += c
        end
    end
    n = length(w.strength)
    n == 0 && return w.chosen
    seed = rand(rng) < exploration || maximum(w.strength) == 0 ? rand(rng, 1:n) : argmax(w.strength)
    for _ in 1:min(depth, n)
        w.selected[seed] = true; push!(w.chosen, seed)
        for k in g.by_variable_terms[seed]
            a, b = g.atoms[g.left[k]][1], g.atoms[g.right[k]][1]
            c = abs(g.coefficient[k]) * (mode == "conditional" ? (w.active[g.left[k]] + w.active[g.right[k]]) / 2 : 1.0)
            if a == seed
                w.affinity[b] += c
            elseif b == seed
                w.affinity[a] += c
            end
        end
        best, quality = 0, -Inf
        for i in 1:n
            w.selected[i] && continue
            q = w.affinity[i] + 1e-6 * w.strength[i]
            if q > quality
                best, quality = i, q
            end
        end
        best == 0 && break
        seed = best
    end
    sort!(w.chosen; alg=Base.Sort.QuickSort)
end

"Bounded sequential field proposal; returned arrays are reused on the next call."
function proposal!(w::Workspace, g::Guide, values, depth, rng; mode="absolute", exploration=0.1, max_candidates=32)
    max_candidates > 0 || throw(ArgumentError("positive guide candidate cap required"))
    scope!(w, g, values, depth, rng; mode, exploration)
    empty!(w.replacements)
    copyto!(w.field, g.linear)
    for k in eachindex(g.coefficient)
        a, b = g.left[k], g.right[k]
        g.atoms[a][1] == g.atoms[b][1] && continue
        w.field[a] += g.coefficient[k] * w.active[b]
        w.field[b] += g.coefficient[k] * w.active[a]
    end
    examined = 0
    for i in w.chosen
        best, oldfield = values[i], 0.0
        for k in g.by_variable[i]
            oldfield += w.active[k] * w.field[k]
        end
        bestfield = oldfield
        for k in g.by_variable[i]
            examined >= max_candidates && break
            candidate = g.atoms[k][2]; examined += 1
            if w.field[k] < bestfield
                best, bestfield = candidate, w.field[k]
            end
        end
        push!(w.replacements, Int(best))
        for a in g.by_variable[i]
            change = Float64(g.atoms[a][2] == best) - w.active[a]
            iszero(change) && continue
            for k in g.by_atom_terms[a]
                b = g.left[k] == a ? g.right[k] : g.left[k]
                g.atoms[b][1] == i && continue
                w.field[b] += g.coefficient[k] * change
            end
            w.active[a] += change
        end
    end
    refresh!(w, g, values)
    (; ids=w.chosen, values=w.replacements, examined,
        delta=delta!(w, g, w.chosen, w.replacements))
end

metadata(g::Guide) = Dict("provenance" => g.provenance, "source_sha256" => g.source_sha256,
    "atoms" => length(g.atoms), "pair_terms" => length(g.coefficient), "authority" => "guidance_only",
    "encoding" => "explicit_one_hot_value_atoms", "convention" => "canonical_polynomial")
end
