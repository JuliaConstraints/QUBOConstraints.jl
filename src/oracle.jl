"""Finite-domain validation evidence bound to a component fingerprint and named oracle.

`pass/exhaustive` concerns only the enumerated instance. `validated` denotes
floating-point/tolerance checks. Counterexamples follow canonical primary-bit
integer order, not a distance-minimal repair order.
"""
struct ExhaustiveReport
    status::Symbol
    profile::Symbol
    proof::Symbol
    component_digest::String
    oracle_id::String
    gap::Real
    atol::Real
    evaluations::Int
    valid_codes::Int
    invalid_codes::Int
    min_violation::Union{Nothing,Real}
    counterexample::Union{Nothing,NamedTuple}
    reason::String
end

"""Extension point: `constraint_oracle(:all_different; kwargs...)` after loading Constraints."""
function constraint_oracle end

_exact_bound(value::Integer) = BigInt(value)
_exact_bound(value::Real) = Rational{BigInt}(value)

mutable struct _OracleCounts{T}
    evaluations::Int
    valid_codes::Int
    invalid_codes::Int
    min_violation::Union{Nothing,T}
end

function _set_bits!(z, positions, mask)
    for (shift, i) in enumerate(positions)
        z[i] = !iszero(mask & (1 << (shift - 1)))
    end
    return z
end

"""Exhaustively minimize auxiliaries for EVERY valid primary code, checking invalid codes too.

`predicate(values)` receives semantic values in `q.codebooks` order and must
return Bool. `gap` also bounds energies of invalid codes below. A budget/time
limit returns `unknown` and never an UNSAT result. Surrogates are not certified.
"""
function exhaustive_check(q::QUBOComponent{T}, predicate;
        oracle_id::AbstractString, profile::Symbol = :zero_set_exact,
        gap::Real = 1, atol::Real = 0, max_states::Integer = 1_048_576,
        time_limit::Real = Inf) where T
    profile in (:zero_set_exact, :indicator_exact, :surrogate) ||
        throw(ArgumentError("unsupported exactness profile"))
    isfinite(atol) && atol >= 0 || throw(ArgumentError("invalid tolerance"))
    isfinite(gap) || throw(ArgumentError("nonfinite gap"))
    # Do threshold arithmetic exactly even when bounds were supplied as floats.
    gap_bound, tolerance = _exact_bound(gap), _exact_bound(atol)
    gap_bound > 2tolerance || throw(ArgumentError("gap must exceed twice the tolerance"))
    lower_bound = gap_bound - tolerance
    profile === :indicator_exact && gap > 1 && throw(ArgumentError("indicator gap cannot exceed one"))
    max_states >= 1 && time_limit > 0 || throw(ArgumentError("invalid enumeration budget"))
    isempty(oracle_id) && throw(ArgumentError("oracle identity is required"))
    digest = _component_digest(q)
    counts = _OracleCounts{T}(0, 0, 0, nothing)
    report(status, proof, reason; counterexample = nothing) = ExhaustiveReport(
        status, profile, proof, digest, String(oracle_id), gap, atol,
        counts.evaluations, counts.valid_codes, counts.invalid_codes,
        counts.min_violation, counterexample, reason)
    profile === :surrogate && return report(:surrogate, :unknown, "surrogates do not receive exactness certificates")
    primary = findall(b -> b.role === :primary, q.bits)
    auxiliaries = findall(b -> b.role !== :primary, q.bits)
    covered = BitID[b for book in q.codebooks for b in book.bits]
    Set(covered) == Set(q.bits[primary]) || throw(ArgumentError("primary bits need complete codebooks"))
    n = length(q.bits)
    (n >= 8sizeof(Int) - 1 || big(2)^n > max_states) &&
        return report(:unknown, :unknown, "enumeration exceeds max_states")
    pos = Dict(b => i for (i, b) in enumerate(q.bits))
    book_positions = Vector{Int}[[pos[b] for b in _book_bits(book)] for book in q.codebooks]
    buffers = BitVector[falses(length(_book_bits(book))) for book in q.codebooks]
    book_values = Vector[_book_values(book) for book in q.codebooks]
    values = Vector{Any}(undef, length(q.codebooks))
    z = falses(n)
    started = time_ns()
    limited = isfinite(time_limit)
    expired() = limited && (time_ns() - started) / 1.0e9 >= time_limit
    for mask in 0:((1 << length(primary)) - 1)
        expired() && return report(:unknown, :unknown, "time limit reached")
        _set_bits!(z, primary, mask)
        valid = true
        for (k, book) in enumerate(q.codebooks)
            for (j, i) in enumerate(book_positions[k])
                buffers[k][j] = z[i]
            end
            index = _value_index(book, buffers[k])
            if index == 0
                valid = false
                break
            end
            values[k] = book_values[k][index]
        end
        valid ? (counts.valid_codes += 1) : (counts.invalid_codes += 1)
        satisfied = false
        if valid
            try
                verdict = predicate(values)
                verdict isa Bool || return report(:failed, :unknown, "oracle must return Bool")
                satisfied = verdict
            catch err
                err isa InterruptException && rethrow()
                return report(:failed, :unknown, "oracle failure: " * sprint(showerror, err))
            end
        end
        best = zero(T)
        witness = 0
        for a in 0:((1 << length(auxiliaries)) - 1)
            expired() && return report(:unknown, :unknown, "time limit reached")
            _set_bits!(z, auxiliaries, a)
            e = _energy(q, z)
            counts.evaluations += 1
            isfinite(e) || return report(:failed, :unknown, "nonfinite energy")
            if a == 0 || e < best
                best, witness = e, a
            end
        end
        if !valid || !satisfied
            counts.min_violation = counts.min_violation === nothing ? best : min(counts.min_violation, best)
        end
        ok = if !valid
            best >= lower_bound
        elseif satisfied
            abs(best) <= tolerance
        elseif profile === :indicator_exact
            abs(best - 1) <= tolerance
        else
            best >= lower_bound
        end
        if !ok
            _set_bits!(z, auxiliaries, witness)
            counterexample = (; primary_mask = mask, bits = collect(z),
                values = valid ? copy(values) : Any[], valid_code = valid,
                satisfied, energy = best)
            return report(:fail, :unknown,
                valid ? "valid-code energy violates requested profile" : "invalid code is not sufficiently penalized";
                counterexample)
        end
    end
    proof = T === Float64 || !iszero(atol) ? :validated : :exhaustive
    return report(:pass, proof, "all primary codes and auxiliary assignments enumerated; representation-uniform")
end
