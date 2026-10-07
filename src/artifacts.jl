_bit_artifact(b::BitID) = Dict("owner" => String(b.owner), "index" => b.index, "role" => String(b.role))

function _value_artifact(v)
    v isa Union{Integer,Rational,AbstractFloat,AbstractString,Symbol,Bool} ||
        v === nothing || throw(ArgumentError("artifact value type not supported: $(typeof(v))"))
    v isa Real && !isfinite(v) && throw(ArgumentError("nonfinite semantic value"))
    return Dict("type" => string(typeof(v)), "value" => v === nothing ? "nothing" : string(v))
end

_book_payload(book::Codebook) = Dict("variable" => String(book.variable),
    "bits" => _bit_artifact.(book.bits), "values" => _value_artifact.(book.values),
    "codes" => [Int.(collect(col)) for col in eachcol(book.codes)],
    "value_indices" => book.value_indices)

_book_payload(book::StructuredCodebook) = Dict("variable"=>String(book.variable),
    "bits"=>_bit_artifact.(book.bits),"values"=>_value_artifact.(book.values),
    "encoding"=>String(book.encoding),"weights"=>string.(book.weights),
    "coefficient_bound"=>string(book.coefficient_bound),"representation"=>"structured/1")

function _component_payload(q::QUBOComponent{T}) where T
    return Dict{String,Any}(
        "schema_version" => any(b->b isa StructuredCodebook,q.codebooks) ? "qubo-component/2" : "qubo-component/1",
        "convention" => "offset+linear+strict-upper-quadratic",
        "coefficient_type" => string(T), "offset" => string(q.offset),
        "bits" => _bit_artifact.(q.bits),
        "linear" => [Dict("i" => i, "coefficient" => string(v)) for (i, v) in _linear_terms(q)],
        "quadratic" => [Dict("i" => i, "j" => j, "coefficient" => string(v))
            for (i, j, v) in _quadratic_terms(q)],
        "codebooks" => [_book_payload(book) for book in q.codebooks],
        "auxiliaries" => [Dict("bit" => _bit_artifact(b), "meaning" => q.auxiliary_meanings[b])
            for b in q.bits if b.role !== :primary],
        "applicability" => q.applicability, "provenance" => q.provenance)
end

_component_digest(q) = bytes2hex(SHA.sha256(sprint(io -> TOML.print(io, _component_payload(q); sorted = true))))

"""Portable TOML-compatible artifact; exact coefficients are decimal/rational strings.

Certificates are bound to the complete component payload, including codebooks.
This is an export format; it does not execute code or import claimed proofs.
"""
function component_artifact(q::QUBOComponent; report::Union{Nothing,ExhaustiveReport} = nothing)
    data = _component_payload(q)
    data["sha256"] = _component_digest(q)
    if report !== nothing
        report.component_digest == data["sha256"] || throw(ArgumentError("certificate does not match component"))
        evidence = Dict{String,Any}("status" => String(report.status),
            "profile" => String(report.profile), "proof" => String(report.proof),
            "oracle_id" => report.oracle_id, "gap" => string(report.gap),
            "atol" => string(report.atol), "evaluations" => report.evaluations,
            "valid_codes" => report.valid_codes, "invalid_codes" => report.invalid_codes,
            "representation_policy" => "every_valid_code",
            "reason" => report.reason)
        report.min_violation === nothing || (evidence["min_violation"] = string(report.min_violation))
        if report.counterexample !== nothing
            c = report.counterexample
            evidence["counterexample"] = Dict("primary_mask" => c.primary_mask,
                "bits" => Int.(c.bits), "values" => _value_artifact.(c.values),
                "valid_code" => c.valid_code, "satisfied" => c.satisfied, "energy" => string(c.energy))
        end
        data["certificate"] = evidence
    end
    return data
end

function write_component(io::IO, q::QUBOComponent; report = nothing)
    TOML.print(io, component_artifact(q; report); sorted = true)
    return nothing
end
