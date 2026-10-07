"""A constructive compiler exceeded its declared budget; never a synthesis UNSAT claim."""
struct CompilationLimit <: Exception
    stage::Symbol
    requested::BigInt
    limit::BigInt
end
Base.showerror(io::IO,e::CompilationLimit) = print(io,"unsupported_under_budget: ",e.stage,
    " requested ",e.requested," > ",e.limit)

"""Exact finite truth-table construction with optional Rosenberg quadratization.

Every invalid primary code is assigned target 1; every valid code is assigned !predicate.
This is STRONGER than the public oracle's invalid-code >=1 contract. Degree>2 of this
fixed extension does not prove that all other valid-code-equivalent QUBOs need auxiliaries.
Returns an exact constructive upper bound, not an auxiliary-minimal learned formulation.

The callback receives a reusable semantic values vector in codebook order and must be pure
and return Bool. Budgets apply before full-cube allocation and during substitution.
Private auxiliaries implement AND gates; all coefficients are BigInt.
"""
function truth_component(books,predicate;oracle_id::AbstractString,
        max_primary_states::Integer=4096,max_terms::Integer=100_000,
        max_auxiliaries::Integer=1024,max_work::Integer=1_000_000,
        zero_aux_only::Bool=false,namespace::Symbol=:truth_quadratization)
    isempty(oracle_id) && throw(ArgumentError("oracle identity is required"))
    books = collect(books)
    all(book->book isa AbstractCodebook,books) || throw(ArgumentError("codebooks required"))
    max_primary_states>=1 && max_terms>=1 && max_auxiliaries>=0 && max_work>=1 || throw(ArgumentError("invalid compilation budget"))
    allunique(book.variable for book in books) || throw(ArgumentError("duplicate codebook variable"))
    bits = reduce(vcat,(_book_bits(book) for book in books);init=BitID[])
    allunique(bits) || throw(ArgumentError("overlapping codebooks"))
    n = length(bits)
    n<8sizeof(Int)-1 || throw(CompilationLimit(:primary_bits,big(n),big(8sizeof(Int)-2)))
    states = 1<<n
    states<=max_primary_states || throw(CompilationLimit(:primary_states,big(states),big(max_primary_states)))
    positions = Dict(b=>i for (i,b) in enumerate(bits))
    book_positions = Vector{Int}[[positions[b] for b in _book_bits(book)] for book in books]
    buffers = BitVector[falses(length(_book_bits(book))) for book in books]
    values = Vector{Any}(undef,length(books))
    book_values = Vector[_book_values(book) for book in books]
    coefficients = fill(big(1),states)
    zero_value,one_value = big(0),big(1)
    for mask in 0:states-1
        valid = true
        for (k,book) in enumerate(books)
            for (j,i) in enumerate(book_positions[k])
                buffers[k][j] = !iszero(mask & (1<<(i-1)))
            end
            index = _value_index(book,buffers[k])
            if index==0
                valid = false
                break
            end
            values[k] = book_values[k][index]
        end
        if valid
            satisfied = predicate(values)
            satisfied isa Bool || throw(ArgumentError("truth predicate must return Bool"))
            coefficients[mask+1] = satisfied ? zero_value : one_value
        end
    end
    truth_digest = bytes2hex(SHA.sha256(join(string.(coefficients))))
    for bit in 0:n-1, mask in 0:states-1
        !iszero(mask & (1<<bit)) && (coefficients[mask+1] -= coefficients[mask-(1<<bit)+1])
    end
    offset = coefficients[1]
    polynomial = Dict{Tuple{Vararg{Int}},BigInt}()
    for mask in 1:states-1
        iszero(coefficients[mask+1]) && continue
        monomial = Tuple(i for i in 1:n if !iszero(mask & (1<<(i-1))))
        polynomial[monomial] = coefficients[mask+1]
    end
    length(polynomial)<=max_terms || throw(CompilationLimit(:polynomial_terms,big(length(polynomial)),big(max_terms)))
    degree = maximum(length,keys(polynomial);init=0)
    zero_aux_only && degree>2 && throw(CompilationLimit(:fixed_extension_degree,big(degree),big(2)))
    magnitude = sum(abs,Base.values(polynomial);init=big(0))
    weight = magnitude+2
    gates = Tuple{Int,Int,Int}[]
    meanings = Dict{BitID,String}()
    work = big(0)
    while true
        higher = [key for key in keys(polynomial) if length(key)>2]
        isempty(higher) && break
        sort!(higher;by=key->(-length(key),key))
        a,b = first(higher)[1:2]
        length(gates)<max_auxiliaries || throw(CompilationLimit(:auxiliaries,big(length(gates)+1),big(max_auxiliaries)))
        work += length(polynomial)
        work<=max_work || throw(CompilationLimit(:substitution_work,work,big(max_work)))
        auxiliary = BitID(namespace,length(gates)+1;role=:quadratization_auxiliary)
        push!(bits,auxiliary)
        y = length(bits)
        push!(gates,(a,b,y))
        meanings[auxiliary] = "AND($(bits[a].role):$(bits[a].owner):$(bits[a].index), $(bits[b].role):$(bits[b].owner):$(bits[b].index))"
        replacement = Dict{Tuple{Vararg{Int}},BigInt}()
        # Replace this pair in ALL terms. An already replaced pair never reappears.
        for (key,c) in polynomial
            updated = if a in key && b in key
                Tuple(sort!([i for i in key if i!=a && i!=b] ∪ [y]))
            else
                key
            end
            replacement[updated] = get(replacement,updated,big(0))+c
        end
        filter!(p->!iszero(last(p)),replacement)
        polynomial = replacement
    end
    length(polynomial)+4big(length(gates))<=max_terms ||
        throw(CompilationLimit(:qubo_terms,big(length(polynomial))+4length(gates),big(max_terms)))
    linear = [(key[1],c) for (key,c) in polynomial if length(key)==1]
    quadratic = [(key[1],key[2],c) for (key,c) in polynomial if length(key)==2]
    for (a,b,y) in gates
        push!(linear,(y,3weight))
        append!(quadratic,[(a,b,weight),(a,y,-2weight),(b,y,-2weight)])
    end
    # Any inconsistent gate costs >=weight; the remaining polynomial is >=-magnitude.
    # Therefore it costs >=2, above every target (0 or 1). Consistent gates recover the table.
    return QUBOComponent(bits;linear,quadratic,offset,codebooks=books,auxiliary_meanings=meanings,
        applicability="finite fixed-cube indicator; invalid codes exactly 1 after minimization; no minimality claim",
        provenance="truth/Rosenberg v1; oracle=$(oracle_id); primary_states=$(states); truth_sha256=$(truth_digest); original_degree=$(degree); auxiliary_count=$(length(gates)); guard_weight=$(weight)")
end
