"""Exact nonnegative one-hot formulations for the five initial research families.

Supported families: `:all_different`, `:ordered`, `:linear_sum`, `:no_overlap`,
`:channel`. Domains are explicit finite lists; channel values are indices 1:n.
`ordered` means nondecreasing. `no_overlap` uses fixed nonnegative integer lengths
(zero-length tasks are not ignored). Coefficients and rhs for linear_sum are integers.
No auxiliary bits or learned/global-minimum estimates are used.
"""
function one_hot_constraint(family::Symbol, domains;
        variables = [Symbol(:x,i) for i in eachindex(domains)],
        rhs::Integer = 0, coefficients = ones(Int,length(domains)),
        lengths = ones(Int,length(domains)))
    family in (:all_different,:ordered,:linear_sum,:no_overlap,:channel) ||
        throw(ArgumentError("unsupported one-hot family"))
    n = length(domains)
    n > 0 || throw(ArgumentError("at least one variable required"))
    length(variables) == n && allunique(variables) || throw(ArgumentError("invalid variable identities"))
    books = [codebook(variables[i],domains[i]) for i in 1:n]
    bits = reduce(vcat,(b.bits for b in books))
    positions = Dict(b=>i for (i,b) in enumerate(bits))
    indices = [[positions[b] for b in book.bits] for book in books]
    linear = Tuple{Int,BigInt}[]
    quadratic = Tuple{Int,Int,BigInt}[]
    offset = big(n)
    # Σ_i(Σ_a z_ia - 1)^2. All later summands are nonnegative on ALL bitstrings.
    for ids in indices
        append!(linear,((i,big(-1)) for i in ids))
        append!(quadratic,((ids[a],ids[b],big(2)) for a in eachindex(ids) for b in a+1:length(ids)))
    end
    if family === :linear_sum
        length(coefficients)==n && all(c->c isa Integer,coefficients) ||
            throw(ArgumentError("one integer coefficient per variable required"))
        all(b->all(v->v isa Integer,b.values),books) || throw(ArgumentError("integer sum domains required"))
        weights = BigInt[big(coefficients[i])*big(v) for (i,b) in enumerate(books) for v in b.values]
        target = big(rhs)
        offset += target^2
        append!(linear,((i,w*w-2target*w) for (i,w) in enumerate(weights)))
        append!(quadratic,((i,j,2weights[i]*weights[j]) for i in eachindex(weights) for j in i+1:length(weights)))
    elseif family === :channel
        all(b->all(v->v isa Integer && 1<=v<=n,b.values),books) ||
            throw(ArgumentError("channel domains must be subsets of 1:n"))
        for i in 1:n, j in i+1:n
            a = findfirst(==(j),books[i].values)
            b = findfirst(==(i),books[j].values)
            # Missing semantic values have indicator zero, not an extra binary bit.
            a === nothing || push!(linear,(indices[i][a],big(1)))
            b === nothing || push!(linear,(indices[j][b],big(1)))
            if a !== nothing && b !== nothing
                push!(quadratic,(indices[i][a],indices[j][b],big(-2)))
            end
        end
    else
        if family === :no_overlap
            length(lengths)==n && all(l->l isa Integer && l>=0,lengths) ||
                throw(ArgumentError("one fixed nonnegative integer length per task required"))
            all(b->all(v->v isa Integer,b.values),books) || throw(ArgumentError("integer time domains required"))
        end
        for i in 1:n, j in i+1:n
            family === :ordered && j != i+1 && continue
            for (a,x) in enumerate(books[i].values), (b,y) in enumerate(books[j].values)
                bad = if family === :all_different
                    isequal(x,y)
                elseif family === :ordered
                    !(x <= y)
                else
                    !(big(x)+big(lengths[i]) <= big(y) || big(y)+big(lengths[j]) <= big(x))
                end
                bad && push!(quadratic,(indices[i][a],indices[j][b],big(1)))
            end
        end
    end
    return QUBOComponent{BigInt}(bits;linear,quadratic,offset,codebooks=books,
        applicability="finite one-hot; nonnegative construction; family=$(family)",
        provenance="Q2 exact construction v1; family=$(family); rhs=$(rhs); coefficients=$(coefficients); lengths=$(lengths)")
end
