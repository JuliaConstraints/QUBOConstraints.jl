"""Exact nonnegative validity QUBO for a structured encoding, with unit minimum gap.

Minimization is over explicitly identified validity witnesses, never over primary
representations. Binary/Fibonacci range constraints use bounded binary slack; Gray
uses prefix-XOR gates and range slack. This proves validity, not arbitrary constraints
or safety when adding an unbounded negative objective. No minimum-auxiliary claim.
"""
function encoding_validity(book::StructuredCodebook;max_terms::Integer=1_000_000)
    max_terms>0 || throw(ArgumentError("positive term budget required"))
    bits = copy(book.bits)
    linear = Tuple{Int,BigInt}[]
    quadratic = Tuple{Int,Int,BigInt}[]
    offset = big(0)
    meanings = Dict{BitID,String}()
    function auxiliary(label)
        bit = BitID(Symbol(book.variable,"/validity"),length(meanings)+1;role=:quadratization_auxiliary)
        push!(bits,bit)
        meanings[bit] = label
        return length(bits)
    end
    function square!(terms,c=0)
        length(linear)+length(quadratic)+length(terms)*(big(length(terms))+1)÷2<=max_terms ||
            throw(ArgumentError("validity generation exceeds max_terms"))
        offset += big(c)^2
        for (i,w) in terms
            push!(linear,(i,big(w)^2+2big(c)*w))
        end
        for j in eachindex(terms), i in 1:j-1
            push!(quadratic,(terms[i][1],terms[j][1],2big(terms[i][2])*terms[j][2]))
        end
    end
    mode = book.encoding
    width = length(book.bits)
    K = length(book.values)-1
    # Range on already locally valid primary codes, not on the unrestricted weight sum.
    max_rank = mode===:fibonacci ? sum(book.weights[width:-2:1];init=big(0)) :
        mode in (:native,:binary,:gray) ? sum(book.weights;init=big(0)) : big(K)
    needs_range = max_rank>K
    terms = [(i,book.weights[i]) for i in 1:width]
    if mode===:one_hot
        square!([(i,big(1)) for i in 1:width],-1)
    elseif mode===:zero_one_hot
        big(width)*(width-1)÷2<=max_terms || throw(ArgumentError("validity generation exceeds max_terms"))
        append!(quadratic,((i,j,big(1)) for i in 1:width for j in i+1:width))
    elseif mode===:domain_wall
        2big(max(0,width-1))<=max_terms || throw(ArgumentError("validity generation exceeds max_terms"))
        for i in 1:width-1
            push!(linear,(i+1,big(1)))
            push!(quadratic,(i,i+1,big(-1)))
        end
    else
        if mode===:fibonacci
            append!(quadratic,((i,i+1,big(1)) for i in 1:width-1))
        elseif mode===:gray && width>0 && needs_range
            # Most significant binary bit equals the most significant Gray bit.
            decoded = collect(1:width)
            for i in width-1:-1:1
                a,b = decoded[i+1],i
                u = auxiliary("AND of decoded bit $(i+1) and Gray bit $(i)")
                d = auxiliary("decoded binary bit $(i)")
                # Rosenberg AND penalty: ab-2au-2bu+3u >= 0 on all bits.
                append!(quadratic,[(a,b,big(1)),(a,u,big(-2)),(b,u,big(-2))])
                push!(linear,(u,big(3)))
                square!([(a,big(1)),(b,big(1)),(u,big(-2)),(d,big(-1))])
                decoded[i] = d
            end
            terms = [(decoded[i],book.weights[i]) for i in 1:width]
        end
        if needs_range
            slack = structured_codebook(:slack,0:K;encoding=:bounded_binary)
            append!(terms,[(auxiliary("range slack weight $(w)"),w) for w in slack.weights])
            square!(terms,-K)
        end
    end
    length(linear)+length(quadratic)<=max_terms || throw(ArgumentError("validity generation exceeds max_terms"))
    return QUBOComponent(bits;linear,quadratic,offset,codebooks=[book],auxiliary_meanings=meanings,
        applicability="validity only; every primary representation; minimize private witnesses",
        provenance="structured validity v1; encoding=$(mode); integer nonnegative factors")
end
