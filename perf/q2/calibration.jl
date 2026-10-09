using QUBOConstraints

function calibration_qubo(kind::Symbol,n::Int)
    n>=1 || throw(ArgumentError("positive arity required"))
    kind in (:positive_monomial,:negative_monomial_shifted,:even_parity) || throw(ArgumentError("unknown calibration"))
    books = [codebook(Symbol(:x,i),0:1;encoding=:domain_wall) for i in 1:n]
    bits = reduce(vcat,(b.bits for b in books))
    linear = Tuple{Int,BigInt}[]
    quadratic = Tuple{Int,Int,BigInt}[]
    offset = big(0)
    if n==1
        offset = kind===:positive_monomial ? big(0) : big(1)
        push!(linear,(1,kind===:positive_monomial ? big(1) : big(-1)))
    elseif kind===:negative_monomial_shifted
        offset = big(1)
        if n==2
            push!(quadratic,(1,2,big(-1)))
        else
            push!(bits,BitID(:y,1;role=:quadratization_auxiliary))
            push!(linear,(n+1,big(n-1)))
            append!(quadratic,((i,n+1,big(-1)) for i in 1:n))
        end
    else
        ell = ndigits(n-1;base=2)
        weights = vcat(ones(BigInt,n),[-big(2)^i for i in 1:ell-1])
        append!(bits,[BitID(:y,i;role=:quadratization_auxiliary) for i in 1:ell-1])
        # Odd roots 1,3,...,2^ell-1 cover every odd weight up to n. This also
        # avoids the odd-n endpoint defect of Eq.27 in the 2018 preprint.
        constant = kind===:even_parity ? big(-1) : big(2)^ell-n
        if kind===:positive_monomial
            offset = div(constant*(constant-1),2)
            append!(linear,((i,div(w*w+(2constant-1)*w,2)) for (i,w) in enumerate(weights)))
            append!(quadratic,((i,j,weights[i]*weights[j]) for i in eachindex(weights) for j in i+1:length(weights)))
        else
            offset = constant^2
            append!(linear,((i,w*w+2constant*w) for (i,w) in enumerate(weights)))
            append!(quadratic,((i,j,2weights[i]*weights[j]) for i in eachindex(weights) for j in i+1:length(weights)))
        end
    end
    return QUBOComponent{BigInt}(bits;linear,quadratic,offset,codebooks=books,
        provenance="Boros Crama Rodriguez-Heck compact quadratizations; $(kind); n=$(n); negative monomial shifted by +1")
end

calibration_truth(kind,v) = kind===:positive_monomial ? !all(==(1),v) :
    kind===:negative_monomial_shifted ? all(==(1),v) : isodd(sum(v))
