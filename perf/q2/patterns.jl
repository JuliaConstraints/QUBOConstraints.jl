using QUBOConstraints, TOML, SHA

function square_pattern(id,a,b,i,j,d,p)
    id==1 && return Int(a==b)
    id==2 && return Int(a!=b)
    id==3 && return Int(a>b)
    id==4 && return Int(a>=b)
    id==5 && return Int(a<b)
    id==6 && return Int(a<=b)
    own = a==i || b==j
    outside = a!=i && a!=j && b!=i && b!=j
    id==7 && return -Int(own)
    id==8 && return Int(own)
    id==9 && return -Int(outside)
    id==10 && return Int(outside)
    id==11 && return a==j && b==i ? -1 : Int(a==j || b==i)
    id==12 && return max(0,p-abs(a-b))
    id==13 && return max(0,abs(a-b)-(d-1-p))
    id==14 && return 2a*b
    error("unknown square pattern")
end

function raw_pattern_component(selection,triangle,n,d,p;offset=0)
    books = [codebook(Symbol(:x,i),1:d) for i in 1:n]
    bits = reduce(vcat,(b.bits for b in books))
    linear = Tuple{Int,BigInt}[]
    quadratic = Tuple{Int,Int,BigInt}[]
    for i in 1:n, a in 1:d
        id = (i-1)*d+a
        value = triangle==1 ? 0 : triangle==2 ? -1 : a*a-2p*a
        push!(linear,(id,big(value)))
        for b in a+1:d
            push!(quadratic,(id,(i-1)*d+b,big(2+(triangle==3 ? 2a*b : 0))))
        end
    end
    for i in 1:n, j in i+1:n, a in 1:d, b in 1:d
        value = sum(square_pattern(k,a,b,i,j,d,p) for k in selection;init=0)
        value==0 || push!(quadratic,((i-1)*d+a,(j-1)*d+b,big(value)))
    end
    return QUBOComponent{BigInt}(bits;linear,quadratic,offset,codebooks=books,
        provenance="Q2 candidate atoms; square=$(selection); triangle=$(triangle); n=$(n); d=$(d); p=$(p)")
end

function selections()
    result = Vector{Int}[]
    for mask in 0:(1<<14)-1
        selected = [i for i in 1:14 if !iszero(mask & (1<<(i-1)))]
        count(i->i<=6,selected)<=1 || continue
        any(a in selected && b in selected for (a,b) in ((7,8),(9,10),(12,13))) && continue
        push!(result,selected)
    end
    return result
end

function pattern_truth(family,x,p)
    family===:all_different && return allunique(x)
    family===:ordered && return issorted(x)
    family===:linear_sum && return sum(x)==p
    family===:no_overlap && return all(abs(x[i]-x[j])>=p for i in eachindex(x) for j in i+1:length(x))
    family===:channel && return all(i->x[x[i]]==i,eachindex(x))
    error("unknown pattern calibration family")
end

function pattern_campaign(family)
    n,d,p = family===:no_overlap ? (3,7,2) : family===:linear_sum ? (4,4,10) : (4,4,0)
    candidates = collect(Iterators.product(ntuple(_->1:d,n)...))[:]
    labels = [pattern_truth(family,x,p) for x in candidates]
    first_positive = findfirst(labels)
    first_positive===nothing && error("no positive calibration example")
    features = [sum(square_pattern(k,x[i],x[j],i,j,d,p) for i in 1:n for j in i+1:n;init=0)
        for x in candidates, k in 1:14]
    triangles = [t==1 ? 0 : t==2 ? -n : sum(a*a-2p*a for a in x) for x in candidates,t in 1:3]
    output = joinpath(@__DIR__,"results","patterns",string(family))
    mkpath(output)
    compositions = Dict{String,Any}[]
    energies = zeros(Int,length(candidates))
    for selected in selections(), triangle in 1:3
        for i in eachindex(candidates)
            energies[i] = triangles[i,triangle]+sum(features[i,k] for k in selected;init=0)
        end
        ground = energies[first_positive]
        all(labels[i] ? energies[i]==ground : energies[i]>ground for i in eachindex(labels)) || continue
        raw = raw_pattern_component(selected,triangle,n,d,p;offset=-ground)
        guarded = guard_one_hot(raw)
        q = guarded.component
        # Independent codebook-wide check after canonical construction and weighting.
        for (i,x) in enumerate(candidates)
            z = [a==x[j] for j in 1:n for a in 1:d]
            energy(q,z)==energies[i]-ground || error("pattern reconstruction mismatch")
        end
        filename = "component-$(lpad(length(compositions)+1,4,'0')).toml"
        open(joinpath(output,filename),"w") do io
            write_component(io,q)
        end
        push!(compositions,Dict("square"=>selected,"triangle"=>triangle,"file"=>filename,
            "validity_weight"=>string(guarded.weight),"global_lower_bound"=>string(guarded.lower_bound),
            "codebook_states"=>length(candidates),"global_exactness"=>"finite exhaustive codebook plus analytic validity bound",
            "minimum_positive_gap"=>minimum(energies[i]-ground for i in eachindex(labels) if !labels[i])))
    end
    isempty(compositions) && error("no exact pattern composition for $(family)")
    open(joinpath(output,"campaign.toml"),"w") do io
        TOML.print(io,Dict("family"=>string(family),"n"=>n,"d"=>d,"p"=>p,
            "grammar_candidates"=>3length(selections()),"compositions"=>compositions,
            "scope"=>"complete finite grammar and complete finite codebook; validity corrected analytically; no parametric generalization claim");sorted=true)
    end
    println(family,": ",length(compositions)," globally exact finite-instance candidates archived")
end

if abspath(PROGRAM_FILE)==@__FILE__
    for family in (:all_different,:ordered,:linear_sum,:no_overlap,:channel)
        pattern_campaign(family)
    end
end
