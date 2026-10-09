using QUBOConstraints, TOML, SparseArrays

function valid_rows(book)
    width = length(book.bits)
    rows,labels,masks = Vector{Bool}[],Int[],Int[]
    for mask in 0:(1<<width)-1
        z = [!iszero(mask & (1<<(i-1))) for i in 1:width]
        d = decode_code(book,z)
        if d.valid
            push!(rows,z); push!(labels,Int(d.value)+1); push!(masks,mask)
        end
    end
    return (;rows,labels,masks)
end

# Exact integer row operations; gcd normalization controls coefficient growth.
function exact_rank(input)
    a = BigInt.(input)
    row = 1
    for col in axes(a,2)
        pivot = findfirst(i->!iszero(a[i,col]),row:size(a,1))
        pivot===nothing && continue
        p = row+pivot-1
        a[row,:],a[p,:] = copy(a[p,:]),copy(a[row,:])
        for i in row+1:size(a,1)
            x = a[i,col]
            iszero(x) && continue
            for j in col+1:size(a,2)
                a[i,j] = a[row,col]*a[i,j]-x*a[row,j]
            end
            a[i,col] = 0
            divisor = foldl(gcd,view(a,i,col+1:size(a,2));init=big(0))
            if divisor>1
                for j in col+1:size(a,2)
                    a[i,j] = div(a[i,j],divisor)
                end
            end
        end
        row += 1
        row>size(a,1) && break
    end
    return row-1
end

function capabilities(rows,labels,nvalues)
    width = length(first(rows))
    design = reduce(vcat,transpose.([[1;Int.(z);[Int(z[i]&&z[j]) for i in 1:width for j in i+1:width]] for z in rows]))
    targets = [Int(labels[i]==j) for i in eachindex(rows), j in 1:nvalues]
    rank = exact_rank(design)
    augmented = exact_rank(hcat(design,targets))
    return Dict("quadratic_rank"=>rank,"target_augmented_rank"=>augmented,
        "all_semantic_functions_on_valid_codes"=>augmented==rank,
        "scope"=>"codebook only, no auxiliary; global invalid-code protection is a separate condition")
end

function graph_metrics(masks,width)
    lookup = Dict(m=>i for (i,m) in enumerate(masks))
    adjacency = [[lookup[m ⊻ (1<<b)] for b in 0:width-1 if haskey(lookup,m ⊻ (1<<b))] for m in masks]
    components,diameter = 0,0
    seen = falses(length(masks))
    for start in eachindex(masks)
        !seen[start] && (components+=1)
        distances = fill(-1,length(masks)); distances[start]=0
        queue = [start]
        head = 1
        while head<=length(queue)
            node = queue[head]; head+=1; seen[node]=true
            for neighbour in adjacency[node]
                if distances[neighbour]<0
                    distances[neighbour]=distances[node]+1
                    push!(queue,neighbour)
                end
            end
        end
        diameter = max(diameter,maximum(distances))
    end
    return Dict("vertices"=>length(masks),"edges"=>sum(length,adjacency)÷2,
        "components"=>components,"largest_component_diameter"=>diameter,
        "global_diameter"=>components==1 ? string(diameter) : "disconnected")
end

function campaign(output=joinpath(@__DIR__,"results"))
    mkpath(output)
    records = Dict{String,Any}[]
    for entry in encoding_registry().entries
        entry.status===:implemented && entry.encoding!==:signed_offset || continue
        mode = entry.encoding
        for K in 0:(mode===:native ? 1 : 8)
            book = structured_codebook(:x,0:K;encoding=mode,coefficient_bound=2)
            data = valid_rows(book)
            counts = [count(==(i),data.labels) for i in 1:K+1]
            all(>(0),counts) || error("not surjective")
            q = encoding_validity(book)
            report = exhaustive_check(q,v->true;oracle_id="q3/$(mode)/K$(K)")
            report.status===:pass || error("validity counterexample")
            filename = "$(mode)-$(K).toml"
            open(joinpath(output,filename),"w") do io
                write_component(io,q;report)
            end
            canonical = [encode_code(book,v) for v in 0:K]
            hamming = [count(j->canonical[i][j]!=canonical[i+1][j],eachindex(book.bits)) for i in 1:K]
            coefficients = vcat(nonzeros(q.linear),nonzeros(q.quadratic))
            record = Dict{String,Any}("encoding"=>string(mode),"id"=>entry.id,"K"=>K,
                "primary_bits"=>length(book.bits),"validity_auxiliaries"=>length(q.bits)-length(book.bits),
                "injective_on_valid_codes"=>all(==(1),counts),"surjective"=>true,
                "multiplicities"=>counts,"mean_multiplicity"=>sum(counts)/length(counts),
                "valid_fraction"=>length(data.rows)/2^length(book.bits),
                "canonical_neighbour_hamming"=>hamming,"graph"=>graph_metrics(data.masks,length(book.bits)),
                "unary"=>capabilities(data.rows,data.labels,K+1),
                "linear_terms"=>nnz(q.linear),"quadratic_terms"=>nnz(q.quadratic),
                "max_abs_coefficient"=>string(maximum(abs,coefficients;init=big(0))),
                "validity_gap"=>report.min_violation===nothing ? "no_invalid_codes" : string(report.min_violation),
                "rank_decoding_polynomial_degree"=>mode===:gray ? length(book.bits) : Int(any(!iszero,book.weights)),
                "rank_decoding_polynomial_terms"=>mode===:gray ? 2^length(book.bits)-1 : count(!iszero,book.weights),
                "artifact"=>filename,"proof"=>string(report.proof),"evaluations"=>report.evaluations)
            if K<=3
                rows = [vcat(a,b) for a in data.rows for b in data.rows]
                labels = [(a-1)*(K+1)+b for a in data.labels for b in data.labels]
                record["binary"] = capabilities(rows,labels,(K+1)^2)
            else
                record["binary"] = Dict("status"=>"not_checked","reason"=>"pair-codebook rank campaign bounded to K<=3")
            end
            push!(records,record)
        end
    end
    open(joinpath(output,"campaign.toml"),"w") do io
        TOML.print(io,Dict("schema_version"=>"qubo-encoding-campaign/1","julia"=>string(VERSION),
            "scope"=>"E00 K0:1; E01:E10 K0:8; E09 mu2; signed offset covered in tests; explicit-code small graphs",
            "instances"=>records);sorted=true)
    end
    println("Q3 campaign: ",length(records)," exact instances")
end

if abspath(PROGRAM_FILE)==@__FILE__
    campaign(isempty(ARGS) ? joinpath(@__DIR__,"results") : joinpath(@__DIR__,"results",ARGS[1]))
end
