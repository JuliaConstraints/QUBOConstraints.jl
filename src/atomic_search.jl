"""Finite discrete DAG model. No constraint name, recipe, target or stored plan.

Each of `node_slots` slots has five integer decisions: operation (0 disables the
slot), three predecessor ports and a literal. The last decision selects the root.
Source ports are fixed by the input codebooks. Constants and topology are learned
decisions, not data hidden in the decoder. Capacity is explicit, not universal.
Arrays are owned by the model and must be treated as read-only.
"""
struct AtomicSearchSpace
    codebooks::Vector{AbstractCodebook}
    node_slots::Int
    constant_bound::BigInt
    operators::Vector{Symbol}
    max_nodes::Int
    max_local_rows::Int
    max_total_rows::Int
    max_value_bits::Int
    max_depth::Int
end

const _SEARCH_ATOMIC_OPERATORS = (:constant, _ATOMIC_UNARY..., _ATOMIC_BINARY..., Symbol("if"))
_search_arity(op) = op===:constant ? 0 : op in _ATOMIC_UNARY ? 1 : op in _ATOMIC_BINARY ? 2 : 3

function AtomicSearchSpace(books; node_slots::Integer, constant_bound::Integer=65535,
        operators=_SEARCH_ATOMIC_OPERATORS, max_nodes::Integer=100000,
        max_local_rows::Integer=100000, max_total_rows::Integer=1000000,
        max_value_bits::Integer=256, max_depth::Integer=1000)
    0<=node_slots<=div(typemax(Int)-1,5) || throw(ArgumentError("invalid slot budget"))
    constant_bound>=0 || throw(ArgumentError("negative literal bound"))
    all(>(0),(max_nodes,max_local_rows,max_total_rows,max_value_bits,max_depth)) ||
        throw(ArgumentError("positive decoding budgets required"))
    ops=Symbol[operators...]
    !isempty(ops) && allunique(ops) && all(in(_SEARCH_ATOMIC_OPERATORS),ops) ||
        throw(ArgumentError("invalid operation alphabet"))
    owned=AbstractCodebook[deepcopy(b) for b in books]
    allunique(b.variable for b in owned) || throw(ArgumentError("duplicate source identities"))
    all(b->all(v->v isa Integer,b.values),owned) || throw(ArgumentError("integer sources required"))
    for book in owned, value in book.values
        width=ndigits(abs(big(value));base=2)
        width<=max_value_bits || throw(CompilationLimit(:search_source_bits,big(width),big(max_value_bits)))
    end
    length(owned)+big(node_slots)>0 || throw(ArgumentError("model has no possible root"))
    AtomicSearchSpace(owned,Int(node_slots),big(constant_bound),ops,Int(max_nodes),
        Int(max_local_rows),Int(max_total_rows),Int(max_value_bits),Int(max_depth))
end

"""Independent integer bounds; dependent port/arity rules are checked by decoding."""
function atomic_decision_domains(space::AtomicSearchSpace)
    domains=Tuple{BigInt,BigInt}[]; n=length(space.codebooks)
    for slot in 1:space.node_slots
        push!(domains,(big(0),big(length(space.operators))))
        append!(domains,fill((big(0),big(n+slot-1)),3))
        push!(domains,(-space.constant_bound,space.constant_bound))
    end
    push!(domains,(big(1),big(n+space.node_slots)))
    domains
end

function atomic_operator_code(space::AtomicSearchSpace, op::Symbol)
    code=findfirst(==(op),space.operators)
    code===nothing && throw(ArgumentError("operation $op not in model alphabet"))
    code
end

"""Encode a manually obtained witness into the SAME decisions exposed to train.

The plan is consumed only here; the returned vector contains all graph choices.
No model field is changed, and no lookup table of successful recipes is installed.
"""
function encode_atomic_weights(space::AtomicSearchSpace, plan::AtomicPlan)
    validate_atomic_plan(plan;max_rows=space.max_total_rows,max_value_bits=space.max_value_bits)
    length(space.codebooks)==length(plan.codebooks) &&
        all(_compatible(a,b) for (a,b) in zip(space.codebooks,plan.codebooks)) ||
        throw(ArgumentError("source codebooks differ from search model"))
    weights=zeros(BigInt,5space.node_slots+1)
    ports=zeros(Int,length(plan.nodes)); slot=0; n=length(space.codebooks)
    for (i,node) in enumerate(plan.nodes)
        if node.operator===:variable
            ports[i]=node.data
            continue
        end
        slot+=1
        slot<=space.node_slots || throw(CompilationLimit(:search_slots,big(slot),big(space.node_slots)))
        at=5(slot-1); weights[at+1]=atomic_operator_code(space,node.operator)
        for (j,input) in enumerate(node.inputs)
            ports[input]>0 || throw(ArgumentError("non-topological witness"))
            weights[at+1+j]=ports[input]
        end
        if node.operator===:constant
            abs(node.data)<=space.constant_bound || throw(CompilationLimit(:search_literal,abs(node.data),space.constant_bound))
            weights[at+5]=node.data
        end
        ports[i]=n+slot
    end
    weights[end]=ports[plan.root]
    weights
end

"""Decode only model + weights. No core recipe or semantic oracle is consulted.

Invalid discrete programs are rejected, never repaired with a constraint-specific
answer. Resource refusal is distinct from a proof of nonrepresentability.
"""
function decode_atomic_weights(space::AtomicSearchSpace, weights::AbstractVector{<:Integer})
    length(weights)==5space.node_slots+1 || throw(DimensionMismatch("wrong decision count"))
    n=length(space.codebooks); available=zeros(Int,n+space.node_slots)
    nodes=AtomicNode[]; heights=Int[]; total=big(0)
    n<=space.max_nodes || throw(CompilationLimit(:atomic_nodes,big(n),big(space.max_nodes)))
    for (i,b) in enumerate(space.codebooks)
        push!(nodes,AtomicNode(:variable,Int[],i,BigInt.(b.values),Vector{Int}[]))
        push!(heights,1); available[i]=i
    end
    for slot in 1:space.node_slots
        at=5(slot-1); code=weights[at+1]
        0<=code<=length(space.operators) || throw(ArgumentError("operation decision out of bounds"))
        literal=weights[at+5]
        abs(big(literal))<=space.constant_bound || throw(ArgumentError("literal decision out of bounds"))
        if code==0
            all(iszero,weights[at+2:at+5]) || throw(ArgumentError("disabled slot must have zero payload"))
            continue
        end
        op=space.operators[Int(code)]; arity=_search_arity(op)
        inputs=Int[]
        for port in 1:3
            predecessor=weights[at+1+port]
            if port<=arity
                1<=predecessor<n+slot || throw(ArgumentError("port must select a predecessor"))
                available[Int(predecessor)]>0 || throw(ArgumentError("port selects a disabled slot"))
                push!(inputs,available[Int(predecessor)])
            else
                iszero(predecessor) || throw(ArgumentError("unused port must be zero"))
            end
        end
        op===:constant || iszero(literal) || throw(ArgumentError("nonconstant literal must be zero"))
        if op===:constant
            width=ndigits(abs(big(literal));base=2)
            width<=space.max_value_bits || throw(CompilationLimit(:search_literal_bits,big(width),big(space.max_value_bits)))
        end
        length(nodes)<space.max_nodes || throw(CompilationLimit(:atomic_nodes,big(length(nodes)+1),big(space.max_nodes)))
        height=isempty(inputs) ? 1 : 1+maximum(heights[i] for i in inputs)
        height<=space.max_depth || throw(CompilationLimit(:atomic_depth,big(height),big(space.max_depth)))
        rows=Vector{Int}[]; values=BigInt[]
        if op===:constant
            push!(values,big(literal))
        else
            count=prod((big(length(nodes[i].domain)) for i in inputs);init=big(1))
            count<=space.max_local_rows || throw(CompilationLimit(:atomic_local_rows,count,big(space.max_local_rows)))
            total+=count
            total<=space.max_total_rows || throw(CompilationLimit(:atomic_total_rows,total,big(space.max_total_rows)))
            for indices in Iterators.product((eachindex(nodes[i].domain) for i in inputs)...)
                arguments=BigInt[nodes[i].domain[k] for (i,k) in zip(inputs,indices)]
                value=_atomic_apply(op,arguments,space.max_value_bits)
                k=findfirst(==(value),values)
                if k===nothing
                    push!(values,value); k=length(values)
                end
                push!(rows,[Int[indices...];k])
            end
        end
        push!(nodes,AtomicNode(op,inputs,op===:constant ? big(literal) : nothing,values,rows))
        push!(heights,height); available[n+slot]=length(nodes)
    end
    root=weights[end]
    1<=root<=length(available) && available[Int(root)]>0 || throw(ArgumentError("invalid root decision"))
    all(in((0,1)),nodes[available[Int(root)]].domain) || throw(ArgumentError("root must be Boolean"))
    plan=AtomicPlan(deepcopy(space.codebooks),nodes,available[Int(root)],"atomic-search/integer-DAG-v1")
    validate_atomic_plan(plan;max_rows=space.max_total_rows,max_value_bits=space.max_value_bits)
    plan
end

compile_atomic_weights(space::AtomicSearchSpace, weights;kwargs...) =
    compile_atomic(decode_atomic_weights(space,weights);kwargs...)

"""Reference exhaustive discrete learner, NOT an efficient large-campaign backend.

It visits the public integer box without dropping any weight assignment, rejecting
invalid programs. A finite budget never implies impossibility. No recipe seed is
required. More efficient optimizers can use exactly the same domains/decoder.
"""
struct AtomicEnumerativeOptimizer <: AbstractOptimizer
    max_candidates::Int
    function AtomicEnumerativeOptimizer(;max_candidates::Integer=10000)
        max_candidates>0 || throw(ArgumentError("positive search budget required"))
        new(Int(max_candidates))
    end
end

struct AtomicTrainingResult
    weights::Union{Nothing,Vector{BigInt}}
    mismatches::Union{Nothing,Int}
    candidates::Int
    rejected::Int
    resource_limited::Int
    status::Symbol
end

function train(space::AtomicSearchSpace, assignments::AbstractVector, labels::AbstractVector{Bool};
        optimizer::AtomicEnumerativeOptimizer=AtomicEnumerativeOptimizer(), initial_weights=nothing)
    !isempty(assignments) && length(assignments)==length(labels) || throw(ArgumentError("nonempty paired training data required"))
    # Validate data OUTSIDE the candidate loop: a bad corpus is not a bad model.
    for assignment in assignments, book in space.codebooks
        haskey(assignment,book.variable) && assignment[book.variable] in book.values ||
            throw(ArgumentError("training assignment outside model sources"))
    end
    domains=atomic_decision_domains(space); origin=BigInt[first(d) for d in domains]
    initial_weights===nothing || initial_weights isa AbstractVector{<:Integer} ||
        throw(ArgumentError("initial decisions must be an integer vector"))
    weights=initial_weights===nothing ? copy(origin) : BigInt.(initial_weights)
    length(weights)==length(domains) && all(lo<=w<=hi for (w,(lo,hi)) in zip(weights,domains)) ||
        throw(ArgumentError("initial decisions outside model domains"))
    complete_box=weights==origin
    best=nothing; loss=nothing; rejected=0; limited=0
    for tried in 1:optimizer.max_candidates
        plan=try
            decode_atomic_weights(space,weights)
        catch error
            if error isa CompilationLimit
                limited+=1
            elseif error isa ArgumentError || error isa DimensionMismatch
                rejected+=1
            else
                rethrow()
            end
            nothing
        end
        if plan!==nothing
            misses=count(((d,label),)->(atomic_values(plan,d)[plan.root]==1)!=label,zip(assignments,labels))
            if loss===nothing || misses<loss
                best=copy(weights); loss=misses
            end
            misses==0 && return AtomicTrainingResult(best,loss,tried,rejected,limited,:training_fit)
        end
        # Mixed-radix successor over BigInt bounds, including every legal witness.
        digit=1
        while digit<=length(weights) && weights[digit]==last(domains[digit])
            weights[digit]=first(domains[digit]); digit+=1
        end
        digit>length(weights) && return AtomicTrainingResult(best,loss,tried,rejected,limited,
            !complete_box ? :suffix_exhausted : limited==0 ? :space_exhausted : :space_exhausted_under_resource_limits)
        weights[digit]+=1
    end
    AtomicTrainingResult(best,loss,optimizer.max_candidates,rejected,limited,:budget_exhausted)
end
