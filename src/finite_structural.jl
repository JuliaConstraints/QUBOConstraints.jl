# Explicit finite-cube routes: their exponential bound is part of the public contract.
function _finite_scope(books,scope)
    books=collect(books)
    all(b->b isa AbstractCodebook,books) || throw(ArgumentError("codebooks required"))
    allunique(b.variable for b in books) || throw(ArgumentError("duplicate codebook variable"))
    positions=Dict(b.variable=>i for (i,b) in enumerate(books))
    scope=Symbol.(collect(scope))
    all(s->haskey(positions,s),scope) || throw(ArgumentError("unknown scope variable"))
    all(b->all(v->v isa Integer,_book_values(b)),books) || throw(ArgumentError("integer structural domains required"))
    return books,[positions[s] for s in scope],positions
end

"""Finite exact lexicographic comparison of equally long lists; not a scalable lex encoding."""
function lex_component(books,left,right;operator::Symbol=:le,kwargs...)
    operator in (:lt,:le,:gt,:ge) || throw(ArgumentError("unsupported lexicographic order"))
    length(left)==length(right) || throw(DimensionMismatch("equal-length lex lists required"))
    books,l,_=_finite_scope(books,left)
    _,r,_=_finite_scope(books,right)
    function predicate(v)
        for (i,j) in zip(l,r)
            v[i]==v[j] && continue
            return operator in (:lt,:le) ? v[i]<v[j] : v[i]>v[j]
        end
        return operator in (:le,:ge)
    end
    return truth_component(books,predicate;oracle_id="Q6/lex/finite-v1/$(operator)/$(l)/$(r)",kwargs...)
end

"""Finite exact element with a constant IntegerCondition and optional variable index.

Without index the condition is existential. With index it addresses one scope position.
Nonzero start_index is an explicit general-format variant, not a core coverage claim.
FIRST/LAST ranks and matrix/variable-condition variants are deliberately not implicit.
"""
function element_component(books,scope,condition::IntegerCondition;index::Union{Nothing,Symbol}=nothing,
        start_index::Integer=0,kwargs...)
    books,positions,byname=_finite_scope(books,scope)
    index===nothing && start_index!=0 && throw(ArgumentError("start_index requires indexed element"))
    index===nothing || haskey(byname,index) || throw(ArgumentError("unknown index variable"))
    function predicate(v)
        index===nothing && return any(i->satisfies(condition,v[i]),positions)
        p=big(v[byname[index]])-start_index
        return 0<=p<length(positions) && satisfies(condition,v[positions[Int(p)+1]])
    end
    return truth_component(books,predicate;
        oracle_id="Q6/element/finite-v1/$(positions)/$(index)/$(start_index)/$(condition)",kwargs...)
end

"""Finite exact ordered first occurrences, with optional covered final value."""
function precedence_component(books,scope,values;covered::Bool=false,kwargs...)
    books,positions,_=_finite_scope(books,scope)
    values=collect(values)
    !isempty(values) && allunique(values) && all(v->v isa Integer,values) ||
        throw(ArgumentError("nonempty distinct integer precedence values required"))
    function predicate(v)
        word=v[positions]
        firsts=[findfirst(==(value),word) for value in values]
        covered && last(firsts)===nothing && return false
        return all(j->firsts[j+1]===nothing || (firsts[j]!==nothing && firsts[j]<firsts[j+1]),1:length(values)-1)
    end
    return truth_component(books,predicate;oracle_id="Q6/precedence/finite-v1/$(positions)/$(values)/covered=$(covered)",kwargs...)
end

"""Finite exact run lengths and optional consecutive-run patterns (historical/full stretch).

Each value has one positive integer interval of allowed run lengths. Every word value
must be listed. Empty words satisfy this mathematical extension. Not a core-3.2 claim.
"""
function stretch_component(books,scope,values,widths;patterns=nothing,kwargs...)
    books,positions,_=_finite_scope(books,scope)
    values=collect(values); widths=collect(widths)
    length(values)==length(widths) && allunique(values) && all(v->v isa Integer,values) ||
        throw(ArgumentError("distinct stretch values and matching width intervals required"))
    all(w->w isa AbstractUnitRange && eltype(w)<:Integer && !isempty(w) && first(w)>0,widths) ||
        throw(ArgumentError("positive run-length intervals required"))
    patterns=patterns===nothing ? nothing : collect(patterns)
    patterns===nothing || all(p->length(p)==2 && all(in(values),p),patterns) || throw(ArgumentError("invalid run pattern"))
    function predicate(v)
        word=v[positions]; previous=nothing; i=1
        while i<=length(word)
            k=findfirst(==(word[i]),values)
            k===nothing && return false
            j=i+1
            while j<=length(word) && word[j]==word[i]; j+=1; end
            j-i in widths[k] || return false
            if previous!==nothing && patterns!==nothing
                any(p->p[1]==previous && p[2]==word[i],patterns) || return false
            end
            previous=word[i]; i=j
        end
        return true
    end
    return truth_component(books,predicate;
        oracle_id="Q6/stretch/full-finite-v1/$(positions)/$(values)/$(widths)/$(patterns)",kwargs...)
end
