predict(x, Q) = transpose(x) * Q * x

loss(x, y, Q) = (predict(x, Q) .-y).^2

function make_df(X, Q, penalty)
	df = DataFrame()
	for (i,x) in enumerate(X)
		if i == 1
			df = DataFrame(transpose(x), :auto)
		else
			push!(df, transpose(x))
		end
	end

	dim = length(df[1,:])

	df[!,:penalty] = map(r -> penalty(Vector(r)), eachrow(df))
	df[!,:predict] = map(r -> predict(Vector(r), Q), eachrow(df[:, 1:dim]))

	min_false = minimum(
        filter(:penalty => >(minimum(df[:,:penalty])), df)[:,:predict];
        init = typemax(Int)
    )
    df[!,:shifted] = df[:,:predict] .- min_false
    df[!,:accurate] = df[:, :penalty] .* df[:,:shifted] .≥ 0.

	return df
end

function oversample(X, f)
    X_true = Vector{eltype(X)}()
    X_false = Vector{eltype(X)}()

    μ = minimum(f, X)

    foreach(x -> push!(f(x) == μ ? X_true : X_false, x), X)

    b = length(X_true) > length(X_false)
    Y = reverse(b ? X_true : X_false)
    it = Iterators.cycle(b ? X_false : X_true)

    Z = Vector{eltype(X)}()
    l = length(Y)
    for (i, x) in enumerate(it)
        push!(Z, x, Y[i])
        i == l && break
    end

    return Z
end

function preliminaries(X, domains, binarization)
    Y = binarization==:none ? X : map(x -> collect(binarize(X, domains; binarization)), X)
    n = length(first(Y))
    return Y, zeros(n,n)
end

function preliminaries(X, _)
    n = length(first(X))
    return X, zeros(n,n)
end

function train!(Q, X, penalty, η, precision, X_check, oversampling)
    θ = params(Q)
    for x in (oversampling ? oversample(X, penalty) : X)
        grads = gradient(() -> loss(x, penalty(x), Q), θ)
        Q .-= η * grads[Q]
    end

    Q[:,:] = round.(precision*Q)

    df = make_df(X_check, Q, penalty)
    return pretty_table(describe(df[!, [:penalty, :predict, :shifted, :accurate]]))
end

function train(
    X,
    domains::Vector{D},
    penalty;
    η = .001,
    precision = 5,
    X_check = X,
    binarization = :none,
    oversampling = false,
) where {D <: DiscreteDomain}
    Y, Q = preliminaries(X, domains, binarization)
    return train!(Q, Y, penalty, η, precision, X_check, oversampling)
end

function train(
    X, penalty;
    η = .001, precision = 5, X_check = X, binarization = :one_hot, oversampling = false,
)
    return train(
        X, to_domains(X), penalty;
        η, precision, X_check, binarization, oversampling,
    )
end

function train(
    X, domain_size::Int, penalty;
    η = .001, precision = 5, X_check = X, binarization = :one_hot, oversampling = false,
)
    return train(
        X, to_domains(X, domain_size), penalty;
        η, precision, X_check, binarization, oversampling
    )
end

function train(
    X, domain_sizes::Vector{Int}, penalty;
    η = .001, precision = 5, X_check = X, binarization = :one_hot, oversampling = false,
)
    return train(
        X, to_domains(X, domain_sizes), penalty;
        η, precision, X_check, binarization, oversampling
    )
end

function train(
    X, d::D, penalty;
    η = .001, precision = 5, X_check = X, binarization = :one_hot, oversampling = false,
) where {D <: DiscreteDomain}
    return train(
        X, to_domains(X, d), penalty;
        η, precision, X_check, binarization, oversampling
    )
end
