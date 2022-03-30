predict(x, Q) = transpose(x) * Q * x

loss(x, y, Q) = (predict(x, Q) .-y).^2

function make_df(X, Q, sat, opt)
	df = DataFrame()
	for (i,x) in enumerate(X)
		if i == 1
			df = DataFrame(transpose(x), :auto)
		else
			push!(df, transpose(x))
		end
	end

	dim = length(df[1,:])

	df[!,:sat] = map(r -> sat(Vector(r)), eachrow(df))
	df[!,:opt] = map(r -> opt(Vector(r)), eachrow(df))
	df[!,:predict] = map(r -> predict(Vector(r), Q), eachrow(df[:, 1:dim]))

	min_false = minimum(
        filter(:sat => >(minimum(df[:,:sat])), df)[:,:predict];
        init = typemax(Int)
    )
    df[!,:shifted] = df[:,:predict] .- min_false
    df[!,:accurate] = df[:, :sat] .* df[:,:shifted] .≥ 0.

	return df
end

function oversampling(X, sat, opt)
    X_true = Vector{eltype(X)}()
    X_false = Vector{eltype(X)}()

    f = x -> sat(x) + opt(x)
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

function train!(Q, X, sat, opt; η = .001, precision = 5)
    θ = params(Q)
    # for x in oversampling(X, sat, opt)
    for x in X
        grads = gradient(() -> loss(x, opt(x) + sat(x), Q), θ)
        update!(Q, η * grads[Q])
    end

    Q[:,:] = round.(precision*Q)

    df = make_df(X, Q, sat, opt)
    return pretty_table(describe(df))
end
