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

	# max_true = maximum(filter(:sat => >=(0.), df)[:,:predict])
	min_false = minimum(filter(:sat => >(0.), df)[:,:predict])

	# shift = 0.
    # if max_true ≥ min_false
    #     shift = (min_false + max_true) / 2

    # min_false > 0. ? (-min_false) : min_false
    df[!,:shifted] = df[:,:predict] .- min_false

    df[!,:accurate] = df[:, :sat] .* df[:,:shifted] .≥ 0.

	return df
end

function train!(Q, X, sat, opt; η = .001, precision = 100)
    θ = params(Q)
    for x in X
        grads = gradient(() -> loss(x, opt(x) + sat(x), Q), θ)
        update!(Q, η * grads[Q])
    end
    display(Q)
    Q[:,:] = round.(precision*Q) + QUBO_base(Int(sqrt(size(Q,1))))

    return pretty_table(make_df(X, Q, sat, opt))
end
