@testset "Linear Regression" begin
    X₃₃ = [rand(0:2, 3) for _ in 1:10]
    X₃₃_test = [rand(0:2, 3) for _ in 1:100]

    B₉ = [rand(Bool, 9) for _ in 1:10]
    B₉_test = [rand(Bool, 9) for _ in 1:2000]

    training_configs = [
        Dict(
            :info => "No binarization on ⟦0,2⟧³",
            :train => X₃₃,
            :test => X₃₃_test,
            :encoding => :none,
            :binarization => :none
        ),
        Dict(
            :info => "Domain Wall binarization on ⟦0,2⟧³",
            :train => X₃₃,
            :test => X₃₃_test,
            :encoding => :none,
            :binarization => :domain_wall
        ),
        Dict(
            :info => "One-Hot binarization on ⟦0,2⟧³",
            :train => X₃₃,
            :test => X₃₃_test,
            :encoding => :none,
            :binarization => :one_hot
        ),
        Dict(
            :info => "Domain Wall pre-encoded on ⟦0,2⟧³",
            :train => B₉,
            :test => B₉_test,
            :encoding => :domain_wall,
            :binarization => :none
        ),
        Dict(
            :info => "One-Hot pre-encoded on ⟦0,2⟧³",
            :train => B₉,
            :test => B₉_test,
            :encoding => :one_hot,
            :binarization => :none
        )
    ]

    function all_different(x, encoding)
        encoding == :none && (return allunique(x))
        isv = if encoding == :one_hot
            mapreduce(i -> is_valid(x[i:(i + 2)], Val(encoding)), *, 1:3:9)
        else
            mapreduce(i -> is_valid(x[i:(i + 1)], Val(encoding)), *, 1:2:6)
        end
        if isv
            b = all_different(debinarize(x; binarization = encoding), :none)
            return b ? -1.0 : 1.0
        else
            return length(x)
        end
    end

    function all_different(x, encoding, binarization)
        return all_different(x, encoding == :none ? binarization : encoding)
    end

    for config in training_configs
        @testset "$(config[:info])" begin
            println("\nTest for $(config[:info])")
            penalty = x -> all_different(x, config[:encoding], config[:binarization])
            optimizer = GradientDescentOptimizer(;
                binarization = config[:binarization]
            )
            train(config[:train], penalty; optimizer, X_test = config[:test])
        end
    end
end
