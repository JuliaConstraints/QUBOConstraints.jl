@testset "Linear Regression" begin

    X₃₃ = [rand(0:2, 3) for _ in 1:100]
    B₉ = [rand(Bool, 9) for _ in 1:100]

    training_configs = [
        Dict(
            :info => "No binarization on ⟦0,2⟧³",
            :train => X₃₃,
            :encoding => :none,
            :binarization => :none,
        ),
        # Dict(
        #     :info => "Domain Wall binarization on ⟦0,2⟧³",
        #     :train => X₃₃,
        #     :encoding => :none,
        #     :binarization => :domain_wall,
        # ),
        # Dict(
        #     :info => "One-Hot binarization on ⟦0,2⟧³",
        #     :train => X₃₃,
        #     :encoding => :none,
        #     :binarization => :one_hot,
        # ),
        # Dict(
        #     :info => "Domain Wall pre-encoded on ⟦0,2⟧³",
        #     :train => B₉,
        #     :encoding => :domain_wall,
        #     :binarization => :none,
        # ),
        # Dict(
        #     :info => "One-Hot pre-encoded on ⟦0,2⟧³",
        #     :train => B₉,
        #     :encoding => :one_hot,
        #     :binarization => :none,
        # ),
    ]

    function all_different(x, encoding)
        encoding == :none && (return allunique(x))
        if isvalid(x; encoding)
            b = all_different(debinarize(x; binarization = encoding), :none)
            return b ? -1. : 1.
        else
            return sqrt(length(x))
        end
    end

    function all_different(x, encoding, binarization)
        return all_different(x, encoding == :none ? binarization : encoding)
    end

    for config in training_configs
        @testset "$(config[:info])" begin
            println("\nTest for $(config[:info])")
            penalty = x -> all_different(x, config[:encoding], config[:binarization])
            train(config[:train], penalty; binarization = config[:binarization])
        end
    end

    # @testset "All Different 3×3" begin
    #     println("\nTest for All Different 3×3")
    #     all_different(x) = allunique(x)

    #     binary_all_different(x) = all_different(integerize(x))

    #     opt_all_different(x) = 0.

    #     function sat_all_different(X)
    #         if count(x -> x == sqrt(length(X)), integerize(X)) > 0
    #             return sqrt(length(X))
    #         elseif binary_all_different(X)
    #             return -1.
    #         else
    #             return 1.
    #         end
    #     end

    #     for ss in set_sizes, ds in domain_sizes, st in set_types
    #         str = "ALL_DIFFERENT $ds×$ds: set size = $ss, set type = $(st[1])"
    #         println("\n\t$str\n")
    #         X_train = st[2](ds, ss)
    #         X_check = first(set_types)[2](ds, last(set_sizes))
    #         train(X_train, sat_all_different, opt_all_different; X_check)
    #     end
    # end

    # @testset "Ordered" begin
    #     println("\nTest for Ordered 3×3")
    #     ordered(x) = issorted(x)

    #     binary_ordered(x) = ordered(integerize(x))

    #     opt_ordered(x) = 0.

    #     function sat_ordered(X)
    #         if count(x -> x == sqrt(length(X)), integerize(X)) > 0
    #             return sqrt(length(X))
    #         elseif binary_ordered(X)
    #             return -1.
    #         else
    #             return 1.
    #         end
    #     end

    #     for ss in set_sizes, ds in domain_sizes, st in set_types
    #         str = "ORDERED $ds×$ds: set size = $ss, set type = $(st[1])"
    #         println("\n\t$str\n")
    #         X_train = st[2](ds, ss)
    #         X_check = first(set_types)[2](ds, last(set_sizes))
    #         train(X_train, sat_ordered, opt_ordered; X_check)
    #     end
    # end

    # @testset "LinearSum" begin
    #     function param_set(ds)
    #         N = ds^2
    #         return [
    #             0,
    #             ds,
    #             round(Int,1.5*ds),
    #             round(Int, .25*N),
    #             round(Int, .5*N),
    #             round(Int, .75*N),
    #             N,
    #         ]
    #     end

    #     for ss in set_sizes, ds in domain_sizes, st in set_types, σ in param_set(ds)
    #         linear_sum(x) = sum(x) == σ

    #         binary_linear_sum(x) = linear_sum(integerize(x))

    #         opt_linear_sum(x) = 0.

    #         function sat_linear_sum(X)
    #             if count(x -> x == sqrt(length(X)), integerize(X)) > 0
    #                 return sqrt(length(X))
    #             elseif binary_linear_sum(X)
    #                 return -1.
    #             else
    #                 return abs(sum(integerize(X)) - σ)
    #             end
    #         end

    #         QP = QUBO_linear_sum(ds, σ)

    #         qubo_handmade(X) = transpose(X) * QP * X + σ^2

    #         str = "LINEAR_SUM(=$σ) $ds×$ds: set size = $ss, set type = $(st[1])"
    #         println("\n\t$str\n")
    #         X_train = st[2](ds, ss)
    #         X_check = first(set_types)[2](ds, last(set_sizes))
    #         train(X_train, qubo_handmade, opt_linear_sum; X_check)
    #     end
    # end
end
