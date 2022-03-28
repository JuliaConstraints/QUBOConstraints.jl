@testset "Linear Regression" begin
    @testset "All Different 3×3" begin
        println("\nTest for All Different 3×3")
        X_train = map(collect ∘ binarize ∘ collect, Iterators.product(0:2,0:2,0:2))

        all_different(x) = allunique(x)

        binary_all_different(x) = all_different(integerize(x))

        opt_all_different(x) = 0.

        function sat_all_different(X)
            if count(x -> x == sqrt(length(X)), integerize(X)) > 0
                return sqrt(length(X))
            elseif binary_all_different(X)
                return -1.
            else
                return 1.
            end
        end

        Q = zeros(9,9)

        @show train!(Q, X_train, sat_all_different, opt_all_different)

        display(Q)
    end

    @testset "Ordered" begin
        println("\nTest for Ordered 3×3")
        X_train = map(collect ∘ binarize ∘ collect, Iterators.product(0:2,0:2,0:2))

        ordered(x) = issorted(x)

        binary_ordered(x) = ordered(integerize(x))

        opt_ordered(x) = 0.

        function sat_ordered(X)
            if count(x -> x == sqrt(length(X)), integerize(X)) > 0
                return sqrt(length(X))
            elseif binary_ordered(X)
                return -1.
            else
                return 1.
            end
        end

        Q = zeros(9,9)

        @show train!(Q, X_train, sat_ordered, opt_ordered)

        display(Q)
    end

    @testset "LinearSum" begin
        println("\nTest for Linear Sum 3×3")
        X_train = map(collect ∘ binarize ∘ collect, Iterators.product(0:2,0:2,0:2))

        σ = 4

        linear_sum(x) = sum(x) == σ

        binary_linear_sum(x) = linear_sum(integerize(x))

        opt_linear_sum(x) = 0.

        function sat_linear_sum(X)
            if count(x -> x == sqrt(length(X)), integerize(X)) > 0
                return sqrt(length(X))
            elseif binary_linear_sum(X)
                return -1.
            else
                return 1. + abs(sum(X) - σ)
            end
        end

        Q = zeros(9,9)

        @show train!(Q, X_train, sat_linear_sum, opt_linear_sum)

        display(Q)
    end
end
