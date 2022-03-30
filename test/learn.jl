@testset "Linear Regression" begin
    B_train = [(rand(0:1,9)) for _ in 1:100000]
    X_train = [collect(binarize(rand(0:2,3))) for _ in 1:100000]
    @testset "All Different 3×3" begin
        println("\nTest for All Different 3×3")
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

        for train in [X_train, B_train]
            str = train == X_train ? "Only valid configurations" : "One-Hot encoding complete space"
            println("\n\t$str\n")
            Q = zeros(9,9)
            train!(Q, train, sat_all_different, opt_all_different)
            display(Q)
        end
    end

    @testset "Ordered" begin
        println("\nTest for Ordered 3×3")
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

        for train in [X_train, B_train]
            str = train == X_train ? "Only valid configurations" : "One-Hot encoding complete space"
            println("\n\t$str\n")
            Q = zeros(9,9)
            train!(Q, train, sat_ordered, opt_ordered)
            display(Q)
        end
    end

    @testset "LinearSum" begin
        for σ in 0:1:6
            println("\nTest for Linear Sum 3×3, σ = $σ")

            QP = QUBO_linear_sum(3, σ)

            linear_sum(x) = sum(x) == σ

            binary_linear_sum(x) = linear_sum(integerize(x))

            opt_linear_sum(x) = 0.

            function sat_linear_sum(X)
                if count(x -> x == sqrt(length(X)), integerize(X)) > 0
                    return sqrt(length(X))
                elseif binary_linear_sum(X)
                    return -1.
                else
                    return abs(sum(integerize(X)) - σ)
                end
            end

            qubo_handmade(X) = transpose(X) * QP * X + σ^2 - 0.5

            for train in [X_train, B_train]
                str = train == X_train ? "Only valid configurations" : "One-Hot encoding complete space"
                println("\n\t$str\n")
                Q = zeros(9,9)
                train!(Q, train, qubo_handmade, opt_linear_sum)
                display(Q)
            end
        end
    end
end
