@testset "Linear Regression" begin

    set_sizes = [100, 1000, 10000, 100000]

    domain_sizes = [3, 5, 10]

    set_types = [
        # Binary representation (either domain wall or one-hot)
        "binary" => (ds, ss) -> [rand(0:1, ds^2) for _ in 1:ss],
        # Integer representation
        "integer" => (ds, ss) -> [collect(binarize(rand(0:ds-1, ds))) for _ in 1:ss]
    ]

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

        for ss in set_sizes, ds in domain_sizes, st in set_types
            str = "ALL_DIFFERENT $ds×$ds: set size = $ss, set type = $(st[1])"
            println("\n\t$str\n")
            X_train = st[2](ds, ss)
            X_check = first(set_types)[2](ds, last(set_sizes))
            N = ds^2
            Q = zeros(N, N)
            train!(Q, X_train, sat_all_different, opt_all_different; X_check)
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

        for ss in set_sizes, ds in domain_sizes, st in set_types
            str = "ORDERED $ds×$ds: set size = $ss, set type = $(st[1])"
            println("\n\t$str\n")
            X_train = st[2](ds, ss)
            X_check = first(set_types)[2](ds, last(set_sizes))
            N = ds^2
            Q = zeros(N, N)
            train!(Q, X_train, sat_ordered, opt_ordered; X_check)
        end
    end

    @testset "LinearSum" begin
        function param_set(ds)
            N = ds^2
            return [
                0,
                ds,
                round(Int,1.5*ds),
                round(Int, .25*N),
                round(Int, .5*N),
                round(Int, .75*N),
                N,
            ]
        end

        for ss in set_sizes, ds in domain_sizes, st in set_types, σ in param_set(ds)
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

            QP = QUBO_linear_sum(ds, σ)

            qubo_handmade(X) = transpose(X) * QP * X + σ^2

            str = "LINEAR_SUM(=$σ) $ds×$ds: set size = $ss, set type = $(st[1])"
            println("\n\t$str\n")
            X_train = st[2](ds, ss)
            X_check = first(set_types)[2](ds, last(set_sizes))
            N = ds^2
            Q = zeros(N, N)
            train!(Q, X_train, qubo_handmade, opt_linear_sum; X_check)
        end
    end
end
