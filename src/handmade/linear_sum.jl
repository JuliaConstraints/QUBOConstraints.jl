function QUBO_linear_sum(n, σ)
    N = n^2
    Q = zeros(N,N)
    for i in 1:n, k in 1:n
        a = (i - 1) * n + k
        Q[a, a] = -2 * σ * (k - 1)
        for j in 1:n, l in 1:n
            b = (j - 1) * n + l
            Q[a, b] += (k - 1) * (l - 1)
        end
    end
    return Q
end
