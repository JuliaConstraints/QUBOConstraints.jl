function QUBO_base(n, weight = 1)
    N = n^2
    M = zeros(N, N)
    for k in 1:n, j in ((k-1)*n+1):(k*n), i in ((k-1)*n+1):(k*n)
        M[j,i] = weight
    end
    return M - 2 * weight * I
end
