"""
    QUBO_linear_sum(n, σ)

Exact upper-triangular BigInt matrix for n one-hot variables with domains 0:n-1,
including one-hot validity. Add the constant offset n+σ² for a zero-set penalty.

The old Float64 matrix omitted validity and is no longer returned. Bit order remains
variable-major then value-major. Prefer `one_hot_constraint(:linear_sum, domains; rhs=σ)`
for explicit domains, offset, metadata and certification.
"""
function QUBO_linear_sum(n::Integer, σ::Integer)
    1<=n<=isqrt(typemax(Int)) || throw(ArgumentError("invalid dense variable count"))
    width = Int(n)
    N = width^2
    Q = fill(big(0),N,N)
    weights = BigInt[mod(i-1,width) for i in 1:N]
    target = big(σ)
    for j in 1:N
        Q[j,j] = weights[j]^2 - 2target*weights[j] - 1
        for i in 1:j-1
            Q[i,j] = 2weights[i]*weights[j] + (div(i-1,width)==div(j-1,width) ? 2 : 0)
        end
    end
    return Q
end
